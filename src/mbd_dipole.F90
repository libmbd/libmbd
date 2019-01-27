! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MBD_TYPE

module mbd_dipole
!! Construction of dipole tensors and dipole matrices.

use mbd_constants
use mbd_matrix, only: matrix_re_t, matrix_cplx_t
use mbd_geom, only: geom_t, supercell_circum
use mbd_damping, only: damping_t, damping_fermi, damping_sqrtfermi, &
    op1minus_grad
use mbd_gradients, only: grad_t, grad_matrix_re_t, grad_matrix_cplx_t, &
    grad_scalar_t, grad_request_t
use mbd_lapack, only: eigvals, inverse
use mbd_linalg, only: outer
use mbd_utils, only: tostr, shift_idx

implicit none

private
public :: dipole_matrix, T_bare, T_erf_coulomb, damping_grad, T_erfc, B_erfc, C_erfc

interface dipole_matrix
    !! Form either a real or a complex dipole matrix.
    !!
    !! The real-typed version is equivalent to \(\mathbf q=0\).
    !!
    !! $$
    !! \boldsymbol{\mathcal A}:=[\mathbf a_1\mathbf a_2\mathbf
    !! a_3],\qquad\boldsymbol{\mathcal B}:=[\mathbf b_1\mathbf b_2\mathbf b_3]
    !! \\ \boldsymbol{\mathcal B}=2\pi(\boldsymbol{\mathcal A}^{-1})^\mathrm
    !! T,\qquad \partial\boldsymbol{\mathcal B}=-\big((\partial\boldsymbol{\mathcal
    !! A})\boldsymbol{\mathcal A}^{-1}\big)^\mathrm T\boldsymbol{\mathcal B}
    !! \\ \mathbf R_\mathbf n=\boldsymbol{\mathcal A}\mathbf
    !! n,\qquad\partial\mathbf R_\mathbf n=(\partial\boldsymbol{\mathcal
    !! A})\mathbf n,
    !! \\ \mathbf G_\mathbf m=\boldsymbol{\mathcal B}\mathbf m,\qquad
    !! \partial\mathbf G_\mathbf m=-\big((\partial\boldsymbol{\mathcal
    !! A})\boldsymbol{\mathcal A}^{-1}\big)^\mathrm T\mathbf G_\mathbf m,
    !! \\ \frac{\partial G_{\mathbf ma}}{\partial A_{bc}}=-\mathcal A^{-1}_{ca}G_{\mathbf mb}
    !! $$
    !!
    !! $$
    !! \begin{gathered}
    !! \mathbf T_{ij}(\mathbf q)=\sum_{\mathbf n}\mathbf T(\mathbf R_{\mathbf
    !! nij})\mathrm e^{-\mathrm i\mathbf q\cdot\mathbf R_{\mathbf nij}},\quad\mathbf
    !! R_{\mathbf nij}=\mathbf R_j+\mathbf R_\mathbf n-\mathbf R_i
    !! \\ \frac{\mathrm d\mathbf R_{\mathbf nij}}{\mathrm d\mathbf
    !! R_k}=(\delta_{jk}-\delta_{ik})\mathbf I
    !! \\ \mathbf{T}_{ij}(\mathbf{q})\approx\mathbf{T}^\text{Ew}_{ij}(\mathbf{q})
    !! =\sum_\mathbf n^{|\mathbf R_{\mathbf nij}|<R_\text c}\mathbf
    !! T^\text{erfc}(\mathbf R_{\mathbf nij};\gamma)\mathrm e^{-\mathrm i\mathbf
    !! q\cdot\mathbf R_{\mathbf nij}} +\frac{4\pi}{V_\text{uc}}\sum_{\mathbf
    !! m}^{0<|\mathbf k_\mathbf m|<k_\text c}\mathbf{\hat k}_\mathbf
    !! m\otimes\mathbf{\hat k}_\mathbf m\,\mathrm e^{-\frac{k_\mathbf
    !! m^2}{4\gamma^2}-\mathrm i\mathbf G_\mathbf m\cdot\mathbf R_{ij}}
    !! \\ -\frac{4\gamma^3}{3\sqrt\pi}\delta_{ij}\mathbf I +\delta(\mathbf q)\frac{4
    !! \pi}{3 V_\text{uc}}\mathbf I,\qquad \mathbf k_\mathbf m=\mathbf G_\mathbf
    !! m+\mathbf q
    !! \end{gathered}
    !! $$
    !!
    !!
    !! $$
    !! \partial\left(\frac{4\pi}{V_\text{uc}}\right)=-\frac{4\pi}{V_\text{uc}}\frac{\partial
    !! V_\text{uc}}{V_\text{uc}},\qquad\frac{\partial
    !! V_\text{uc}}{\partial\boldsymbol{\mathcal
    !! A}}=V_\text{uc}(\boldsymbol{\mathcal A}^{-1})^\mathrm T
    !! \\ \partial(k^2)=2\mathbf k\cdot\partial\mathbf k
    !! \\ \mathbf{\hat k}\otimes\partial\mathbf{\hat k}=\frac{\mathbf
    !! k\otimes\partial\mathbf k}{k^2}-\frac{\mathbf k\otimes\mathbf
    !! k}{2k^4}\partial(k^2)
    !! $$
    !!
    !! $$
    !! \gamma:=\frac{2.5}{\sqrt[3]{V_\text{uc}}},\quad R_\text
    !! c:=\frac6\gamma,\quad k_\text c:=10\gamma
    !! $$
    module procedure dipole_matrix_real
    module procedure dipole_matrix_complex
end interface

contains

#   define MBD_TYPE 0
#endif

#if MBD_TYPE == 0
type(matrix_re_t) function dipole_matrix_real( &
        geom, damp, ddipmat, grad) result(dipmat)
    use mbd_constants, only: ZERO => ZERO_REAL
#elif MBD_TYPE == 1
type(matrix_cplx_t) function dipole_matrix_complex( &
        geom, damp, ddipmat, grad, q) result(dipmat)
    use mbd_constants, only: ZERO => ZERO_COMPLEX
#endif

    type(geom_t), intent(inout) :: geom
    type(damping_t), intent(in) :: damp
    type(grad_request_t), intent(in), optional :: grad
#if MBD_TYPE == 0
    type(grad_matrix_re_t), intent(out), optional :: ddipmat
#elif MBD_TYPE == 1
    type(grad_matrix_cplx_t), intent(out), optional :: ddipmat
    real(dp), intent(in) :: q(3)
#endif

    real(dp) :: Rn(3), Rnij(3), Rnij_norm, T(3, 3), f_damp, &
        sigma_ij, T0(3, 3), beta_R_vdw
    integer :: i_atom, j_atom, i_cell, n(3), range_n(3), i, j, &
        n_atoms, my_i_atom, my_j_atom, i_latt
    logical :: do_ewald, is_periodic
    type(grad_matrix_re_t) :: dT, dT0, dTew
    type(grad_scalar_t) :: df
    type(grad_request_t) :: grad_ij
#if MBD_TYPE == 0
    real(dp) :: Tij(3, 3)
    type(grad_matrix_re_t) :: dTij
#elif MBD_TYPE == 1
    complex(dp) :: Tij(3, 3), exp_qR
    type(grad_matrix_cplx_t) :: dTij
#endif

    do_ewald = .false.
    is_periodic = allocated(geom%lattice)
    n_atoms = geom%siz()
    if (present(grad)) then
        grad_ij = grad
        grad_ij%dcoords = grad%dcoords .or. grad%dlattice
    end if
#ifdef WITH_SCALAPACK
    call dipmat%init(geom%idx, geom%blacs)
#else
    call dipmat%init(geom%idx)
#endif
    if (is_periodic) then
        do_ewald = geom%gamm > 0d0
        range_n = supercell_circum(geom%lattice, geom%real_space_cutoff)
    else
        range_n(:) = 0
    end if
    if (grad_ij%dcoords) allocate (dTij%dr(3, 3, 3))
    associate (my_nr => size(dipmat%idx%i_atom), my_nc => size(dipmat%idx%j_atom))
        allocate (dipmat%val(3*my_nr, 3*my_nc), source=ZERO)
        if (present(grad)) then
            if (grad%dcoords) allocate (ddipmat%dr(3*my_nr, 3*my_nc, 3), source=ZERO)
            if (grad%dlattice) then
                allocate (ddipmat%dlattice(3*my_nr, 3*my_nc, 3, 3), source=ZERO)
            end if
            if (grad%dr_vdw) then
                allocate (ddipmat%dvdw(3*my_nr, 3*my_nc), source=ZERO)
                allocate (dTij%dvdw(3, 3))
            end if
            if (grad%dsigma) then
                allocate (ddipmat%dsigma(3*my_nr, 3*my_nc), source=ZERO)
                allocate (dTij%dsigma(3, 3))
            end if
#if MBD_TYPE == 1
            if (grad%dq) then
                allocate (ddipmat%dq(3*my_nr, 3*my_nc, 3), source=ZERO)
                allocate (dTij%dq(3, 3, 3))
            end if
#endif
        end if
    end associate
    call geom%clock(11)
    n = [0, 0, -1]
    each_cell: do i_cell = 1, product(1+2*range_n)
        call shift_idx(n, -range_n, range_n)
        if (is_periodic) then
            Rn = matmul(geom%lattice, n)
        else
            Rn(:) = 0d0
        end if
        each_atom: do my_i_atom = 1, size(dipmat%idx%i_atom)
            i_atom = dipmat%idx%i_atom(my_i_atom)
            each_atom_pair: do my_j_atom = 1, size(dipmat%idx%j_atom)
                j_atom = dipmat%idx%j_atom(my_j_atom)
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                Rnij = geom%coords(:, i_atom)-geom%coords(:, j_atom)-Rn
                Rnij_norm = sqrt(sum(Rnij**2))
                if (is_periodic .and. Rnij_norm > geom%real_space_cutoff) cycle
                if (allocated(damp%R_vdw)) then
                    beta_R_vdw = damp%beta*sum(damp%R_vdw([i_atom, j_atom]))
                end if
                if (allocated(damp%sigma)) then
                    sigma_ij = damp%mayer_scaling &
                        * sqrt(sum(damp%sigma([i_atom, j_atom])**2))
                end if
                call geom%clock(13)
                select case (damp%version)
                    case ("bare")
                        T = T_bare(Rnij, dT, grad_ij%dcoords)
                    case ("dip,1mexp")
                        T = T_1mexp_coulomb(Rnij, beta_R_vdw, damp%a)
                    case ("fermi,dip")
                        f_damp = damping_fermi(Rnij, beta_R_vdw, damp%a, df, grad_ij)
                        T0 = T_bare(Rnij, dT0, grad_ij%dcoords)
                        T = damping_grad(f_damp, df, T0, dT0, dT, grad_ij)
                    case ("sqrtfermi,dip")
                        T = damping_sqrtfermi(Rnij, beta_R_vdw, damp%a)*T_bare(Rnij)
                    case ("custom,dip")
                        T = damp%damping_custom(i_atom, j_atom)*T_bare(Rnij)
                    case ("dip,custom")
                        T = damp%potential_custom(:, :, i_atom, j_atom)
                    case ("dip,gg")
                        T = T_erf_coulomb(Rnij, sigma_ij, dT, grad_ij)
                    case ("fermi,dip,gg")
                        f_damp = damping_fermi(Rnij, beta_R_vdw, damp%a, df, grad_ij)
                        call op1minus_grad(f_damp, df)
                        T0 = T_erf_coulomb(Rnij, sigma_ij, dT0, grad_ij)
                        T = damping_grad(f_damp, df, T0, dT0, dT, grad_ij)
                        do_ewald = .false.
                    case ("sqrtfermi,dip,gg")
                        T = (1d0-damping_sqrtfermi(Rnij, beta_R_vdw, damp%a)) * &
                            T_erf_coulomb(Rnij, sigma_ij)
                        do_ewald = .false.
                    case ("custom,dip,gg")
                        T = (1d0-damp%damping_custom(i_atom, j_atom)) * &
                            T_erf_coulomb(Rnij, sigma_ij)
                        do_ewald = .false.
                end select
                call geom%clock(-13)
                if (grad_ij%dr_vdw) dT%dvdw = damp%beta*dT%dvdw
                if (do_ewald) then
                    T = T &
                        + T_erfc(Rnij, geom%gamm, dTew, grad_ij) &
                        - T_bare(Rnij, dT0, grad_ij%dcoords)
                    if (grad_ij%dcoords) dT%dr = dT%dr + dTew%dr - dT0%dr
                end if
                Tij = T
                if (grad_ij%dcoords) dTij%dr = dT%dr
                if (grad_ij%dr_vdw) dTij%dvdw = dT%dvdw
                if (grad_ij%dsigma) dTij%dsigma = dT%dsigma
#if MBD_TYPE == 1
                exp_qR = exp(-IMI*(dot_product(q, Rnij)))
                Tij = T*exp_qR
                if (grad_ij%dcoords) then
                    forall (i = 1:3)
                        dTij%dr(:, :, i) = dT%dr(:, :, i)*exp_qR - IMI*q(i)*Tij
                    end forall
                end if
                if (grad_ij%dsigma) dTij%dsigma = dT%dsigma*exp_qR
                if (grad_ij%dr_vdw) dTij%dvdw = dT%dvdw*exp_qR
                if (grad_ij%dq) then
                    forall (i = 1:3)
                        dTij%dq(:, :, i) = -IMI*Rnij(i)*Tij
                    end forall
                end if
#endif
                i = 3*(my_i_atom-1)
                j = 3*(my_j_atom-1)
                associate (T_sub => dipmat%val(i+1:i+3, j+1:j+3))
                    T_sub = T_sub + Tij
                end associate
                if (.not. present(grad)) cycle
                if (grad%dcoords .and. i_atom /= j_atom) then
                    associate (dTdR_sub => ddipmat%dr(i+1:i+3, j+1:j+3, :))
                        dTdR_sub = dTdR_sub + dTij%dr
                    end associate
                end if
                if (grad%dlattice) then
                     do i_latt = 1, 3
                        associate ( &
                            dTda_sub => ddipmat%dlattice(i+1:i+3, j+1:j+3, i_latt, :) &
                        )
                            dTda_sub = dTda_sub - dTij%dr*(n(i_latt))
                        end associate
                    end do
                end if
                if (grad%dr_vdw) then
                    associate (dTdRvdw_sub => ddipmat%dvdw(i+1:i+3, j+1:j+3))
                        dTdRvdw_sub = dTdRvdw_sub + dTij%dvdw
                    end associate
                end if
                if (grad%dsigma) then
                    associate (dTdsigma_sub => ddipmat%dsigma(i+1:i+3, j+1:j+3))
                        dTdsigma_sub = dTdsigma_sub + dTij%dsigma
                    end associate
                end if
#if MBD_TYPE == 1
                if (grad%dq) then
                    associate (dTdq_sub => ddipmat%dq(i+1:i+3, j+1:j+3, :))
                        dTdq_sub = dTdq_sub + dTij%dq
                    end associate
                end if
#endif
            end do each_atom_pair
        end do each_atom
    end do each_cell
    call geom%clock(-11)
    if (do_ewald) then
#if MBD_TYPE == 0
        call add_ewald_dipole_parts_real(geom, dipmat, ddipmat, grad)
#elif MBD_TYPE == 1
        call add_ewald_dipole_parts_complex(geom, dipmat, ddipmat, grad, q)
#endif
    end if
end function

#if MBD_TYPE == 0
subroutine add_ewald_dipole_parts_real(geom, dipmat, ddipmat, grad)
    type(matrix_re_t), intent(inout) :: dipmat
    type(grad_matrix_re_t), intent(inout) :: ddipmat
#elif MBD_TYPE == 1
subroutine add_ewald_dipole_parts_complex(geom, dipmat, ddipmat, grad, q)
    type(matrix_cplx_t), intent(inout) :: dipmat
    type(grad_matrix_cplx_t), intent(inout) :: ddipmat
#endif
    type(geom_t), intent(inout) :: geom
    type(grad_request_t), intent(in), optional :: grad
#if MBD_TYPE == 1
    real(dp), intent(in) :: q(3)
#endif

    logical :: do_surface
    real(dp) :: rec_latt(3, 3), volume, G(3), Rij(3), k(3), &
        k_sq, G_Rij, latt_inv(3, 3), &
        dGdA(3), dk_sqdA, dkk_dA(3, 3), vol_prefactor, &
        k_otimes_k(3, 3), exp_k_sq_gamma, vol_kk_exp_ksq(3, 3)
    integer :: &
        i_atom, j_atom, i, j, i_xyz, m(3), i_m, &
        range_m(3), my_i_atom, my_j_atom, i_latt, a, b
#if MBD_TYPE == 0
    real(dp) :: Tij(3, 3), exp_GR, vol_exp
#elif MBD_TYPE == 1
    complex(dp) :: Tij(3, 3), exp_GR, vol_exp
    integer :: c
    real(dp) :: dkk_dq(3, 3, 3)
#endif

    latt_inv = inverse(geom%lattice)
    rec_latt = 2*pi*transpose(latt_inv)
    volume = abs(dble(product(eigvals(geom%lattice))))
    vol_prefactor = 4*pi/volume
    range_m = supercell_circum(rec_latt, geom%rec_space_cutoff)
    call geom%clock(12)
    m = [0, 0, -1]
    each_recip_vec: do i_m = 1, product(1+2*range_m)
        call shift_idx(m, -range_m, range_m)
        G = matmul(rec_latt, m)
#if MBD_TYPE == 1
        k = G + q
#elif MBD_TYPE == 0
        k = G
#endif
        k_sq = sum(k**2)
        if (sqrt(k_sq) > geom%rec_space_cutoff .or. sqrt(k_sq) < 1d-15) cycle
        exp_k_sq_gamma = exp(-k_sq/(4*geom%gamm**2))
        forall (a = 1:3, b = 1:3) k_otimes_k(a, b) = k(a)*k(b)/k_sq
        each_atom: do my_i_atom = 1, size(dipmat%idx%i_atom)
            i_atom = dipmat%idx%i_atom(my_i_atom)
            each_atom_pair: do my_j_atom = 1, size(dipmat%idx%j_atom)
                j_atom = dipmat%idx%j_atom(my_j_atom)
                Rij = geom%coords(:, i_atom)-geom%coords(:, j_atom)
                G_Rij = dot_product(G, Rij)
#if MBD_TYPE == 1
                exp_GR = exp(IMI*G_Rij)
#elif MBD_TYPE == 0
                exp_GR = cos(G_Rij)
#endif
                vol_kk_exp_ksq = vol_prefactor*k_otimes_k*exp_k_sq_gamma
                Tij = vol_kk_exp_ksq*exp_GR
                i = 3*(my_i_atom-1)
                j = 3*(my_j_atom-1)
                associate (T_sub => dipmat%val(i+1:i+3, j+1:j+3))
                    T_sub = T_sub + Tij
                end associate
                if (.not. present(grad)) cycle
                vol_exp = vol_prefactor*exp_k_sq_gamma*exp_GR
                if (grad%dcoords .and. i_atom /= j_atom) then
                    associate (dTdR_sub => ddipmat%dr(i+1:i+3, j+1:j+3, :))
                        forall (i_xyz = 1:3)
                            dTdR_sub(:, :, i_xyz) = dTdR_sub(:, :, i_xyz) &
#if MBD_TYPE == 1
                                + Tij*IMI*G(i_xyz)
#elif MBD_TYPE == 0
                                - vol_kk_exp_ksq*sin(G_Rij)*G(i_xyz)
#endif
                        end forall
                    end associate
                end if
                if (grad%dlattice) then
                    do i_latt = 1, 3
                        do i_xyz = 1, 3
                            dGdA = -latt_inv(i_latt, :)*G(i_xyz)
                            dk_sqdA = 2*dot_product(k, dGdA)
                            forall (a = 1:3, b = 1:3)
                                dkk_dA(a, b) = k(a)*dGdA(b)/k_sq &
                                    - k(a)*k(b)*dk_sqdA/(2*k_sq**2)
                            end forall
                            dkk_dA = dkk_dA + transpose(dkk_dA)
                            ! Using associate here was causing weird seg faults
                            ! with some Intel compilers, reporting i_xyz being
                            ! zero in the index, even though it printed as 1
                            ! associate ( &
                            !     dTda_sub => ddipmat%dlattice(i+1:i+3, j+1:j+3, i_latt, i_xyz) &
                            ! )
                            ddipmat%dlattice(i+1:i+3, j+1:j+3, i_latt, i_xyz) = &
                                ddipmat%dlattice(i+1:i+3, j+1:j+3, i_latt, i_xyz) &
                                - Tij*latt_inv(i_latt, i_xyz) &
                                + vol_exp*dkk_dA &
                                - Tij*dk_sqdA/(4*geom%gamm**2) &
#if MBD_TYPE == 1
                                + Tij*IMI*dot_product(dGdA, Rij)
#elif MBD_TYPE == 0
                                - vol_kk_exp_ksq*sin(G_Rij)*dot_product(dGdA, Rij)
#endif
                            ! end associate
                        end do
                    end do
                end if
#if MBD_TYPE == 1
                if (grad%dq) then
                    forall (b = 1:3, a = 1:3)
                        forall (c = 1:3) dkk_dq(b, c, a) = -2*k(a)*k(b)*k(c)/k_sq**2
                        dkk_dq(b, a, a) = dkk_dq(b, a, a) + k(b)/k_sq
                        dkk_dq(a, b, a) = dkk_dq(a, b, a) + k(b)/k_sq
                    end forall
                    associate (dTdq_sub => ddipmat%dq(i+1:i+3, j+1:j+3, :))
                        dTdq_sub = dTdq_sub + vol_exp*dkk_dq
                        forall (a = 1:3)
                            dTdq_sub(:, :, a) = dTdq_sub(:, :, a) &
                                - Tij*k(a)/(2*geom%gamm**2)
                        end forall
                    end associate
                end if
#endif
            end do each_atom_pair
        end do each_atom
    end do each_recip_vec
    ! self energy
    call dipmat%add_diag_scalar(-4*geom%gamm**3/(3*sqrt(pi)))
    ! surface term
#if MBD_TYPE == 1
    do_surface = sqrt(sum(q**2)) < 1d-15
#else
    do_surface = .true.
#endif
    if (do_surface) then
        do my_i_atom = 1, size(dipmat%idx%i_atom)
            do my_j_atom = 1, size(dipmat%idx%j_atom)
                do i_xyz = 1, 3
                    i = 3*(my_i_atom-1)+i_xyz
                    j = 3*(my_j_atom-1)+i_xyz
                    dipmat%val(i, j) = dipmat%val(i, j) + vol_prefactor/3
                    if (.not. present(grad)) cycle
                    if (grad%dlattice) then
                        ddipmat%dlattice(i, j, :, :) = ddipmat%dlattice(i, j, :, :) &
                            - vol_prefactor/3*latt_inv
                    end if
                end do
            end do
        end do
    end if
    call geom%clock(-12)
end subroutine

#if MBD_TYPE == 0
#   undef MBD_TYPE
#   define MBD_TYPE 1
#   include "mbd_dipole.F90"

function T_bare(r, dT, grad) result(T)
    !! $$
    !! T_{ab}(\mathbf r)=\frac{\partial^2}{\partial r_a\partial r_b}\frac1r=
    !! \frac{-3r_ar_b+r^2\delta_{ab}}{r^5},\qquad
    !! \frac{\partial T_{ab}(\mathbf r)}{\partial r_c}=-3\left(
    !! \frac{r_a\delta_{bc}+r_b\delta_{ca}+r_c\delta_{ab}}{r^5}-
    !! \frac{5r_ar_br_c}{r^7}
    !! \right)
    !! $$
    real(dp), intent(in) :: r(3)
    type(grad_matrix_re_t), intent(out), optional :: dT
    logical, intent(in), optional :: grad
    real(dp) :: T(3, 3)

    integer :: a, b, c
    real(dp) :: r_1, r_2, r_5, r_7

    r_2 = sum(r**2)
    r_1 = sqrt(r_2)
    r_5 = r_1**5
    forall (a = 1:3)
        T(a, a) = (-3*r(a)**2+r_2)/r_5
        forall (b = a+1:3)
            T(a, b) = -3*r(a)*r(b)/r_5
            T(b, a) = T(a, b)
        end forall
    end forall
    if (.not. present(grad)) return
    if (.not. grad) return
    allocate (dT%dr(3, 3, 3))
    r_7 = r_1**7
    forall (a = 1:3)
        dT%dr(a, a, a) = -3*(3*r(a)/r_5-5*r(a)**3/r_7)
        forall (b = a+1:3)
            dT%dr(a, a, b) = -3*(r(b)/r_5-5*r(a)**2*r(b)/r_7)
            dT%dr(a, b, a) = dT%dr(a, a, b)
            dT%dr(b, a, a) = dT%dr(a, a, b)
            dT%dr(b, b, a) = -3*(r(a)/r_5-5*r(b)**2*r(a)/r_7)
            dT%dr(b, a, b) = dT%dr(b, b, a)
            dT%dr(a, b, b) = dT%dr(b, b, a)
            forall (c = b+1:3)
                dT%dr(a, b, c) = 15*r(a)*r(b)*r(c)/r_7
                dT%dr(a, c, b) = dT%dr(a, b, c)
                dT%dr(b, a, c) = dT%dr(a, b, c)
                dT%dr(b, c, a) = dT%dr(a, b, c)
                dT%dr(c, a, b) = dT%dr(a, b, c)
                dT%dr(c, b, a) = dT%dr(a, b, c)
            end forall
        end forall
    end forall
end function

real(dp) function B_erfc(r, gamm, dB, grad) result(B)
    !! $$
    !! \begin{aligned}
    !! B(R,\gamma)
    !! &=\operatorname{erfc}(\gamma R)
    !! +\frac{2\gamma R}{\sqrt\pi}\mathrm e^{-(\gamma R)^2}
    !! \\ \partial B(R,\gamma)
    !! &=-\frac4{\sqrt\pi}(\gamma R)^2\mathrm e^{-(\gamma R)^2}
    !! (R\partial\gamma+\gamma\partial R)
    !! \end{aligned}
    !! $$
    real(dp), intent(in) :: r
    real(dp), intent(in) :: gamm
    type(grad_scalar_t), intent(out), optional :: dB
    type(grad_request_t), intent(in), optional :: grad

    real(dp) :: tmp, gamma_r_sq

    gamma_r_sq = (gamm*r)**2
    B = (erfc(gamm*r)+(2*gamm*r/sqrt(pi))*exp(-gamma_r_sq))
    if (.not. present(grad)) return
    tmp = -4d0/sqrt(pi)*gamma_r_sq*exp(-gamma_r_sq)
    if (grad%dcoords) dB%dr_1 = tmp*gamm
    if (grad%dgamma) dB%dgamma = tmp*r
end function

real(dp) function C_erfc(r, gamm, dC, grad) result(C)
    !! $$
    !! \begin{aligned}
    !! C(r,\gamma)
    !! &=3\operatorname{erfc}(\gamma R)
    !! +\frac{2\gamma R}{\sqrt\pi}(3+2(\gamma R)^2)\mathrm e^{-(\gamma R)^2}
    !! \\ \partial C(R,\gamma)
    !! &=-\frac8{\sqrt\pi}(\gamma R)^4\mathrm e^{-(\gamma R)^2}
    !! (R\partial\gamma+\gamma\partial R)
    !! \end{aligned}
    !! $$
    real(dp), intent(in) :: r
    real(dp), intent(in) :: gamm
    type(grad_scalar_t), intent(out), optional :: dC
    type(grad_request_t), intent(in), optional :: grad

    real(dp) :: tmp, gamma_r_sq

    gamma_r_sq = (gamm*r)**2
    C = (3*erfc(gamm*r)+(2*gamm*r/sqrt(pi))*(3d0+2*gamma_r_sq)*exp(-gamma_r_sq))
    if (.not. present(grad)) return
    tmp = -8d0/sqrt(pi)*gamma_r_sq**2*exp(-gamma_r_sq)
    if (grad%dcoords) dC%dr_1 = tmp*gamm
    if (grad%dgamma) dC%dgamma = tmp*r
end function

function T_erfc(r, gamm, dT, grad) result(T)
    !! $$
    !! T_{ab}^\text{erfc}(\mathbf r,\gamma)
    !! =-3\frac{r_ar_b}{r^5}C(r,\gamma)+\frac{\delta_{ab}}{r^3}B(r,\gamma)
    !! $$
    !!
    !! $$
    !! \begin{aligned}
    !! \frac{\partial T_{ab}^\text{erfc}(\mathbf r,\gamma)}{\partial r_c}
    !! &=-\left(
    !! \frac{r_a\delta_{bc}+r_b\delta_{ca}}{r^5}-
    !! 5\frac{r_ar_br_c}{r^7}
    !! \right)C(r,\gamma)-3\frac{r_c\delta_{ab}}{r^5}B(r,\gamma)
    !! \\ &-\frac{r_ar_br_c}{r^6}\frac{\partial C(r,\gamma)}{\partial
    !! r}+\frac{r_c\delta_{ab}}{r^4}\frac{\partial B(r,\gamma)}{\partial r}
    !! \end{aligned}
    !! $$
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: gamm
    type(grad_matrix_re_t), intent(out), optional :: dT
    type(grad_request_t), intent(in), optional :: grad
    real(dp) :: T(3, 3)

    integer :: a, b, c
    real(dp) :: r_1, r_2, r_3, r_4, r_5, r_6, r_7, B_ew, C_ew
    type(grad_scalar_t) :: dB, dC

    r_2 = sum(r**2)
    r_1 = sqrt(r_2)
    r_3 = r_1*r_2
    r_5 = r_3*r_2
    B_ew = B_erfc(r_1, gamm, dB, grad)
    C_ew = C_erfc(r_1, gamm, dC, grad)
    forall (a = 1:3)
        T(a, a) = -C_ew*r(a)**2/r_5+B_ew/r_3
        forall (b = a+1:3)
            T(a, b) = -C_ew*r(a)*r(b)/r_5
            T(b, a) = T(a, b)
        end forall
    end forall
    if (.not. present(grad)) return
    if (grad%dcoords) then
        allocate (dT%dr(3, 3, 3))
        r_7 = r_1**7
        r_4 = r_2**2
        r_6 = r_4*r_2
        forall (c = 1:3)
            dT%dr(c, c, c) = &
                -(2*r(c)/r_5-5*r(c)**3/r_7)*C_ew - 3*r(c)/r_5*B_ew &
                - r(c)**3/r_6*dC%dr_1 + r(c)/r_4*dB%dr_1
            forall (a = 1:3, a /= c)
                dT%dr(a, c, c) = &
                    -(r(a)/r_5-5*r(a)*r(c)**2/r_7)*C_ew &
                    - r(a)*r(c)**2/r_6*dC%dr_1
                dT%dr(c, a, c) = dT%dr(a, c, c)
                dT%dr(a, a, c) = &
                    5*r(a)**2*r(c)/r_7*C_ew - 3*r(c)/r_5*B_ew &
                    - r(a)**2*r(c)/r_6*dC%dr_1 + r(c)/r_4*dB%dr_1
                forall (b = a+1:3, b /= c)
                    dT%dr(a, b, c) = &
                        5*r(a)*r(b)*r(c)/r_7*C_ew - r(a)*r(b)*r(c)/r_6*dC%dr_1
                    dT%dr(b, a, c) = dT%dr(a, b, c)
                end forall
            end forall
        end forall
    end if
    if (grad%dgamma) then
        allocate (dT%dgamma(3, 3))
        forall (a = 1:3)
            dT%dgamma(a, a) = -dC%dgamma*r(a)**2/r_5+dB%dgamma/r_3
            forall (b = a+1:3)
                dT%dgamma(a, b) = -dC%dgamma*r(a)*r(b)/r_5
                dT%dgamma(b, a) = dT%dgamma(a, b)
            end forall
        end forall
    end if
end function

function T_erf_coulomb(r, sigma, dT, grad) result(T)
    !! $$
    !! \begin{aligned}
    !! T^\text{GG}_{ab}(\mathbf r,\sigma)&=
    !! \frac{\partial^2}{\partial r_a\partial r_b}\frac{\operatorname{erf}(\zeta)}r=
    !! \big(\operatorname{erf}(\zeta)-\Theta(\zeta)\big)T_{ab}(\mathbf r)+
    !! 2\zeta^2\Theta(\zeta)\frac{r_ar_b}{r^5}
    !! \\ \Theta(\zeta)&=\frac{2\zeta}{\sqrt\pi}\exp(-\zeta^2),\qquad
    !! \zeta=\frac r\sigma
    !! \\ \frac{\mathrm d T_{ab}^\text{GG}(\mathbf r,\sigma)}{\mathrm dr_c}&=
    !! 2\zeta\Theta(\zeta)\left(T_{ab}(\mathbf r)+(3-2\zeta^2)\frac{r_ar_b}{r^5}\right)
    !! \frac{\mathrm d\zeta}{\mathrm dr_c}
    !! \\ &+\big(\operatorname{erf}(\zeta)-\Theta(\zeta)\big)
    !! \frac{\partial T_{ab}(\mathbf r)}{\partial r_c}-
    !! 2\zeta^2\Theta(\zeta)\left(
    !! \frac13\frac{\partial T_{ab}(\mathbf r)}{\partial r_c}+
    !! \frac{r_c\delta_{ab}}{r^5}
    !! \right)
    !! \\ \qquad\frac{\mathrm d\zeta}{\mathrm dr_c}&=
    !! \frac{r_c}{r\sigma}-\frac r{\sigma^2}\frac{\mathrm d\sigma}{\mathrm dr_c}
    !! \end{aligned}
    !! $$
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: sigma
    type(grad_matrix_re_t), intent(out), optional :: dT
    type(grad_request_t), intent(in), optional :: grad
    real(dp) :: T(3, 3)

    real(dp) :: theta, erf_theta, r_5, r_1, zeta, bare(3, 3)
    type(grad_matrix_re_t) :: dbare
    real(dp) :: tmp33(3, 3), tmp333(3, 3, 3), rr_r5(3, 3)
    integer :: a, c

    bare = T_bare(r, dbare, grad%dcoords)
    r_1 = sqrt(sum(r**2))
    r_5 = r_1**5
    rr_r5 = outer(r, r)/r_5
    zeta = r_1/sigma
    theta = 2*zeta/sqrt(pi)*exp(-zeta**2)
    erf_theta = erf(zeta)-theta
    T = erf_theta*bare+2*(zeta**2)*theta*rr_r5
    if (.not. present(grad)) return
    tmp33 = 2*zeta*theta*(bare+(3-2*zeta**2)*rr_r5)
    if (grad%dcoords) then
        allocate (dT%dr(3, 3, 3))
        forall (c = 1:3) dT%dr(:, :, c) = tmp33*r(c)/(r_1*sigma)
        tmp333 = dbare%dr/3
        forall (a = 1:3, c = 1:3) tmp333(a, a, c) = tmp333(a, a, c) + r(c)/r_5
        dT%dr = dT%dr + erf_theta*dbare%dr-2*(zeta**2)*theta*tmp333
    end if
    if (grad%dsigma) dT%dsigma = -tmp33*r_1/sigma**2
end function

function T_1mexp_coulomb(rxyz, sigma, a) result(T)
    real(dp), intent(in) :: rxyz(3), sigma, a
    real(dp) :: T(3, 3)

    real(dp) :: r_sigma, zeta_1, zeta_2

    r_sigma = (sqrt(sum(rxyz**2))/sigma)**a
    zeta_1 = 1d0-exp(-r_sigma)-a*r_sigma*exp(-r_sigma)
    zeta_2 = -r_sigma*a*exp(-r_sigma)*(1+a*(-1+r_sigma))
    T = zeta_1*T_bare(rxyz)-zeta_2*outer(rxyz, rxyz)/sqrt(sum(rxyz**2))**5
end function

function damping_grad(f, df, T, dT, dfT, grad) result(fT)
    real(dp), intent(in) :: f
    type(grad_scalar_t), intent(in) :: df
    real(dp), intent(in) :: T(3, 3)
    type(grad_matrix_re_t), intent(in) :: dT
    type(grad_matrix_re_t), intent(out) :: dfT
    type(grad_request_t), intent(in) :: grad
    real(dp) :: fT(3, 3)

    integer :: c

    fT = f*T
    if (grad%dcoords) then
        allocate (dfT%dr(3, 3, 3), source=0d0)
        if (allocated(df%dr)) forall (c = 1:3) dfT%dr(:, :, c) = df%dr(c)*T
        if (allocated(dT%dr)) dfT%dr = dfT%dr + f*dT%dr
    end if
    if (grad%dr_vdw) then
        allocate (dfT%dvdw(3, 3), source=0d0)
        if (allocated(df%dvdw)) dfT%dvdw = df%dvdw*T
        if (allocated(dT%dvdw)) dfT%dvdw = dfT%dvdw + f*dT%dvdw
    end if
    if (grad%dsigma) dfT%dsigma = f*dT%dsigma
end function

end module

#endif
