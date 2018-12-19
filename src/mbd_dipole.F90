! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MBD_TYPE

module mbd_dipole
!! Construction of dipole tensors and dipole matrices.

use mbd_constants
use mbd_matrix, only: matrix_re_t, matrix_cplx_t
use mbd_geom, only: geom_t
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
    !! \begin{gathered}
    !! \mathbf T_{ij}(\mathbf q)=\sum_{\mathbf m}\mathbf T(\mathbf R_{\mathbf
    !! mij})\mathrm e^{-\mathrm i\mathbf q\cdot\mathbf R_{\mathbf mij}},\quad\mathbf
    !! R_{\mathbf mij}=\mathbf R_j+\mathbf R_\mathbf m-\mathbf R_i
    !! \\ \frac{\mathrm d\mathbf R_{\mathbf mij}}{\mathrm d\mathbf R_k}=(\delta_{jk}-\delta_{ik})\mathbf I
    !! \\ \mathbf{T}_{ij}(\mathbf{q})\approx\mathbf{T}^\text{Ew}_{ij}(\mathbf{q})
    !! =\sum_\mathbf m^{|\mathbf R_{\mathbf mij}|<R_\text c}\mathbf
    !! T^\text{erfc}(\mathbf R_{\mathbf mij},\gamma)\mathrm e^{-\mathrm i\mathbf
    !! q\cdot\mathbf R_{\mathbf mij}} +\frac{4\pi}{V_\text{uc}}\sum_{\mathbf
    !! n}^{0<|\mathbf k_\mathbf n|<k_\text c}\mathbf{\hat k}_\mathbf
    !! n\otimes\mathbf{\hat k}_\mathbf n\,\mathrm e^{-\frac{k_\mathbf
    !! n^2}{4\gamma^2}-\mathrm i\mathbf G_\mathbf n\cdot\mathbf R_{ij}}
    !! \\ -\frac{4\gamma^3}{3\sqrt\pi}\delta_{ij}\mathbf I +\delta(\mathbf q)\frac{4
    !! \pi}{3 V_\text{uc}}\mathbf I,\qquad \mathbf k_\mathbf n=\mathbf G_\mathbf
    !! n+\mathbf q
    !! \end{gathered}
    !! $$
    !!
    !! $$
    !! \mathbf a:=(\mathbf a_1,\mathbf a_2,\mathbf a_3),\quad\mathbf b:=(\mathbf
    !! b_1,\mathbf b_2,\mathbf b_3)
    !! \\ \mathbf b=2\pi(\mathbf a^{-1})^\mathrm T,\qquad \partial\mathbf
    !! b=\big((\partial\mathbf a)\mathbf a^{-1}\big)^\mathrm T\mathbf b
    !! \\ \mathbf R_\mathbf m=\mathbf a\mathbf m,\qquad\partial\mathbf R_\mathbf
    !! m=(\partial\mathbf a)\mathbf m,
    !! \\ \mathbf G_\mathbf n=\mathbf b\mathbf n,\qquad \partial\mathbf
    !! G_\mathbf n=\big((\partial\mathbf a)\mathbf a^{-1}\big)^\mathrm T\mathbf
    !! G_\mathbf n
    !! $$
    !!
    !! $$
    !! \partial\left(\frac{4\pi}{V_\text{uc}}\right)=-\frac{4\pi}{V_\text{uc}^2}\partial
    !! V_\text{uc},
    !! \\ \partial\mathbf{\hat k}_\mathbf n=\frac1{k_\mathbf
    !! n}\left(\partial\mathbf k_\mathbf n-\frac{\mathbf k_\mathbf n}{k_\mathbf
    !! n^2}\mathbf k_\mathbf n\cdot\partial\mathbf k_\mathbf n\right)
    !! \\ \partial\left(-\frac{k_\mathbf n^2}{4\gamma^2}\right)=-\frac{\mathbf
    !! k_\mathbf n\cdot\partial\mathbf k_\mathbf n}{2\gamma^2}+\frac{k_\mathbf
    !! n^2}{2\gamma^3}\partial\gamma
    !! \\ \partial\left(\frac{4\gamma^3}{3\sqrt\pi}\right)=\frac{4\gamma^2}{9\sqrt\pi}
    !! $$
    !!
    !! $$
    !! \gamma:=\frac{2.5}{\sqrt[3]{V_\text{uc}}},
    !! \qquad\partial\gamma=-\frac{2.5}3(V_\text{uc})^{-\frac43}\partial
    !! V_\text{uc}
    !! \\ R_\text c:=\frac6\gamma,\quad k_\text c:=10\gamma
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
        geom, damp, ddipmat, grad, k_point) result(dipmat)
    use mbd_constants, only: ZERO => ZERO_COMPLEX
#endif

    type(geom_t), intent(inout) :: geom
    type(damping_t), intent(in) :: damp
    type(grad_request_t), intent(in), optional :: grad
#if MBD_TYPE == 0
    type(grad_matrix_re_t), intent(out), optional :: ddipmat
#elif MBD_TYPE == 1
    type(grad_matrix_cplx_t), intent(out), optional :: ddipmat
    real(dp), intent(in) :: k_point(3)
#endif

    real(dp) :: R_cell(3), r(3), r_norm, R_vdw_ij, Tpp(3, 3), f_damp, &
        sigma_ij, volume, ewald_alpha, real_space_cutoff, T0pp(3, 3)
    integer :: i_atom, j_atom, i_cell, idx_cell(3), range_cell(3), i, j, &
        n_atoms, my_i_atom, my_j_atom
    logical :: do_ewald, is_periodic
    type(grad_matrix_re_t) :: dTpp, dT0pp
    type(grad_scalar_t) :: df
    type(grad_request_t) :: grad_
#if MBD_TYPE == 0
    real(dp) :: Tpp_ij(3, 3)
#elif MBD_TYPE == 1
    complex(dp) :: Tpp_ij(3, 3)
#endif

    do_ewald = .false.
    is_periodic = allocated(geom%lattice)
    n_atoms = geom%siz()
    if (present(grad)) grad_ = grad
#ifdef WITH_SCALAPACK
    call dipmat%init(geom%idx, geom%blacs)
#else
    call dipmat%init(geom%idx)
#endif
    if (is_periodic) then
        if (any(geom%vacuum_axis)) then
            real_space_cutoff = geom%calc%param%dipole_low_dim_cutoff
        else if (geom%calc%param%ewald_on) then
            if (grad%dcoords .or. grad%dr_vdw .or. grad%dsigma) then
                geom%calc%exc%code = MBD_EXC_UNIMPL
                geom%calc%exc%msg = 'Forces not implemented for periodic systems'
                return
            end if
            do_ewald = .true.
            volume = max(abs(dble(product(eigvals(geom%lattice)))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1d0/3)
            real_space_cutoff = &
                6d0/ewald_alpha*geom%calc%param%ewald_real_cutoff_scaling
            call geom%calc%print( &
                'Ewald: using alpha = ' // trim(tostr(ewald_alpha)) &
                // ', real cutoff = ' // trim(tostr(real_space_cutoff)) &
            )
        else
            real_space_cutoff = geom%calc%param%dipole_cutoff
        end if
        range_cell = geom%supercell_circum(geom%lattice, real_space_cutoff)
        call geom%calc%print( &
            'Ewald: summing real part in cell vector range of ' &
            // trim(tostr(1+2*range_cell(1))) // 'x' &
            // trim(tostr(1+2*range_cell(2))) // 'x' &
            // trim(tostr(1+2*range_cell(3))) &
        )
    else
        range_cell(:) = 0
    end if
    associate (my_nr => size(dipmat%idx%i_atom), my_nc => size(dipmat%idx%j_atom))
        allocate (dipmat%val(3*my_nr, 3*my_nc), source=ZERO)
        if (grad_%dcoords) allocate (ddipmat%dr(3*my_nr, 3*my_nc, 3), source=ZERO)
        if (grad_%dr_vdw) allocate (ddipmat%dvdw(3*my_nr, 3*my_nc), source=ZERO)
        if (grad_%dsigma) allocate (ddipmat%dsigma(3*my_nr, 3*my_nc), source=ZERO)
    end associate
    call geom%clock(11)
    idx_cell = [0, 0, -1]
    each_cell: do i_cell = 1, product(1+2*range_cell)
        call shift_idx(idx_cell, -range_cell, range_cell)
        if (is_periodic) then
            R_cell = matmul(geom%lattice, idx_cell)
        else
            R_cell(:) = 0d0
        end if
        each_atom: do my_i_atom = 1, size(dipmat%idx%i_atom)
            i_atom = dipmat%idx%i_atom(my_i_atom)
            each_atom_pair: do my_j_atom = 1, size(dipmat%idx%j_atom)
                j_atom = dipmat%idx%j_atom(my_j_atom)
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                r = geom%coords(:, i_atom)-geom%coords(:, j_atom)-R_cell
                r_norm = sqrt(sum(r**2))
                if (is_periodic .and. r_norm > real_space_cutoff) cycle
                if (allocated(damp%R_vdw)) then
                    R_vdw_ij = sum(damp%R_vdw([i_atom, j_atom]))
                end if
                if (allocated(damp%sigma)) then
                    sigma_ij = damp%mayer_scaling * &
                        sqrt(sum(damp%sigma([i_atom, j_atom])**2))
                end if
                select case (damp%version)
                    case ("bare")
                        Tpp = T_bare(r, dTpp, grad_%dcoords)
                    case ("dip,1mexp")
                        Tpp = T_1mexp_coulomb(r, damp%beta*R_vdw_ij, damp%a)
                    case ("fermi,dip")
                        f_damp = damping_fermi(r, damp%beta*R_vdw_ij, damp%a, df, grad_)
                        T0pp = T_bare(r, dT0pp, grad_%dcoords)
                        Tpp = damping_grad(f_damp, df, T0pp, dT0pp, dTpp, grad_)
                    case ("sqrtfermi,dip")
                        Tpp = damping_sqrtfermi(r, damp%beta*R_vdw_ij, damp%a)*T_bare(r)
                    case ("custom,dip")
                        Tpp = damp%damping_custom(i_atom, j_atom)*T_bare(r)
                    case ("dip,custom")
                        Tpp = damp%potential_custom(:, :, i_atom, j_atom)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma_ij, dTpp, grad_)
                    case ("fermi,dip,gg")
                        f_damp = damping_fermi(r, damp%beta*R_vdw_ij, damp%a, df, grad_)
                        call op1minus_grad(f_damp, df)
                        T0pp = T_erf_coulomb(r, sigma_ij, dT0pp, grad_)
                        Tpp = damping_grad(f_damp, df, T0pp, dT0pp, dTpp, grad_)
                        do_ewald = .false.
                    case ("sqrtfermi,dip,gg")
                        Tpp = (1d0-damping_sqrtfermi(r, damp%beta*R_vdw_ij, damp%a)) * &
                            T_erf_coulomb(r, sigma_ij)
                        do_ewald = .false.
                    case ("custom,dip,gg")
                        Tpp = (1d0-damp%damping_custom(i_atom, j_atom)) * &
                            T_erf_coulomb(r, sigma_ij)
                        do_ewald = .false.
                end select
                if (grad_%dr_vdw) dTpp%dvdw = damp%beta*dTpp%dvdw
                if (do_ewald) Tpp = Tpp+T_erfc(r, ewald_alpha)-T_bare(r)
#if MBD_TYPE == 0
                Tpp_ij = Tpp
#elif MBD_TYPE == 1
                Tpp_ij = Tpp*exp(-cmplx(0d0, 1d0, 8)*(dot_product(k_point, r)))
#endif
                i = 3*(my_i_atom-1)
                j = 3*(my_j_atom-1)
                associate (T => dipmat%val(i+1:i+3, j+1:j+3))
                    T = T + Tpp_ij
                end associate
                if (.not. present(grad)) cycle
#if MBD_TYPE == 0
                if (grad%dcoords) then
                    associate (T => ddipmat%dr(i+1:i+3, j+1:j+3, :))
                        T = T + dTpp%dr
                    end associate
                end if
                if (grad%dr_vdw) then
                    associate (dTdRvdw => ddipmat%dvdw(i+1:i+3, j+1:j+3))
                        dTdRvdw = dTdRvdw + dTpp%dvdw
                    end associate
                end if
                if (grad%dsigma) then
                    associate (dTdsigma => ddipmat%dsigma(i+1:i+3, j+1:j+3))
                        dTdsigma = dTdsigma + dTpp%dsigma
                    end associate
                end if
#endif
            end do each_atom_pair
        end do each_atom
    end do each_cell
    call geom%clock(-11)
    if (do_ewald) then
#if MBD_TYPE == 0
        call add_ewald_dipole_parts_real(geom, ewald_alpha, dipmat)
#elif MBD_TYPE == 1
        call add_ewald_dipole_parts_complex(geom, ewald_alpha, dipmat, k_point)
#endif
    end if
end function

#if MBD_TYPE == 0
subroutine add_ewald_dipole_parts_real(geom, alpha, dipmat)
    type(matrix_re_t), intent(inout) :: dipmat
#elif MBD_TYPE == 1
subroutine add_ewald_dipole_parts_complex(geom, alpha, dipmat, k_point)
    type(matrix_cplx_t), intent(inout) :: dipmat
#endif
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha
#if MBD_TYPE == 1
    real(dp), intent(in) :: k_point(3)
#endif

    logical :: do_surface
    real(dp) :: rec_unit_cell(3, 3), volume, G_vector(3), r(3), k_total(3), &
        k_sq, rec_space_cutoff, k_prefactor(3, 3)
    integer :: &
        i_atom, j_atom, i, j, i_xyz, j_xyz, idx_G_vector(3), i_G_vector, &
        range_G_vector(3), my_i_atom, my_j_atom
#if MBD_TYPE == 0
    real(dp) :: Tpp(3, 3)
#elif MBD_TYPE == 1
    complex(dp) :: Tpp(3, 3)
#endif

    rec_unit_cell = 2*pi*inverse(transpose(geom%lattice))
    volume = abs(dble(product(eigvals(geom%lattice))))
    rec_space_cutoff = 10d0*alpha*geom%calc%param%ewald_rec_cutoff_scaling
    range_G_vector = geom%supercell_circum(rec_unit_cell, rec_space_cutoff)
    call geom%calc%print( &
        'Ewald: using reciprocal cutoff = ' // trim(tostr(rec_space_cutoff)) &
    )
    call geom%calc%print( &
        'Ewald: summing reciprocal part in G vector range of ' &
        // trim(tostr(1+2*range_G_vector(1))) // 'x' &
        // trim(tostr(1+2*range_G_vector(2))) // 'x' &
        // trim(tostr(1+2*range_G_vector(3))) &
    )
    call geom%clock(12)
    idx_G_vector = [0, 0, -1]
    each_recip_cell: do i_G_vector = 1, product(1+2*range_G_vector)
        call shift_idx(idx_G_vector, -range_G_vector, range_G_vector)
        G_vector = matmul(rec_unit_cell, idx_G_vector)
#if MBD_TYPE == 1
        k_total = G_vector + k_point
#elif MBD_TYPE == 0
        k_total = G_vector
#endif
        k_sq = sum(k_total**2)
        if (sqrt(k_sq) > rec_space_cutoff .or. sqrt(k_sq) < 1d-15) cycle
        k_prefactor = 4*pi/volume*exp(-k_sq/(4*alpha**2))
        forall (i_xyz = 1:3, j_xyz = 1:3) &
                k_prefactor(i_xyz, j_xyz) = k_prefactor(i_xyz, j_xyz) &
                *k_total(i_xyz)*k_total(j_xyz)/k_sq
        each_atom: do my_i_atom = 1, size(dipmat%idx%i_atom)
            i_atom = dipmat%idx%i_atom(my_i_atom)
            each_atom_pair: do my_j_atom = 1, size(dipmat%idx%j_atom)
                j_atom = dipmat%idx%j_atom(my_j_atom)
                r = geom%coords(:, i_atom)-geom%coords(:, j_atom)
#if MBD_TYPE == 1
                Tpp = k_prefactor*exp(cmplx(0d0, 1d0, 8)*dot_product(G_vector, r))
#elif MBD_TYPE == 0
                Tpp = k_prefactor*cos(dot_product(G_vector, r))
#endif
                i = 3*(my_i_atom-1)
                j = 3*(my_j_atom-1)
                associate (T => dipmat%val(i+1:i+3, j+1:j+3))
                    T = T + Tpp
                end associate
            end do each_atom_pair
        end do each_atom
    end do each_recip_cell
    ! self energy
    call dipmat%add_diag_scalar(-4*alpha**3/(3*sqrt(pi)))
    ! surface term
#if MBD_TYPE == 1
    do_surface = sqrt(sum(k_point**2)) < 1d-15
#else
    do_surface = .true.
#endif
    if (do_surface) then
        do my_i_atom = 1, size(dipmat%idx%i_atom)
            do my_j_atom = 1, size(dipmat%idx%j_atom)
                do i_xyz = 1, 3
                    i = 3*(my_i_atom-1)+i_xyz
                    j = 3*(my_j_atom-1)+i_xyz
                    dipmat%val(i, j) = dipmat%val(i, j) + 4*pi/(3*volume)
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
