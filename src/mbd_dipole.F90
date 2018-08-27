! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MBD_INCLUDED
module mbd_dipole

use mbd_constants
use mbd_matrix_type, only: mbd_matrix_real, mbd_matrix_complex
use mbd_system_type, only: mbd_system
use mbd_damping_type, only: mbd_damping, damping_fermi, damping_sqrtfermi, &
    op1minus_grad
use mbd_gradients_type, only: mbd_gradients, mbd_grad_matrix_real, &
    mbd_grad_matrix_complex, mbd_grad_scalar, mbd_grad_switch
use mbd_lapack, only: eigvals, inverse
use mbd_linalg, only: outer
use mbd_common, only: tostr, shift_cell

implicit none

#ifndef MODULE_UNIT_TESTS
private
public :: dipole_matrix
#endif

interface dipole_matrix
    module procedure dipole_matrix_real
    module procedure dipole_matrix_complex
end interface

contains

#endif

#ifndef MBD_TYPE
#define MBD_TYPE 0
#endif

!> \f[
!> T_{(ij)ab}\equiv\frac{\partial^2}{\partial R_a\partial R_b}\frac1R=
!> \frac{-3R_aR_b+R^2\delta_{ab}}{R^5},\qquad
!> \mathbf R\equiv \mathbf R_{ij}\equiv\mathbf R_j-\mathbf R_i
!> \f]
!>
!> \f[
!> \frac{\partial\mathbf T_{ij}}{\partial\mathbf R_k}=
!> \frac{\partial\mathbf T}{\partial\mathbf R}(\delta_{jk}-\delta_{ik}),\qquad
!> \frac{\partial T_{ab}}{\partial R_c}=-3\left(
!> \frac{R_a\delta_{bc}+R_b\delta_{ca}+R_c\delta_{ab}}{R^5}-
!> \frac{5R_aR_bR_c}{R^7}
!> \right)
!> \f]
#if MBD_TYPE == 0
type(mbd_matrix_real) function dipole_matrix_real( &
        sys, damp, ddipmat, grad) result(dipmat)
    use mbd_constants, only: ZERO => ZERO_REAL
#elif MBD_TYPE == 1
type(mbd_matrix_complex) function dipole_matrix_complex( &
        sys, damp, ddipmat, grad, k_point) result(dipmat)
    use mbd_constants, only: ZERO => ZERO_COMPLEX
#endif

    type(mbd_system), intent(inout) :: sys
    type(mbd_damping), intent(in) :: damp
    type(mbd_grad_switch), intent(in), optional :: grad
#if MBD_TYPE == 0
    type(mbd_grad_matrix_real), intent(out), optional :: ddipmat
#elif MBD_TYPE == 1
    type(mbd_grad_matrix_complex), intent(out), optional :: ddipmat
    real(dp), intent(in) :: k_point(3)
#endif

    real(dp) :: R_cell(3), r(3), r_norm, R_vdw_ij, Tpp(3, 3), f_damp, &
        sigma_ij, volume, ewald_alpha, real_space_cutoff, T0pp(3, 3)
    integer :: i_atom, j_atom, i_cell, idx_cell(3), range_cell(3), i, j, &
        n_atoms, my_i_atom, my_j_atom
    logical :: do_ewald, is_periodic
    type(mbd_grad_matrix_real) :: dTpp, dT0pp
    type(mbd_grad_scalar) :: df
    type(mbd_grad_switch) :: grad_
#if MBD_TYPE == 0
    real(dp) :: Tpp_final(3, 3)
#elif MBD_TYPE == 1
    complex(dp) :: Tpp_final(3, 3)
#endif

    do_ewald = .false.
    is_periodic = allocated(sys%lattice)
    n_atoms = sys%siz()
    if (present(grad)) grad_ = grad
#ifdef WITH_SCALAPACK
    call dipmat%init(sys%idx, sys%blacs)
#else
    call dipmat%init(sys%idx)
#endif
    if (is_periodic) then
        if (any(sys%vacuum_axis)) then
            real_space_cutoff = sys%calc%param%dipole_low_dim_cutoff
        else if (sys%calc%param%ewald_on) then
            if (grad%dcoords .or. grad%dr_vdw .or. grad%dsigma) then
                sys%calc%exc%code = MBD_EXC_UNIMPL
                sys%calc%exc%msg = 'Forces not implemented for periodic systems'
                return
            end if
            do_ewald = .true.
            volume = max(abs(dble(product(eigvals(sys%lattice)))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1d0/3)
            real_space_cutoff = &
                6d0/ewald_alpha*sys%calc%param%ewald_real_cutoff_scaling
            sys%calc%info%ewald_alpha = &
                'Ewald: using alpha = ' // trim(tostr(ewald_alpha)) // &
                ', real cutoff = ' // trim(tostr(real_space_cutoff))
        else
            real_space_cutoff = sys%calc%param%dipole_cutoff
        end if
        range_cell = sys%supercell_circum(sys%lattice, real_space_cutoff)
    else
        range_cell(:) = 0
    end if
    if (is_periodic) then
        sys%calc%info%ewald_rsum = &
            'Ewald: summing real part in cell vector range of ' // &
            trim(tostr(1+2*range_cell(1))) // 'x' // &
            trim(tostr(1+2*range_cell(2))) // 'x' // &
            trim(tostr(1+2*range_cell(3)))
    end if
    associate (my_nr => size(dipmat%idx%i_atom), my_nc => size(dipmat%idx%j_atom))
        allocate (dipmat%val(3*my_nr, 3*my_nc), source=ZERO)
        if (grad_%dcoords) allocate (ddipmat%dr(3*my_nr, 3*my_nc, 3), source=ZERO)
        if (grad_%dr_vdw) allocate (ddipmat%dvdw(3*my_nr, 3*my_nc), source=ZERO)
        if (grad_%dsigma) allocate (ddipmat%dsigma(3*my_nr, 3*my_nc), source=ZERO)
    end associate
    call sys%clock(11)
    idx_cell = [0, 0, -1]
    do i_cell = 1, product(1+2*range_cell)
        call shift_cell(idx_cell, -range_cell, range_cell)
        if (is_periodic) then
            R_cell = matmul(sys%lattice, idx_cell)
        else
            R_cell(:) = 0d0
        end if
        do my_i_atom = 1, size(dipmat%idx%i_atom)
            i_atom = dipmat%idx%i_atom(my_i_atom)
            do my_j_atom = 1, size(dipmat%idx%j_atom)
                j_atom = dipmat%idx%j_atom(my_j_atom)
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                r = sys%coords(:, i_atom)-sys%coords(:, j_atom)-R_cell
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
                if (grad_%dr_vdw) then
                    dTpp%dvdw = damp%beta*dTpp%dvdw
                end if
                if (do_ewald) then
                    Tpp = Tpp+T_erfc(r, ewald_alpha)-T_bare(r)
                end if
#if MBD_TYPE == 0
                Tpp_final = Tpp
#elif MBD_TYPE == 1
                Tpp_final = Tpp*exp(-cmplx(0d0, 1d0, 8)*(dot_product(k_point, r)))
#endif
                i = 3*(my_i_atom-1)
                j = 3*(my_j_atom-1)
                associate (T => dipmat%val(i+1:i+3, j+1:j+3))
                    T = T + Tpp_final
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
            end do ! j_atom
        end do ! i_atom
    end do ! i_cell
    call sys%clock(-11)
    if (do_ewald) then
#if MBD_TYPE == 0
        call add_ewald_dipole_parts_real(sys, ewald_alpha, dipmat)
#elif MBD_TYPE == 1
        call add_ewald_dipole_parts_complex(sys, ewald_alpha, dipmat, k_point)
#endif
    end if
end function

#if MBD_TYPE == 0
subroutine add_ewald_dipole_parts_real(sys, alpha, dipmat)
    type(mbd_matrix_real), intent(inout) :: dipmat
#elif MBD_TYPE == 1
subroutine add_ewald_dipole_parts_complex(sys, alpha, dipmat, k_point)
    type(mbd_matrix_complex), intent(inout) :: dipmat
#endif
    type(mbd_system), intent(inout) :: sys
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
    real(dp) :: elem
#endif

    rec_unit_cell = 2*pi*inverse(transpose(sys%lattice))
    volume = abs(dble(product(eigvals(sys%lattice))))
    rec_space_cutoff = 10d0*alpha*sys%calc%param%ewald_rec_cutoff_scaling
    range_G_vector = sys%supercell_circum(rec_unit_cell, rec_space_cutoff)
    sys%calc%info%ewald_cutoff = 'Ewald: using reciprocal cutoff = ' // &
        trim(tostr(rec_space_cutoff))
    sys%calc%info%ewald_recsum = &
        'Ewald: summing reciprocal part in G vector range of ' // &
        trim(tostr(1+2*range_G_vector(1))) // 'x' // &
        trim(tostr(1+2*range_G_vector(2))) // 'x' // &
        trim(tostr(1+2*range_G_vector(3)))
    call sys%clock(12)
    idx_G_vector = [0, 0, -1]
    do i_G_vector = 1, product(1+2*range_G_vector)
        call shift_cell(idx_G_vector, -range_G_vector, range_G_vector)
        if (i_G_vector == 1) cycle
        G_vector = matmul(rec_unit_cell, idx_G_vector)
#if MBD_TYPE == 1
        k_total = G_vector + k_point
#elif MBD_TYPE == 0
        k_total = G_vector
#endif
        k_sq = sum(k_total**2)
        if (sqrt(k_sq) > rec_space_cutoff) cycle
        k_prefactor = 4*pi/volume*exp(-k_sq/(4*alpha**2))
        forall (i_xyz = 1:3, j_xyz = 1:3) &
                k_prefactor(i_xyz, j_xyz) = k_prefactor(i_xyz, j_xyz) &
                *k_total(i_xyz)*k_total(j_xyz)/k_sq
        do my_i_atom = 1, size(dipmat%idx%i_atom)
            i_atom = dipmat%idx%i_atom(my_i_atom)
            do my_j_atom = 1, size(dipmat%idx%j_atom)
                j_atom = dipmat%idx%j_atom(my_j_atom)
                r = sys%coords(:, i_atom)-sys%coords(:, j_atom)
#if MBD_TYPE == 1
                Tpp = k_prefactor*exp(cmplx(0d0, 1d0, 8) &
                    *dot_product(G_vector, r))
#elif MBD_TYPE == 0
                Tpp = k_prefactor*cos(dot_product(G_vector, r))
#endif
                i = 3*(my_i_atom-1)
                j = 3*(my_j_atom-1)
                associate (T => dipmat%val(i+1:i+3, j+1:j+3))
                    T = T + Tpp
                end associate
            end do ! j_atom
        end do ! i_atom
    end do ! i_G_vector
    call dipmat%add_diag_scalar(-4*alpha**3/(3*sqrt(pi))) ! self energy
    do_surface = .true.
#if MBD_TYPE == 1
    k_sq = sum(k_point**2)
    if (sqrt(k_sq) > 1.d-15) then
        do_surface = .false.
        do my_i_atom = 1, size(dipmat%idx%i_atom)
        do my_j_atom = 1, size(dipmat%idx%j_atom)
            do i_xyz = 1, 3
            do j_xyz = 1, 3
                i = 3*(my_i_atom-1)+i_xyz
                j = 3*(my_j_atom-1)+j_xyz
                elem = 4*pi/volume*k_point(i_xyz)*k_point(j_xyz)/k_sq &
                    *exp(-k_sq/(4*alpha**2))
                dipmat%val(i, j) = dipmat%val(i, j) + elem
            end do ! j_xyz
            end do ! i_xyz
        end do ! j_atom
        end do ! i_atom
    end if ! k_sq >
#endif
    if (do_surface) then ! surface energy
        do my_i_atom = 1, size(dipmat%idx%i_atom)
        do my_j_atom = 1, size(dipmat%idx%j_atom)
            do i_xyz = 1, 3
                i = 3*(my_i_atom-1)+i_xyz
                j = 3*(my_j_atom-1)+i_xyz
                dipmat%val(i, j) = dipmat%val(i, j) + 4*pi/(3*volume)
            end do ! i_xyz
        end do ! j_atom
        end do ! i_atom
    end if
    call sys%clock(-12)
end subroutine

#undef MBD_TYPE
#ifndef MBD_INCLUDED
#define MBD_INCLUDED
#define MBD_TYPE 1
#include "mbd_dipole.F90"
#undef MBD_INCLUDED

function T_bare(r, dT, grad) result(T)
    real(dp), intent(in) :: r(3)
    type(mbd_grad_matrix_real), intent(out), optional :: dT
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

real(dp) function B_erfc(r, a) result(B)
    real(dp), intent(in) :: r, a

    B = (erfc(a*r)+(2*a*r/sqrt(pi))*exp(-(a*r)**2))/r**3
end function

real(dp) elemental function C_erfc(r, a) result(C)
    real(dp), intent(in) :: r, a

    C = (3*erfc(a*r)+(2*a*r/sqrt(pi))*(3d0+2*(a*r)**2)*exp(-(a*r)**2))/r**5
end function

function T_erfc(rxyz, alpha) result(T)
    real(dp), intent(in) :: rxyz(3), alpha
    real(dp) :: T(3, 3)

    integer :: i, j
    real(dp) :: r, B, C

    r = sqrt(sum(rxyz(:)**2))
    B = B_erfc(r, alpha)
    C = C_erfc(r, alpha)
    do i = 1, 3
        do j = i, 3
            T(i, j) = -C*rxyz(i)*rxyz(j)
            if (i /= j) T(j, i) = T(i, j)
        end do
        T(i, i) = T(i, i)+B
    end do
end function

!> \f[
!> \begin{gathered}
!> T^\text{GG}_{(ij)ab}\equiv
!> \frac{\partial^2}{\partial R_a\partial R_b}\frac{\operatorname{erf}(\zeta)}R=
!> \big(\operatorname{erf}(\zeta)-\Theta(\zeta)\big)T_{ab}+
!> 2\zeta^2\Theta(\zeta)\frac{R_aR_b}{R^5}, \\\\
!> \Theta(\zeta)=\frac{2\zeta}{\sqrt\pi}\exp(-\zeta^2),\qquad
!> \zeta=\frac{R_{(ij)}}{\sigma_{(ij)}}
!> \end{gathered}
!> \f]
!>
!> \f[
!> \begin{aligned}
!> \frac{\mathrm d T_{ab}^\text{GG}}{\mathrm dR_c}&=
!> 2\zeta\Theta(\zeta)\left(T_{ab}+(3-2\zeta^2)\frac{R_aR_b}{R^5}\right)
!> \frac{\mathrm d\zeta}{\mathrm dR_c} \\\\
!> &\qquad+\big(\operatorname{erf}(\zeta)-\Theta(\zeta)\big)
!> \frac{\partial T_{ab}}{\partial R_c}-
!> 2\zeta^2\Theta(\zeta)\left(
!> \frac13\frac{\partial T_{ab}}{\partial R_c}+
!> \frac{R_c\delta_{ab}}{R^5}
!> \right) \\\\
!> \qquad\frac{\mathrm d\zeta}{\mathrm dR_c}&=
!> \frac{R_c}{R\sigma}-\frac R{\sigma^2}\frac{\mathrm d\sigma}{\mathrm dR_c}
!> \end{aligned}
!> \f]
function T_erf_coulomb(r, sigma, dT, grad) result(T)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: sigma
    type(mbd_grad_matrix_real), intent(out), optional :: dT
    type(mbd_grad_switch), intent(in), optional :: grad
    real(dp) :: T(3, 3)

    real(dp) :: theta, erf_theta, r_5, r_1, zeta, bare(3, 3)
    type(mbd_grad_matrix_real) :: dbare
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
    type(mbd_grad_scalar), intent(in) :: df
    real(dp), intent(in) :: T(3, 3)
    type(mbd_grad_matrix_real), intent(in) :: dT
    type(mbd_grad_matrix_real), intent(out) :: dfT
    type(mbd_grad_switch), intent(in) :: grad
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
