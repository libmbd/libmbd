! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MBD_INCLUDED
module mbd_core

use mbd_constants
use mbd_dipole, only: dipole_matrix
use mbd_matrix_type, only: mbd_matrix_real, mbd_matrix_complex, contract_cross_33
use mbd_system_type, only: mbd_system, mbd_calc
use mbd_gradients_type, only: mbd_gradients, mbd_grad_matrix_real, &
    mbd_grad_matrix_complex, mbd_grad => mbd_grad_switch
use mbd_damping_type, only: mbd_damping
use mbd_lapack, only: inverse
use mbd_common, only: tostr, findval, shift_cell
#ifdef WITH_SCALAPACK
use mbd_blacs, only: all_reduce
#endif

implicit none

#ifndef MODULE_UNIT_TESTS
private
public :: mbd_result, mbd_energy, mbd_scs_energy, scale_TS, test_frequency_grid
#endif

type :: mbd_result
    real(dp) :: energy
    real(dp), allocatable :: k_pts(:, :)
    real(dp), allocatable :: mode_eigs(:)
    real(dp), allocatable :: modes(:, :)
    real(dp), allocatable :: rpa_orders(:)
    real(dp), allocatable :: mode_eigs_k(:, :)
    complex(dp), allocatable :: modes_k(:, :, :)
    complex(dp), allocatable :: modes_k_single(:, :)
    real(dp), allocatable :: rpa_orders_k(:, :)
end type

contains

type(mbd_result) function mbd_scs_energy( &
        sys, variant, alpha_0, C6, damp, dene, grad) result(res)
    type(mbd_system), intent(inout) :: sys
    character(len=*), intent(in) :: variant
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp
    type(mbd_gradients), intent(out) :: dene
    type(mbd_grad), intent(in) :: grad

    real(dp), allocatable :: alpha_dyn(:, :), alpha_dyn_scs(:, :), &
        C6_scs(:), dC6_scs_dalpha_dyn_scs(:, :), &
        dene_dalpha_scs_dyn(:, :), freq_w(:)
    type(mbd_gradients), allocatable :: dalpha_dyn(:), dalpha_dyn_scs(:, :)
    type(mbd_gradients) :: dene_mbd, dr_vdw_scs
    type(mbd_grad) :: grad_scs
    type(mbd_damping) :: damp_scs, damp_mbd
    integer :: n_freq, i_freq, n_atoms, i_atom, my_i_atom
    character(len=15) :: damping_types(2)

    select case (variant)
    case ('scs')
        damping_types = [character(len=15) :: 'dip,gg', 'dip,1mexp']
    case ('rsscs')
        damping_types = [character(len=15) :: 'fermi,dip,gg', 'fermi,dip']
    end select
    n_freq = ubound(sys%calc%omega_grid, 1)
    n_atoms = sys%siz()
    allocate (alpha_dyn(n_atoms, 0:n_freq))
    allocate (alpha_dyn_scs(n_atoms, 0:n_freq))
    allocate (dalpha_dyn_scs(size(sys%idx%i_atom), 0:n_freq))
    if (grad%any()) allocate (dene_dalpha_scs_dyn(n_atoms, 0:n_freq))
    grad_scs = mbd_grad( &
        dcoords=grad%dcoords, dalpha=grad%dalpha .or. grad%dC6, &
        dr_vdw=grad%dr_vdw &
    )
    alpha_dyn = alpha_dynamic_ts(sys%calc, alpha_0, C6, dalpha_dyn, grad)
    damp_scs = damp
    damp_scs%version = damping_types(1)
    do i_freq = 0, n_freq
        alpha_dyn_scs(:, i_freq) = run_scs( &
            sys, alpha_dyn(:, i_freq), damp_scs, dalpha_dyn_scs(:, i_freq), grad_scs &
        )
        if (sys%has_exc()) return
    end do
    C6_scs = C6_from_alpha( &
        sys%calc, alpha_dyn_scs, dC6_scs_dalpha_dyn_scs, grad%any() &
    )
    damp_mbd = damp
    damp_mbd%r_vdw = scale_TS( &
        damp%r_vdw, alpha_dyn_scs(:, 0), alpha_dyn(:, 0), 1d0/3, dr_vdw_scs, &
        mbd_grad(dV=grad%any(), dV_free=grad%dalpha, dX_free=grad%dr_vdw) &
    )
    damp_mbd%version = damping_types(2)
    res = mbd_energy(sys, alpha_dyn_scs(:, 0), C6_scs, damp_mbd, dene_mbd, &
        mbd_grad( &
            dcoords=grad%dcoords, &
            dalpha=grad%any(), dC6=grad%any(), dr_vdw=grad%any() &
        ) &
    )
    if (sys%has_exc()) return
    if (.not. grad%any()) return
    freq_w = sys%calc%omega_grid_w
    freq_w(0) = 1d0
    dene_dalpha_scs_dyn(:, 0) = dene_mbd%dalpha + dene_mbd%dr_vdw*dr_vdw_scs%dV
    do i_freq = 1, n_freq
        dene_dalpha_scs_dyn(:, i_freq) = &
            dene_mbd%dC6*dC6_scs_dalpha_dyn_scs(:, i_freq)
    end do
    if (grad%dcoords) then
        allocate (dene%dcoords(n_atoms, 3), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = sys%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dcoords(sys%idx%j_atom, :) = &
                    dene%dcoords(sys%idx%j_atom, :) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dcoords
            end do
        end do
#ifdef WITH_SCALAPACK
        if (sys%idx%parallel) call all_reduce(dene%dcoords, sys%blacs)
#endif
        dene%dcoords = dene%dcoords + dene_mbd%dcoords
    end if
    if (grad%dalpha) then
        allocate (dene%dalpha(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = sys%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dalpha(sys%idx%j_atom) = dene%dalpha(sys%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dalpha * &
                    dalpha_dyn(i_freq)%dalpha(sys%idx%j_atom)
            end do
        end do
#ifdef WITH_SCALAPACK
        if (sys%idx%parallel) call all_reduce(dene%dalpha, sys%blacs)
#endif
        dene%dalpha = dene%dalpha + dene_mbd%dr_vdw*dr_vdw_scs%dV_free
    end if
    if (grad%dC6) then
        allocate (dene%dC6(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = sys%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dC6(sys%idx%j_atom) = dene%dC6(sys%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dalpha * &
                    dalpha_dyn(i_freq)%dC6(sys%idx%j_atom)
            end do
        end do
#ifdef WITH_SCALAPACK
        if (sys%idx%parallel) call all_reduce(dene%dC6, sys%blacs)
#endif
    end if
    if (grad%dr_vdw) then
        allocate (dene%dr_vdw(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = sys%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dr_vdw(sys%idx%j_atom) = dene%dr_vdw(sys%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dr_vdw
            end do
        end do
#ifdef WITH_SCALAPACK
        if (sys%idx%parallel) call all_reduce(dene%dr_vdw, sys%blacs)
#endif
        dene%dr_vdw = dene%dr_vdw + dene_mbd%dr_vdw*dr_vdw_scs%dX_free
    end if
end function mbd_scs_energy

#endif

#ifndef MBD_TYPE
#define MBD_TYPE 0
#endif

!> \f[
!> \begin{gathered}
!> E_\text{MBD}=\frac12\operatorname{Tr}\big(\sqrt{\mathbf Q})-
!> 3\sum_i\frac{\omega_i}2,\qquad
!> \mathbf Q_{ij}=\omega_i^2\delta_{ij}\mathbf I+
!> \omega_i\omega_j\sqrt{\alpha_{0,i}\alpha_{0,j}}\mathbf T_{ij} \\\\
!> \mathbf Q\equiv\mathbf C\boldsymbol\Lambda\mathbf C^\text T,\qquad
!> \boldsymbol\Lambda\equiv\operatorname{diag}(\{\tilde\omega_i^2\}),\qquad
!> \operatorname{Tr}\big(\sqrt{\mathbf Q}\big)=\sum_i\tilde\omega_i
!> \end{gathered}
!> \f]
!>
!> \f[
!> \begin{aligned}
!> \partial E_\text{MBD}&=\frac14\operatorname{Tr}\big(
!> \mathbf C\boldsymbol\Lambda^{-\frac12}\mathbf C^\text T
!> \partial\mathbf Q
!> \big)-
!> 3\sum_i\frac{\partial\omega_i}2 \\\\
!> \frac{\partial E_\text{MBD}}{\partial X_i}&=
!> \frac12\sum_{p\zeta}(
!> \mathbf C\boldsymbol\Lambda^{-\frac12}\mathbf C^\mathrm T
!> )_{p,i\zeta}
!> \frac{\partial Q_{p,i\zeta}}{\partial X_i}-
!> \frac32\frac{\partial\omega_i}{\partial X_i}
!> \end{aligned}
!> \f]
!>
!> \f[
!> \begin{aligned}
!> \partial\mathbf Q_{ij}=&
!> 2\delta_{ij}\omega_i\partial\omega_i\mathbf I+
!> \omega_i\omega_j\sqrt{\alpha_{0,i}\alpha_{0,j}}\mathbf T_{ij}\left(
!> \frac{\partial\omega_i}{\omega_i}+
!> \frac{\partial\omega_j}{\omega_j}+
!> \frac12\frac{\partial\alpha_{0,i}}{\alpha_{0,i}}+
!> \frac12\frac{\partial\alpha_{0,j}}{\alpha_{0,j}}
!> \right) \\\\
!> &+\omega_i\omega_j\sqrt{\alpha_{0,i}\alpha_{0,j}}
!> \partial\mathbf T_{ij}
!> \end{aligned}
!> \f]
#if MBD_TYPE == 0
type(mbd_result) function mbd_energy_single_real( &
        sys, alpha_0, C6, damp, dene, grad) result(res)
#elif MBD_TYPE == 1
type(mbd_result) function mbd_energy_single_complex( &
        sys, alpha_0, C6, damp, dene, grad, k_point) result(res)
#endif
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp
    type(mbd_gradients), intent(out) :: dene
    type(mbd_grad), intent(in) :: grad
#if MBD_TYPE == 1
    real(dp), intent(in) :: k_point(3)
#endif

#if MBD_TYPE == 0
    type(mbd_matrix_real) :: relay, dQ, T, modes, c_lambda12i_c
    type(mbd_grad_matrix_real) :: dT
    integer :: i_xyz
#elif MBD_TYPE == 1
    type(mbd_matrix_complex) :: relay, T, modes
    type(mbd_grad_matrix_complex) :: dT
#endif
    real(dp), allocatable :: eigs(:), omega(:)
    type(mbd_gradients) :: domega
    integer :: n_negative_eigs, n_atoms
    character(120) :: msg

    n_atoms = sys%siz()
#if MBD_TYPE == 0
    T = dipole_matrix(sys, damp, dT, grad)
#elif MBD_TYPE == 1
    T = dipole_matrix(sys, damp, dT, grad, k_point)
#endif
    if (sys%has_exc()) return
    if (grad%any()) then
        call relay%copy_from(T)
    else
        call relay%move_from(T)
    end if
    omega = omega_eff(C6, alpha_0, domega, grad)
    call relay%mult_cross(omega*sqrt(alpha_0))
    call relay%add_diag(omega**2)
    call sys%clock(21)
    if (sys%get_modes .or. grad%any()) then
        call modes%alloc_from(relay)
        allocate (eigs(3*n_atoms))
        call modes%eigh(eigs, sys%calc%exc, src=relay)
        if (sys%get_modes) then
#if MBD_TYPE == 0
            call move_alloc(modes%val, res%modes)
#elif MBD_TYPE == 1
            call move_alloc(modes%val, res%modes_k_single)
#endif
        end if
    else
        eigs = relay%eigvalsh(sys%calc%exc, destroy=.true.)
    end if
    if (sys%has_exc()) return
    call sys%clock(-21)
    if (sys%get_eigs) res%mode_eigs = eigs
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        msg = "CDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues"
        if (sys%calc%param%zero_negative_eigs) then
            where (eigs < 0) eigs = 0d0
            sys%calc%info%neg_eigvals = msg
        else
            sys%calc%exc%code = MBD_EXC_NEG_EIGVALS
            sys%calc%exc%msg = msg
        end if
    end if
    res%energy = 1d0/2*sum(sqrt(eigs))-3d0/2*sum(omega)
    ! TODO finish implementation for complex
#if MBD_TYPE == 0
    if (.not. grad%any()) return
    call c_lambda12i_c%copy_from(modes)
    call c_lambda12i_c%mult_cols_3n(eigs**(-1d0/4))
    c_lambda12i_c = c_lambda12i_c%mmul(c_lambda12i_c, transB=.true.)
    call dQ%init_from(T)
    if (grad%dcoords) then
        allocate (dene%dcoords(n_atoms, 3))
        do i_xyz = 1, 3
            dQ%val = dT%dr(:, :, i_xyz)
            call dQ%mult_cross(omega*sqrt(alpha_0))
            dQ%val = c_lambda12i_c%val*dQ%val
            dene%dcoords(:, i_xyz) = 1d0/2*dQ%contract_n33_rows()
        end do
    end if
    if (grad%dalpha) then
        dQ%val = T%val
        call dQ%mult_cross(omega*sqrt(alpha_0))
        call dQ%mult_rows(1d0/(2*alpha_0)+domega%dalpha/omega)
        call dQ%add_diag(omega*domega%dalpha)
        dQ%val = c_lambda12i_c%val*dQ%val
        dene%dalpha = 1d0/2*dQ%contract_n33_rows()-3d0/2*domega%dalpha
    end if
    if (grad%dC6) then
        dQ%val = T%val
        call dQ%mult_cross(omega*sqrt(alpha_0))
        call dQ%mult_rows(domega%dC6/omega)
        call dQ%add_diag(omega*domega%dC6)
        dQ%val = c_lambda12i_c%val*dQ%val
        dene%dC6 = 1d0/2*dQ%contract_n33_rows()-3d0/2*domega%dC6
    end if
    if (grad%dr_vdw) then
        dQ%val = dT%dvdw
        call dQ%mult_cross(omega*sqrt(alpha_0))
        dQ%val = c_lambda12i_c%val*dQ%val
        dene%dr_vdw = 1d0/2*dQ%contract_n33_rows()
    end if
#endif
end function

#if MBD_TYPE == 0
type(mbd_result) function rpa_energy_single_real( &
        sys, alpha, damp) result(res)
#elif MBD_TYPE == 1
type(mbd_result) function rpa_energy_single_complex( &
        sys, alpha, damp, k_point) result(res)
#endif
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha(:, 0:)
    type(mbd_damping), intent(in) :: damp
#if MBD_TYPE == 1
    real(dp), intent(in) :: k_point(3)
#endif

#if MBD_TYPE == 0
    type(mbd_matrix_real) :: relay, AT
#elif MBD_TYPE == 1
    type(mbd_matrix_complex) :: relay, AT
#endif
    complex(dp), allocatable :: eigs(:)
    integer :: i_freq, i, my_i_atom, n_order, n_negative_eigs
    type(mbd_damping) :: damp_alpha

    res%energy = 0d0
    damp_alpha = damp
    allocate (eigs(3*sys%siz()))
    do i_freq = 0, ubound(sys%calc%omega_grid, 1)
        damp_alpha%sigma = sigma_selfint(alpha(:, i_freq))
        ! relay = T
#if MBD_TYPE == 0
        relay = dipole_matrix(sys, damp_alpha)
#elif MBD_TYPE == 1
        relay = dipole_matrix(sys, damp_alpha, k_point=k_point)
#endif
        do my_i_atom = 1, size(relay%idx%i_atom)
            associate ( &
                    i_atom => relay%idx%i_atom(my_i_atom), &
                    relay_sub => relay%val(3*(my_i_atom-1)+1:, :) &
            )
                relay_sub(:3, :) = relay_sub(:3, :)*alpha(i_atom, i_freq)
            end associate
        end do
        ! relay = alpha*T
        if (sys%get_rpa_orders) AT = relay
        ! relay = 1+alpha*T
        call relay%add_diag_scalar(1d0)
        call sys%clock(23)
        eigs = relay%eigvals(sys%calc%exc, destroy=.true.)
        call sys%clock(-23)
        if (sys%has_exc()) return
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            sys%calc%exc%code = MBD_EXC_NEG_EIGVALS
            sys%calc%exc%msg = "1+AT matrix has " // &
                trim(tostr(n_negative_eigs)) // " negative eigenvalues"
            return
        end if
        res%energy = res%energy + &
            1d0/(2*pi)*sum(log(dble(eigs)))*sys%calc%omega_grid_w(i_freq)
        if (sys%get_rpa_orders) then
            call sys%clock(24)
            eigs = AT%eigvals(sys%calc%exc, destroy=.true.)
            call sys%clock(-24)
            if (sys%has_exc()) return
            allocate (res%rpa_orders(sys%calc%param%rpa_order_max))
            do n_order = 2, sys%calc%param%rpa_order_max
                res%rpa_orders(n_order) = res%rpa_orders(n_order) &
                    +(-1d0/(2*pi)*(-1)**n_order &
                    *sum(dble(eigs)**n_order)/n_order) &
                    *sys%calc%omega_grid_w(i_freq)
            end do
        end if
    end do
end function

#undef MBD_TYPE
#ifndef MBD_INCLUDED
#define MBD_INCLUDED
#define MBD_TYPE 1
#include "mbd_core.F90"
#undef MBD_INCLUDED

subroutine test_frequency_grid(calc)
    type(mbd_calc), intent(inout) :: calc

    real(dp) :: alpha(1, 0:ubound(calc%omega_grid, 1)), C6(1), error
    type(mbd_gradients), allocatable :: dalpha(:)
    type(mbd_grad) :: grad

    alpha = alpha_dynamic_ts(calc, [21d0], [99.5d0], dalpha, grad)
    C6 = C6_from_alpha(calc, alpha)
    error = abs(C6(1)/99.5d0-1d0)
    calc%info%freq_error = &
        "Relative quadrature error in C6 of carbon atom: " // &
        trim(tostr(error))
end subroutine

!> \f[
!> \bar\alpha_i=\tfrac13\operatorname{Tr}
!> \big(\textstyle\sum_j\boldsymbol{\bar\alpha}_{ij}\big),\qquad
!> \boldsymbol{\bar\alpha}=(\boldsymbol\alpha^{-1}+\mathbf T_\text{GG})^{-1}
!> \f]
!>
!> \f[
!> \begin{gathered}
!> \partial\boldsymbol{\bar\alpha}=
!> -\boldsymbol{\bar\alpha}(
!> \partial\boldsymbol\alpha^{-1}+\partial\mathbf T_\text{GG}
!> )\boldsymbol{\bar\alpha},\qquad
!> \frac{\partial\bar\alpha_i}{\partial X_j}=
!> -\frac13\sum_{\zeta\eta}\big(
!> B_{i\zeta,j\eta}\bar\alpha'_{\zeta,j\eta}+
!> B'_{j\eta,\zeta}\bar\alpha_{j\eta,i\zeta}
!> \big) \\\\
!> \mathbf B=\boldsymbol{\bar\alpha}\mathbf A,
!> \quad A_{i\zeta,j\eta}=
!> \frac{\partial(\alpha_i^{-1})}{\partial X_i}
!> \delta_{ij}\delta_{\zeta\eta}+
!> \frac{\partial T^\text{GG}_{i\zeta,j\eta}}{\partial X_i},\quad
!> \bar\alpha'_{\zeta,p}=\sum_i\bar\alpha_{i\zeta,p},\quad
!> B'_{p,\zeta}=\sum_iB_{p,i\zeta}
!> \end{gathered}
!> \f]
function run_scs(sys, alpha, damp, dalpha_scs, grad) result(alpha_scs)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha(:)
    type(mbd_damping), intent(in) :: damp
    type(mbd_gradients), intent(out) :: dalpha_scs(:)
    type(mbd_grad), intent(in) :: grad
    real(dp) :: alpha_scs(size(alpha))

    type(mbd_matrix_real) :: alpha_full, dQ, T
    integer :: n_atoms, i_xyz, i_atom, my_i_atom
    type(mbd_damping) :: damp_local
    real(dp), allocatable :: dsij_dsi(:), dsigma_dalpha(:), &
        alpha_prime(:, :), B_prime(:, :), grads_i(:)
    type(mbd_grad_matrix_real) :: dT
    type(mbd_grad) :: grad_T

    n_atoms = sys%siz()
    damp_local = damp
    damp_local%sigma = sigma_selfint(alpha, dsigma_dalpha, grad%dalpha)
    grad_T = mbd_grad(dcoords=grad%dcoords, dsigma=grad%dalpha, dr_vdw=grad%dr_vdw)
    T = dipole_matrix(sys, damp_local, dT, grad_T)
    if (sys%has_exc()) return
    if (grad%any()) then
        call alpha_full%copy_from(T)
    else
        call alpha_full%move_from(T)
    end if
    call alpha_full%add_diag(1d0/alpha)
    call sys%clock(32)
    call alpha_full%invh(sys%calc%exc)
    if (sys%has_exc()) return
    call sys%clock(-32)
    alpha_scs = alpha_full%contract_n33diag_cols()
    if (any(alpha_scs < 0)) then
        sys%calc%exc%code = MBD_EXC_NEG_POL
        sys%calc%exc%msg = 'Screening leads to negative polarizability'
        return
    end if
    if (.not. grad%any()) return
    allocate (alpha_prime(3, 3*n_atoms), source=0d0)
    allocate (B_prime(3*n_atoms, 3), source=0d0)
    allocate (grads_i(n_atoms))
    call alpha_full%contract_n_transp('R', alpha_prime)
    call dQ%init_from(T)
    if (grad%dcoords) then
        do my_i_atom = 1, size(sys%idx%i_atom)
            allocate (dalpha_scs(my_i_atom)%dcoords(size(sys%idx%j_atom), 3))
        end do
        do i_xyz = 1, 3
            dQ%val = -dT%dr(:, :, i_xyz)
            dQ = alpha_full%mmul(dQ)
            call dQ%contract_n_transp('C', B_prime)
            do i_atom = 1, n_atoms
                grads_i = contract_cross_33( &
                    i_atom, dQ, alpha_prime, alpha_full, B_prime &
                )
                my_i_atom = findval(sys%idx%i_atom, i_atom)
                if (my_i_atom > 0) then
                    dalpha_scs(my_i_atom)%dcoords(:, i_xyz) = &
                        grads_i(sys%idx%j_atom)
                end if
            end do
        end do
    end if
    if (grad%dalpha) then
        dQ%val = dT%dsigma
        do i_atom = 1, n_atoms
            dsij_dsi = damp_local%sigma(i_atom)*dsigma_dalpha(i_atom) / &
                sqrt(damp_local%sigma(i_atom)**2+damp_local%sigma**2)
            call dQ%mult_col(i_atom, dsij_dsi)
        end do
        call dQ%add_diag(-0.5d0/alpha**2)
        dQ = alpha_full%mmul(dQ)
        call dQ%contract_n_transp('C', B_prime)
        do i_atom = 1, n_atoms
            grads_i = contract_cross_33( &
                i_atom, dQ, alpha_prime, alpha_full, B_prime &
            )
            my_i_atom = findval(sys%idx%i_atom, i_atom)
            if (my_i_atom > 0) then
                dalpha_scs(my_i_atom)%dalpha = grads_i(sys%idx%j_atom)
            end if
        end do
    end if
    if (grad%dr_vdw) then
        dQ%val = dT%dvdw
        dQ = alpha_full%mmul(dQ)
        call dQ%contract_n_transp('C', B_prime)
        do i_atom = 1, n_atoms
            grads_i = contract_cross_33( &
                i_atom, dQ, alpha_prime, alpha_full, B_prime &
            )
            my_i_atom = findval(sys%idx%i_atom, i_atom)
            if (my_i_atom > 0) then
                dalpha_scs(my_i_atom)%dr_vdw = grads_i(sys%idx%j_atom)
            end if
        end do
    end if
end function run_scs

type(mbd_result) function mbd_energy( &
        sys, alpha_0, C6, damp, dene, grad) result(res)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp
    type(mbd_gradients), intent(out) :: dene
    type(mbd_grad), intent(in) :: grad

    real(dp), allocatable :: alpha(:, :)
    type(mbd_gradients), allocatable :: dalpha(:)
    integer :: n_atoms, n_kpts, i_kpt
    real(dp) :: k_point(3)
    type(mbd_result) :: res_k
    type(mbd_gradients) :: dene_k

    n_atoms = sys%siz()
    if (sys%do_rpa) then
        alpha = alpha_dynamic_ts(sys%calc, alpha_0, C6, dalpha, mbd_grad())
    end if
    if (.not. allocated(sys%lattice)) then
        if (.not. sys%do_rpa) then
            res = mbd_energy_single_real(sys, alpha_0, C6, damp, dene, grad)
        else
            res = rpa_energy_single_real(sys, alpha, damp)
            ! TODO gradients
        end if
    else
        res%k_pts = make_k_grid(make_g_grid( &
            sys%calc, sys%k_grid(1), sys%k_grid(2), sys%k_grid(3) &
        ), sys%lattice)
        n_kpts = size(res%k_pts, 2)
        res%energy = 0d0
        if (sys%get_eigs) &
            allocate (res%mode_eigs_k(3*n_atoms, n_kpts), source=0d0)
        if (sys%get_modes) &
            allocate (res%modes_k(3*n_atoms, 3*n_atoms, n_kpts), source=(0d0, 0d0))
        if (sys%get_rpa_orders) allocate ( &
            res%rpa_orders_k(sys%calc%param%rpa_order_max, n_kpts), source=0d0 &
        )
        do i_kpt = 1, n_kpts
            k_point = res%k_pts(:, i_kpt)
            if (.not. sys%do_rpa) then
                res_k = mbd_energy_single_complex( &
                    sys, alpha_0, C6, damp, dene_k, grad, k_point &
                )
                if (sys%get_eigs) res%mode_eigs_k(:, i_kpt) = res_k%mode_eigs
                if (sys%get_modes) res%modes_k(:, :, i_kpt) = res_k%modes_k_single
            else
                res_k = rpa_energy_single_complex(sys, alpha, damp, k_point)
                if (sys%get_rpa_orders) then
                    res%rpa_orders_k(:, i_kpt) = res_k%rpa_orders
                end if
            end if
            if (sys%has_exc()) return
            res%energy = res%energy + res_k%energy
        end do ! k_point loop
        res%energy = res%energy/size(res%k_pts, 2)
        if (sys%get_rpa_orders) res%rpa_orders = res%rpa_orders/n_kpts
    end if
end function mbd_energy

function alpha_dynamic_ts(calc, alpha_0, C6, dalpha, grad) result(alpha)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_gradients), allocatable, intent(out) :: dalpha(:)
    type(mbd_grad), intent(in) :: grad
    real(dp) :: alpha(size(alpha_0), 0:ubound(calc%omega_grid, 1))

    integer :: i_freq, n_atoms
    real(dp), allocatable :: omega(:)
    type(mbd_gradients) :: domega

    n_atoms = size(alpha_0)
    omega = omega_eff(C6, alpha_0, domega, grad)
    allocate (dalpha(0:ubound(alpha, 2)))
    do i_freq = 0, ubound(alpha, 2)
        alpha(:, i_freq) = alpha_osc(&
            alpha_0, omega, calc%omega_grid(i_freq), dalpha(i_freq), &
            mbd_grad(dalpha=grad%dalpha, domega=grad%dalpha .or. grad%dC6) &
        )
        if (grad%dalpha) then
            dalpha(i_freq)%dalpha = dalpha(i_freq)%dalpha + &
                dalpha(i_freq)%domega*domega%dalpha
        end if
        if (grad%dC6) then
            dalpha(i_freq)%dC6 = dalpha(i_freq)%domega*domega%dC6
        end if
        if (allocated(dalpha(i_freq)%domega)) deallocate (dalpha(i_freq)%domega)
    end do
end function

!> \f[
!> \alpha(\mathrm iu)=\frac{\alpha_0}{1+u^2/\omega^2},\qquad
!> \partial\alpha(\mathrm iu)=\alpha(\mathrm iu)\left(
!> \frac{\partial\alpha_0}{\alpha_0}+
!> \frac2\omega\frac{\partial\omega}{1+\omega^2/u^2}
!> \right)
!> \f]
function alpha_osc(alpha_0, omega, u, dalpha, grad) result(alpha)
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: omega(:)
    real(dp), intent(in) :: u
    type(mbd_gradients), intent(out), optional :: dalpha
    type(mbd_grad), intent(in), optional :: grad
    real(dp) :: alpha(size(alpha_0))

    alpha = alpha_0/(1+(u/omega)**2)
    if (.not. present(grad)) return
    if (grad%dalpha) dalpha%dalpha = alpha/alpha_0
    if (grad%domega) dalpha%domega = alpha*2d0/omega/(1d0+(omega/u)**2)
end function

!> \f[
!> \bar X=X\left(\frac{\bar\alpha_0}{\alpha_0}\right)^q,\qquad
!> \partial\bar X=\bar X\left(
!> \frac{\partial X}{X}+
!> q\frac{\partial\bar\alpha_0}{\bar\alpha_0}-
!> q\frac{\partial\alpha_0}{\alpha_0}
!> \right)
!> \f]
function scale_TS(X_free, V, V_free, q, dX, grad) result(X)
    real(dp), intent(in) :: X_free(:), V(:), V_free(:)
    real(dp), intent(in) :: q
    type(mbd_gradients), intent(out), optional :: dX
    type(mbd_grad), intent(in), optional :: grad
    real(dp) :: X(size(X_free))

    X = X_free*(V/V_free)**q
    if (.not. present(grad)) return
    if (grad%dX_free) dX%dX_free = X/X_free
    if (grad%dV) dX%dV = X*q/V
    if (grad%dV_free) dX%dV_free = -X*q/V_free
end function

!> \f[
!> \omega=\frac{4C_6}{3\alpha_{0}^2},\qquad
!> \partial\omega=\omega\left(
!> \frac{\partial C_6}{C_6}-\frac{2\partial\alpha_0}{\alpha_0}
!> \right)
!> \f]
function omega_eff(C6, alpha, domega, grad) result(omega)
    real(dp), intent(in) :: C6(:)
    real(dp), intent(in) :: alpha(:)
    type(mbd_gradients), intent(out), optional :: domega
    type(mbd_grad), intent(in), optional :: grad
    real(dp) :: omega(size(C6))

    omega = 4d0/3*C6/alpha**2
    if (.not. present(grad)) return
    if (grad%dC6) domega%dC6 = omega/C6
    if (grad%dalpha) domega%dalpha = -2*omega/alpha
end function

!> \f[
!> \sigma_i(u)=\left(\frac13\sqrt{\frac2\pi}\alpha_i(u)\right)^{\frac13},\qquad
!> \partial\sigma_i=\sigma_i\frac{\partial\alpha_i}{3\alpha_i}
!> \f]
!>
!> \f[
!> \sigma_{ij}(u)=\sqrt{\sigma_i(u)^2+\sigma_j(u)^2},\qquad
!> \partial\sigma_{ij}=
!> \frac{\sigma_i\partial\sigma_i+\sigma_j\partial\sigma_j}{\sigma_{ij}}
!> \f]
function sigma_selfint(alpha, dsigma_dalpha, grad) result(sigma)
    real(dp), intent(in) :: alpha(:)
    real(dp), allocatable, intent(out), optional :: dsigma_dalpha(:)
    logical, intent(in), optional :: grad
    real(dp) :: sigma(size(alpha))

    sigma = (sqrt(2d0/pi)*alpha/3d0)**(1d0/3)
    if (.not. present(grad)) return
    if (grad) dsigma_dalpha = sigma/(3*alpha)
end function

!> \f[
!> \bar C_6=\frac3\pi\int_0^\infty\mathrm du\,\bar\alpha(u)^2,\qquad
!> \partial\bar C_6=\frac6\pi\int_0^\infty\mathrm du
!> \bar\alpha(u)\partial\bar\alpha(u)
!> \f]
function C6_from_alpha(calc, alpha, dC6_dalpha, grad) result(C6)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(in) :: alpha(:, 0:)
    real(dp), allocatable, intent(out), optional :: dC6_dalpha(:, :)
    logical, intent(in), optional :: grad
    real(dp) :: C6(size(alpha, 1))

    integer :: i_freq, n_atoms

    n_atoms = size(alpha, 1)
    C6 = 0d0
    do i_freq = 0, ubound(alpha, 2)
        C6 = C6 + 3d0/pi*alpha(:, i_freq)**2*calc%omega_grid_w(i_freq)
    end do
    if (.not. present(grad)) return
    if (.not. grad) return
    allocate (dC6_dalpha(n_atoms, 0:ubound(alpha, 2)), source=0d0)
    do i_freq = 0, ubound(alpha, 2)
        dC6_dalpha(:, i_freq) = dC6_dalpha(:, i_freq) + 6d0/pi*alpha(:, i_freq)
    end do
end function

function make_g_grid(calc, n1, n2, n3) result(g_grid)
    type(mbd_calc), intent(in) :: calc
    integer, intent(in) :: n1, n2, n3
    real(dp) :: g_grid(3, n1*n2*n3)

    integer :: g_kpt(3), i_kpt, kpt_range(3)
    real(dp) :: g_kpt_shifted(3)

    g_kpt = [0, 0, -1]
    kpt_range = [n1, n2, n3]
    do i_kpt = 1, n1*n2*n3
        call shift_cell (g_kpt, [0, 0, 0], kpt_range-1)
        g_kpt_shifted = dble(g_kpt)+calc%param%k_grid_shift
        where (2*g_kpt_shifted > kpt_range)
            g_kpt_shifted = g_kpt_shifted-dble(kpt_range)
        end where
        g_grid(:, i_kpt) = g_kpt_shifted/kpt_range
    end do
end function make_g_grid

function make_k_grid(g_grid, uc) result(k_grid)
    real(dp), intent(in) :: g_grid(:, :), uc(3, 3)
    real(dp) :: k_grid(3, size(g_grid, 2))

    integer :: i_kpt
    real(dp) :: ruc(3, 3)

    ruc = 2*pi*inverse(transpose(uc))
    do i_kpt = 1, size(g_grid, 2)
        k_grid(:, i_kpt) = matmul(ruc, g_grid(:, i_kpt))
    end do
end function make_k_grid

end module

#endif
