! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MBD_INCLUDED
module mbd_core

use mbd_constants
use mbd_common, only: omega_eff, sigma_selfint
use mbd_dipole, only: dipole_matrix
use mbd_matrix, only: matrix_re_t, matrix_cplx_t, contract_cross_33
use mbd_calc, only: calc_t
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_t, grad_matrix_re_t, grad_matrix_cplx_t, &
    grad_scalar_t, grad_request_t
use mbd_hamiltonian, only: get_mbd_hamiltonian_energy
use mbd_damping, only: damping_t
use mbd_lapack, only: inverse
use mbd_rpa, only: get_mbd_rpa_energy
use mbd_utils, only: tostr, findval, shift_cell, result_t
#ifdef WITH_SCALAPACK
use mbd_blacs, only: all_reduce
#endif

implicit none

#ifndef MODULE_UNIT_TESTS
private
public :: mbd_energy, mbd_scs_energy, scale_TS, test_frequency_grid
#endif

contains

type(result_t) function mbd_scs_energy( &
        geom, variant, alpha_0, C6, damp, dene, grad) result(res)
    type(geom_t), intent(inout) :: geom
    character(len=*), intent(in) :: variant
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    type(grad_t), intent(out) :: dene
    type(grad_request_t), intent(in) :: grad

    real(dp), allocatable :: alpha_dyn(:, :), alpha_dyn_scs(:, :), &
        C6_scs(:), dC6_scs_dalpha_dyn_scs(:, :), &
        dene_dalpha_scs_dyn(:, :), freq_w(:)
    type(grad_t), allocatable :: dalpha_dyn(:), dalpha_dyn_scs(:, :)
    type(grad_t) :: dene_mbd, dr_vdw_scs
    type(grad_request_t) :: grad_scs
    type(damping_t) :: damp_scs, damp_mbd
    integer :: n_freq, i_freq, n_atoms, i_atom, my_i_atom
    character(len=15) :: damping_types(2)

    select case (variant)
    case ('scs')
        damping_types = [character(len=15) :: 'dip,gg', 'dip,1mexp']
    case ('rsscs')
        damping_types = [character(len=15) :: 'fermi,dip,gg', 'fermi,dip']
    end select
    n_freq = ubound(geom%calc%omega_grid, 1)
    n_atoms = geom%siz()
    allocate (alpha_dyn(n_atoms, 0:n_freq))
    allocate (alpha_dyn_scs(n_atoms, 0:n_freq))
    allocate (dalpha_dyn_scs(size(geom%idx%i_atom), 0:n_freq))
    if (grad%any()) allocate (dene_dalpha_scs_dyn(n_atoms, 0:n_freq))
    grad_scs = grad_request_t( &
        dcoords=grad%dcoords, dalpha=grad%dalpha .or. grad%dC6, &
        dr_vdw=grad%dr_vdw &
    )
    alpha_dyn = alpha_dynamic_ts(geom%calc, alpha_0, C6, dalpha_dyn, grad)
    damp_scs = damp
    damp_scs%version = damping_types(1)
    do i_freq = 0, n_freq
        alpha_dyn_scs(:, i_freq) = run_scs( &
            geom, alpha_dyn(:, i_freq), damp_scs, dalpha_dyn_scs(:, i_freq), grad_scs &
        )
        if (geom%has_exc()) return
    end do
    C6_scs = C6_from_alpha( &
        geom%calc, alpha_dyn_scs, dC6_scs_dalpha_dyn_scs, grad%any() &
    )
    damp_mbd = damp
    damp_mbd%r_vdw = scale_TS( &
        damp%r_vdw, alpha_dyn_scs(:, 0), alpha_dyn(:, 0), 1d0/3, dr_vdw_scs, &
        grad_request_t(dV=grad%any(), dV_free=grad%dalpha, dX_free=grad%dr_vdw) &
    )
    damp_mbd%version = damping_types(2)
    res = mbd_energy(geom, alpha_dyn_scs(:, 0), C6_scs, damp_mbd, dene_mbd, &
        grad_request_t( &
            dcoords=grad%dcoords, &
            dalpha=grad%any(), dC6=grad%any(), dr_vdw=grad%any() &
        ) &
    )
    if (geom%has_exc()) return
    if (.not. grad%any()) return
    freq_w = geom%calc%omega_grid_w
    freq_w(0) = 1d0
    dene_dalpha_scs_dyn(:, 0) = dene_mbd%dalpha + dene_mbd%dr_vdw*dr_vdw_scs%dV
    do i_freq = 1, n_freq
        dene_dalpha_scs_dyn(:, i_freq) = &
            dene_mbd%dC6*dC6_scs_dalpha_dyn_scs(:, i_freq)
    end do
    if (grad%dcoords) then
        allocate (dene%dcoords(n_atoms, 3), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dcoords(geom%idx%j_atom, :) = &
                    dene%dcoords(geom%idx%j_atom, :) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dcoords
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call all_reduce(dene%dcoords, geom%blacs)
#endif
        dene%dcoords = dene%dcoords + dene_mbd%dcoords
    end if
    if (grad%dalpha) then
        allocate (dene%dalpha(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dalpha(geom%idx%j_atom) = dene%dalpha(geom%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dalpha * &
                    dalpha_dyn(i_freq)%dalpha(geom%idx%j_atom)
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call all_reduce(dene%dalpha, geom%blacs)
#endif
        dene%dalpha = dene%dalpha + dene_mbd%dr_vdw*dr_vdw_scs%dV_free
    end if
    if (grad%dC6) then
        allocate (dene%dC6(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dC6(geom%idx%j_atom) = dene%dC6(geom%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dalpha * &
                    dalpha_dyn(i_freq)%dC6(geom%idx%j_atom)
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call all_reduce(dene%dC6, geom%blacs)
#endif
    end if
    if (grad%dr_vdw) then
        allocate (dene%dr_vdw(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dr_vdw(geom%idx%j_atom) = dene%dr_vdw(geom%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dr_vdw
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call all_reduce(dene%dr_vdw, geom%blacs)
#endif
        dene%dr_vdw = dene%dr_vdw + dene_mbd%dr_vdw*dr_vdw_scs%dX_free
    end if
end function mbd_scs_energy

#endif

#ifndef MBD_TYPE
#define MBD_TYPE 0
#endif


#undef MBD_TYPE
#ifndef MBD_INCLUDED
#define MBD_INCLUDED
#define MBD_TYPE 1
#include "mbd_core.F90"
#undef MBD_INCLUDED

subroutine test_frequency_grid(calc)
    type(calc_t), intent(inout) :: calc

    real(dp) :: alpha(1, 0:ubound(calc%omega_grid, 1)), C6(1), error
    type(grad_t), allocatable :: dalpha(:)
    type(grad_request_t) :: grad

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
function run_scs(geom, alpha, damp, dalpha_scs, grad) result(alpha_scs)
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha(:)
    type(damping_t), intent(in) :: damp
    type(grad_t), intent(out) :: dalpha_scs(:)
    type(grad_request_t), intent(in) :: grad
    real(dp) :: alpha_scs(size(alpha))

    type(matrix_re_t) :: alpha_full, dQ, T
    integer :: n_atoms, i_xyz, i_atom, my_i_atom
    type(damping_t) :: damp_local
    real(dp), allocatable :: dsij_dsi(:), dsigma_dalpha(:), &
        alpha_prime(:, :), B_prime(:, :), grads_i(:)
    type(grad_matrix_re_t) :: dT
    type(grad_request_t) :: grad_req

    n_atoms = geom%siz()
    damp_local = damp
    damp_local%sigma = sigma_selfint(alpha, dsigma_dalpha, grad%dalpha)
    grad_req = grad_request_t(dcoords=grad%dcoords, dsigma=grad%dalpha, dr_vdw=grad%dr_vdw)
    T = dipole_matrix(geom, damp_local, dT, grad_req)
    if (geom%has_exc()) return
    if (grad%any()) then
        call alpha_full%copy_from(T)
    else
        call alpha_full%move_from(T)
    end if
    call alpha_full%add_diag(1d0/alpha)
    call geom%clock(32)
    call alpha_full%invh(geom%calc%exc)
    if (geom%has_exc()) return
    call geom%clock(-32)
    alpha_scs = alpha_full%contract_n33diag_cols()
    if (any(alpha_scs < 0)) then
        geom%calc%exc%code = MBD_EXC_NEG_POL
        geom%calc%exc%msg = 'Screening leads to negative polarizability'
        return
    end if
    if (.not. grad%any()) return
    allocate (alpha_prime(3, 3*n_atoms), source=0d0)
    allocate (B_prime(3*n_atoms, 3), source=0d0)
    allocate (grads_i(n_atoms))
    call alpha_full%contract_n_transp('R', alpha_prime)
    call dQ%init_from(T)
    if (grad%dcoords) then
        do my_i_atom = 1, size(geom%idx%i_atom)
            allocate (dalpha_scs(my_i_atom)%dcoords(size(geom%idx%j_atom), 3))
        end do
        do i_xyz = 1, 3
            dQ%val = -dT%dr(:, :, i_xyz)
            dQ = alpha_full%mmul(dQ)
            call dQ%contract_n_transp('C', B_prime)
            do i_atom = 1, n_atoms
                grads_i = contract_cross_33( &
                    i_atom, dQ, alpha_prime, alpha_full, B_prime &
                )
                my_i_atom = findval(geom%idx%i_atom, i_atom)
                if (my_i_atom > 0) then
                    dalpha_scs(my_i_atom)%dcoords(:, i_xyz) = &
                        grads_i(geom%idx%j_atom)
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
            my_i_atom = findval(geom%idx%i_atom, i_atom)
            if (my_i_atom > 0) then
                dalpha_scs(my_i_atom)%dalpha = grads_i(geom%idx%j_atom)
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
            my_i_atom = findval(geom%idx%i_atom, i_atom)
            if (my_i_atom > 0) then
                dalpha_scs(my_i_atom)%dr_vdw = grads_i(geom%idx%j_atom)
            end if
        end do
    end if
end function run_scs

type(result_t) function mbd_energy( &
        geom, alpha_0, C6, damp, dene, grad) result(res)
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    type(grad_t), intent(out) :: dene
    type(grad_request_t), intent(in) :: grad

    real(dp), allocatable :: alpha(:, :)
    type(grad_t), allocatable :: dalpha(:)
    integer :: n_atoms, n_kpts, i_kpt
    real(dp) :: k_point(3)
    type(result_t) :: res_k
    type(grad_t) :: dene_k

    n_atoms = geom%siz()
    if (geom%calc%do_rpa) then
        alpha = alpha_dynamic_ts(geom%calc, alpha_0, C6, dalpha, grad_request_t())
    end if
    if (.not. allocated(geom%lattice)) then
        if (.not. geom%calc%do_rpa) then
            res = get_mbd_hamiltonian_energy(geom, alpha_0, C6, damp, dene, grad)
        else
            res = get_mbd_rpa_energy(geom, alpha, damp)
            ! TODO gradients
        end if
    else
        call geom%ensure_k_pts()
        n_kpts = size(geom%k_pts, 2)
        res%energy = 0d0
        if (geom%calc%get_eigs) &
            allocate (res%mode_eigs_k(3*n_atoms, n_kpts), source=0d0)
        if (geom%calc%get_modes) &
            allocate (res%modes_k(3*n_atoms, 3*n_atoms, n_kpts), source=(0d0, 0d0))
        if (geom%calc%get_rpa_orders) allocate ( &
            res%rpa_orders_k(geom%calc%param%rpa_order_max, n_kpts), source=0d0 &
        )
        do i_kpt = 1, n_kpts
            k_point = geom%k_pts(:, i_kpt)
            if (.not. geom%calc%do_rpa) then
                res_k = get_mbd_hamiltonian_energy( &
                    geom, alpha_0, C6, damp, dene_k, grad, k_point &
                )
                if (geom%calc%get_eigs) res%mode_eigs_k(:, i_kpt) = res_k%mode_eigs
                if (geom%calc%get_modes) res%modes_k(:, :, i_kpt) = res_k%modes_k_single
            else
                res_k = get_mbd_rpa_energy(geom, alpha, damp, k_point)
                if (geom%calc%get_rpa_orders) then
                    res%rpa_orders_k(:, i_kpt) = res_k%rpa_orders
                end if
            end if
            if (geom%has_exc()) return
            res%energy = res%energy + res_k%energy
        end do ! k_point loop
        res%energy = res%energy/size(geom%k_pts, 2)
        if (geom%calc%get_rpa_orders) res%rpa_orders = res%rpa_orders/n_kpts
    end if
end function mbd_energy

function alpha_dynamic_ts(calc, alpha_0, C6, dalpha, grad) result(alpha)
    type(calc_t), intent(in) :: calc
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(grad_t), allocatable, intent(out) :: dalpha(:)
    type(grad_request_t), intent(in) :: grad
    real(dp) :: alpha(size(alpha_0), 0:ubound(calc%omega_grid, 1))

    integer :: i_freq, n_atoms
    real(dp), allocatable :: omega(:)
    type(grad_t) :: domega

    n_atoms = size(alpha_0)
    omega = omega_eff(C6, alpha_0, domega, grad)
    allocate (dalpha(0:ubound(alpha, 2)))
    do i_freq = 0, ubound(alpha, 2)
        alpha(:, i_freq) = alpha_osc(&
            alpha_0, omega, calc%omega_grid(i_freq), dalpha(i_freq), &
            grad_request_t(dalpha=grad%dalpha, domega=grad%dalpha .or. grad%dC6) &
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
    type(grad_t), intent(out), optional :: dalpha
    type(grad_request_t), intent(in), optional :: grad
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
    type(grad_t), intent(out), optional :: dX
    type(grad_request_t), intent(in), optional :: grad
    real(dp) :: X(size(X_free))

    X = X_free*(V/V_free)**q
    if (.not. present(grad)) return
    if (grad%dX_free) dX%dX_free = X/X_free
    if (grad%dV) dX%dV = X*q/V
    if (grad%dV_free) dX%dV_free = -X*q/V_free
end function

!> \f[
!> \bar C_6=\frac3\pi\int_0^\infty\mathrm du\,\bar\alpha(u)^2,\qquad
!> \partial\bar C_6=\frac6\pi\int_0^\infty\mathrm du
!> \bar\alpha(u)\partial\bar\alpha(u)
!> \f]
function C6_from_alpha(calc, alpha, dC6_dalpha, grad) result(C6)
    type(calc_t), intent(in) :: calc
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

end module

#endif
