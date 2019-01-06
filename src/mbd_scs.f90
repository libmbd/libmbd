! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_scs
!! Performing self-consistent screening.

use mbd_constants
use mbd_damping, only: damping_t
use mbd_dipole, only: dipole_matrix
use mbd_formulas, only: sigma_selfint
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_t, grad_matrix_re_t, grad_request_t
use mbd_matrix, only: matrix_re_t, contract_cross_33
use mbd_utils, only: findval

implicit none

private
public :: run_scs

contains

function run_scs(geom, alpha, damp, dalpha_scs, grad) result(alpha_scs)
    !! $$
    !! \begin{gathered}
    !! \bar\alpha_i=\tfrac13\operatorname{Tr}
    !! \big(\textstyle\sum_j\boldsymbol{\bar\alpha}_{ij}\big),\qquad
    !! \boldsymbol{\bar\alpha}=(\boldsymbol\alpha^{-1}+\mathbf T_\text{GG})^{-1}
    !! \\ \partial\boldsymbol{\bar\alpha}=
    !! -\boldsymbol{\bar\alpha}(
    !! \partial\boldsymbol\alpha^{-1}+\partial\mathbf T_\text{GG}
    !! )\boldsymbol{\bar\alpha},\qquad
    !! \frac{\partial\bar\alpha_i}{\partial X_j}=
    !! -\frac13\sum_{\zeta\eta}\big(
    !! B_{i\zeta,j\eta}\bar\alpha'_{\zeta,j\eta}+
    !! B'_{j\eta,\zeta}\bar\alpha_{j\eta,i\zeta}
    !! \big)
    !! \\ \mathbf B=\boldsymbol{\bar\alpha}\mathbf A,
    !! \quad A_{i\zeta,j\eta}=
    !! \frac{\partial(\alpha_i^{-1})}{\partial X_i}
    !! \delta_{ij}\delta_{\zeta\eta}+
    !! \frac{\partial T^\text{GG}_{i\zeta,j\eta}}{\partial X_i},\quad
    !! \bar\alpha'_{\zeta,p}=\sum_i\bar\alpha_{i\zeta,p},\quad
    !! B'_{p,\zeta}=\sum_iB_{p,i\zeta}
    !! \end{gathered}
    !! $$
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha(:)
    type(damping_t), intent(in) :: damp
    type(grad_t), intent(out) :: dalpha_scs(:)
    type(grad_request_t), intent(in) :: grad
    real(dp) :: alpha_scs(size(alpha))

    type(matrix_re_t) :: alpha_full, dQ, T
    integer :: n_atoms, i_xyz, i_atom, my_i_atom, i_latt
    type(damping_t) :: damp_local
    real(dp), allocatable :: dsij_dsi(:), dsigma_dalpha(:), &
        alpha_prime(:, :), B_prime(:, :), grads_i(:), dalphadA(:)
    type(grad_matrix_re_t) :: dT
    type(grad_request_t) :: grad_req

    n_atoms = geom%siz()
    damp_local = damp
    damp_local%sigma = sigma_selfint(alpha, dsigma_dalpha, grad%dalpha)
    grad_req = grad_request_t( &
        dcoords=grad%dcoords, &
        dlattice=grad%dlattice, &
        dsigma=grad%dalpha, &
        dr_vdw=grad%dr_vdw &
    )
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
    if (grad%dlattice) then
        do my_i_atom = 1, size(geom%idx%i_atom)
            allocate (dalpha_scs(my_i_atom)%dlattice(3, 3))
        end do
        do i_latt = 1, 3
            do i_xyz = 1, 3
                dQ%val = -dT%dlattice(:, :, i_latt, i_xyz)
                dQ = alpha_full%mmul(dQ)
                dQ = dQ%mmul(alpha_full)
                dalphadA = dQ%contract_n33diag_cols()
                forall (my_i_atom = 1:size(geom%idx%i_atom))
                    dalpha_scs(my_i_atom)%dlattice(i_latt, i_xyz) &
                        = dalphadA(geom%idx%i_atom(my_i_atom))
                end forall
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
end function

end module
