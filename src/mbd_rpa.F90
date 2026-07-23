! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DO_COMPLEX_TYPE
module mbd_rpa

use mbd_constants
use mbd_damping, only: damping_t
use mbd_dipole, only: dipole_matrix
use mbd_formulas, only: sigma_selfint, rpa_rescale_eigval
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_t, grad_matrix_re_t, grad_matrix_cplx_t, &
    grad_request_t
use mbd_matrix, only: matrix_re_t, matrix_cplx_t
use mbd_utils, only: result_t, tostr

implicit none

private
public :: get_mbd_rpa_energy

interface get_mbd_rpa_energy
    !! Get the MBD energy by the RPA method.
    !!
    !! The energy is obtained from the eigenvalues \(\mu_k(\mathrm iu)\) of the
    !! symmetrized dipole matrix \(\mathbf M(\mathrm iu)=\mathbf A(\mathrm
    !! iu)^\frac12\mathbf T\mathbf A(\mathrm iu)^\frac12\), with \(A_{ij}(\mathrm
    !! iu)=\alpha_i(\mathrm iu)\delta_{ij}\),
    !!
    !! $$
    !! E=\frac1{2\pi}\int_0^\infty\mathrm du\sum_k\log(1+\mu_k(\mathrm iu))
    !! =\frac1{2\pi}\int_0^\infty\mathrm du\,\log\det(\mathbf 1+\mathbf M(\mathrm iu))
    !! $$
    !!
    !! The derivatives use \(\partial\mu_k=\mathbf c_k^\dagger(\partial\mathbf
    !! M)\mathbf c_k\), which upon summation over the eigenvalues gives the
    !! weighted resolvent \(\mathbf B(\mathrm iu)=\mathbf C\operatorname{diag}(
    !! g'(\mu_k))\mathbf C^\dagger\), where \(g(\mu)=\log(1+\mu)\) is the
    !! per-eigenvalue energy contribution,
    !!
    !! $$
    !! \partial E=\frac1{2\pi}\int_0^\infty\mathrm du
    !! \operatorname{Tr}\big(\mathbf B(\mathrm iu)\partial\mathbf M(\mathrm iu)\big)
    !! $$
    !!
    !! Both the explicit derivatives (coordinates, lattice vectors, van der
    !! Waals radii, \(\mathbf q\)) and the implicit derivatives with respect to
    !! the dynamic polarizabilities are returned. The latter propagate both
    !! through the \(\sqrt{\alpha_i\alpha_j}\) prefactor and, when the damping
    !! depends on it, through \(\sigma_{ij}(\alpha)\). With eigenvalue rescaling
    !! ([[mbd_calc_params_t:rpa_rescale_eigs]]) the contribution becomes
    !! \(g(\mu)=\log(1+\lambda(\mu))-\lambda(\mu)\), with
    !! \(\lambda=-\operatorname{erf}(\tfrac{\sqrt\pi}2\mu^4)^\frac14\) for
    !! \(\mu<0\), and the same weighted resolvent is used.
    !!
    !! The real-typed version is equivalent to \(\mathbf q=0\).
    module procedure get_mbd_rpa_energy_real
    module procedure get_mbd_rpa_energy_complex
end interface

contains

#endif

#ifndef DO_COMPLEX_TYPE
type(result_t) function get_mbd_rpa_energy_real( &
        geom, alpha, damp, grad) result(res)
#else
type(result_t) function get_mbd_rpa_energy_complex( &
        geom, alpha, damp, q, grad) result(res)
#endif
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha(:, 0:)
    type(damping_t), intent(in) :: damp
#ifdef DO_COMPLEX_TYPE
    real(dp), intent(in) :: q(3)
#endif
    type(grad_request_t), intent(in), optional :: grad

#ifndef DO_COMPLEX_TYPE
    type(matrix_re_t) :: relay, AT, Mmat, modes, B, dQ
    type(grad_matrix_re_t) :: dT
#else
    type(matrix_cplx_t) :: relay, AT, Mmat, modes, B, dQ
    type(grad_matrix_cplx_t) :: dT
#endif
    real(dp), allocatable :: eigs(:), log_eigs(:), sqrt_alpha(:), &
        g_prime(:), contr(:), dxr(:)
    integer :: i_freq, my_i_atom, n_order, n_negative_eigs, my_j_atom, &
        n_atoms, i_xyz, i_latt
    real(dp) :: freq_w, sigma_ij
    type(damping_t) :: damp_alpha
    type(grad_request_t) :: grad_dip
    logical :: do_grad

    do_grad = .false.
    if (present(grad)) do_grad = grad%any()
    n_atoms = geom%siz()
    if (do_grad) then
        allocate (g_prime(3 * n_atoms))
        grad_dip%dcoords = grad%dcoords
        grad_dip%dlattice = grad%dlattice
        grad_dip%dr_vdw = grad%dr_vdw
        ! the dipole matrix depends on the dynamic polarizability through the
        ! self-consistent-screening width sigma only for the gg dampings
        grad_dip%dsigma = grad%dalpha_dyn .and. index(damp%version, 'gg') > 0
#ifdef DO_COMPLEX_TYPE
        grad_dip%dq = grad%dq
#endif
        if (grad%dcoords) allocate (res%dE%dcoords(n_atoms, 3), source=0d0)
        if (grad%dr_vdw) allocate (res%dE%dr_vdw(n_atoms), source=0d0)
        if (grad%dalpha_dyn) &
            allocate (res%dE%dalpha_dyn(n_atoms, 0:ubound(alpha, 2)), source=0d0)
#ifndef DO_COMPLEX_TYPE
        if (grad%dlattice) allocate (res%dE%dlattice(3, 3), source=0d0)
#else
        if (grad%dlattice) allocate (res%dE%dlattice(3, 3), source=0d0)
        if (grad%dq) allocate (res%dE%dq(3), source=0d0)
#endif
    end if

    res%energy = 0d0
    damp_alpha = damp
    ! implicit allocation doesn't work here in gfortran 4.9
    allocate (eigs(3 * geom%siz()), log_eigs(3 * geom%siz()))
    if (geom%get_rpa_orders) allocate (res%rpa_orders(geom%param%rpa_order_max), source=0d0)
    do i_freq = 0, ubound(geom%freq, 1)
        damp_alpha%sigma = sigma_selfint(alpha(:, i_freq))
        sqrt_alpha = sqrt(alpha(:, i_freq))
        if (do_grad) then
#ifndef DO_COMPLEX_TYPE
            relay = dipole_matrix(geom, damp_alpha, dT, grad_dip)
#else
            relay = dipole_matrix(geom, damp_alpha, dT, grad_dip, q=q)
#endif
        else
#ifndef DO_COMPLEX_TYPE
            relay = dipole_matrix(geom, damp_alpha)
#else
            relay = dipole_matrix(geom, damp_alpha, q=q)
#endif
        end if
        if (geom%has_exc()) return
        if (do_grad) then
            ! keep relay = T, build M = A^1/2 T A^1/2 separately
            call Mmat%copy_from(relay)
            call Mmat%mult_cross(sqrt_alpha)
            call modes%alloc_from(Mmat)
            call geom%clock(23)
            call modes%eigh(eigs, geom%exc, src=Mmat, clock=geom%timer)
            call geom%clock(-23)
        else
            call relay%mult_cross(sqrt_alpha)
            call AT%move_from(relay)
            call geom%clock(23)
            eigs = AT%eigvalsh(geom%exc, destroy=.true.)
            call geom%clock(-23)
        end if
        if (geom%has_exc()) return
        if (geom%param%rpa_rescale_eigs) then
            ! dxr (dlambda/dmu) is filled only when do_grad, reused below
            eigs = rpa_rescale_eigval(eigs, dxr, grad=do_grad)
        end if
        n_negative_eigs = count(eigs(:) <= -1)
        if (n_negative_eigs > 0) then
            geom%exc%code = MBD_EXC_NEG_EIGVALS
            geom%exc%msg = "1+AT matrix has "// &
                trim(tostr(n_negative_eigs))//" negative eigenvalues"
            return
        end if
        log_eigs = log(1 + eigs)
        if (geom%param%rpa_rescale_eigs) then
            log_eigs = log_eigs - eigs
        end if
        res%energy = res%energy + &
            1d0 / (2 * pi) * sum(log_eigs) * geom%freq(i_freq)%weight
        if (geom%get_rpa_orders) then
            do n_order = 2, geom%param%rpa_order_max
                res%rpa_orders(n_order) = res%rpa_orders(n_order) &
                    + (-1d0 / (2 * pi) * (-1)**n_order &
                    * sum(eigs**n_order) / n_order) &
                    * geom%freq(i_freq)%weight
            end do
        end if
        if (.not. do_grad) cycle
        freq_w = geom%freq(i_freq)%weight
        ! The per-eigenvalue energy contribution is g(mu_k), whose derivative
        ! with respect to a parameter is g'(mu_k) times the derivative of the
        ! raw eigenvalue mu_k of M. Summed, this gives the weighted resolvent
        ! B = C diag(g'(mu_k)) C^dagger, contracted below with dM.
        if (.not. geom%param%rpa_rescale_eigs) then
            ! g(mu) = log(1 + mu)
            g_prime = 1d0 / (1d0 + eigs)
        else
            ! g(mu) = log(1 + lambda) - lambda with lambda = eigs (rescaled) and
            ! dlambda/dmu = dxr; chain rule g'(mu) = (1/(1 + lambda) - 1) dlambda/dmu
            g_prime = (1d0 / (1d0 + eigs) - 1d0) * dxr
        end if
        call B%copy_from(modes)
        call B%mult_cols_3n(g_prime)
        B = B%mmul(modes, transB='C')
#ifdef DO_COMPLEX_TYPE
        B%val = conjg(B%val)
#endif
        call dQ%init_from(relay)
        if (grad%dcoords) then
            do i_xyz = 1, 3
                dQ%val = dT%dr(:, :, i_xyz)
                call dQ%mult_cross(sqrt_alpha)
                dQ%val = B%val * dQ%val
                contr = freq_w / pi * dble(dQ%contract_n33_rows())
                res%dE%dcoords(:, i_xyz) = res%dE%dcoords(:, i_xyz) + contr
            end do
        end if
        if (grad%dlattice) then
            do i_latt = 1, 3
                do i_xyz = 1, 3
                    dQ%val = dT%dlattice(:, :, i_latt, i_xyz)
                    call dQ%mult_cross(sqrt_alpha)
                    dQ%val = B%val * dQ%val
                    res%dE%dlattice(i_latt, i_xyz) = res%dE%dlattice(i_latt, i_xyz) + &
                        freq_w / (2 * pi) * dble(dQ%sum_all())
                end do
            end do
        end if
        if (grad%dr_vdw) then
            dQ%val = dT%dvdw
            call dQ%mult_cross(sqrt_alpha)
            dQ%val = B%val * dQ%val
            contr = freq_w / pi * dble(dQ%contract_n33_rows())
            res%dE%dr_vdw = res%dE%dr_vdw + contr
        end if
#ifdef DO_COMPLEX_TYPE
        if (grad%dq) then
            do i_latt = 1, 3
                dQ%val = dT%dq(:, :, i_latt)
                call dQ%mult_cross(sqrt_alpha)
                dQ%val = B%val * dQ%val
                res%dE%dq(i_latt) = res%dE%dq(i_latt) + &
                    freq_w / (2 * pi) * dble(dQ%sum_all())
            end do
        end if
#endif
        if (grad%dalpha_dyn) then
            ! channel through the sqrt(alpha_i alpha_j) prefactor
            dQ%val = relay%val
            call dQ%mult_cross(sqrt_alpha)
            call dQ%mult_rows(1d0 / (2 * alpha(:, i_freq)))
            dQ%val = B%val * dQ%val
            contr = freq_w / pi * dble(dQ%contract_n33_rows())
            res%dE%dalpha_dyn(:, i_freq) = res%dE%dalpha_dyn(:, i_freq) + contr
            ! channel through sigma_ij(alpha), if the damping uses it
            if (grad_dip%dsigma) then
                dQ%val = dT%dsigma
                call dQ%mult_cross(sqrt_alpha)
                dQ%val = B%val * dQ%val
                ! scale each block (i, j) by 1 / sigma_ij
                do my_i_atom = 1, size(dQ%idx%i_atom)
                    do my_j_atom = 1, size(dQ%idx%j_atom)
                        associate ( &
                                i_atom => dQ%idx%i_atom(my_i_atom), &
                                j_atom => dQ%idx%j_atom(my_j_atom), &
                                dQ_sub => dQ%val( &
                                    3 * (my_i_atom - 1) + 1:, &
                                    3 * (my_j_atom - 1) + 1: &
                                ) &
                        )
                            sigma_ij = damp%mayer_scaling * sqrt(sum( &
                                damp_alpha%sigma([i_atom, j_atom])**2))
                            dQ_sub(:3, :3) = dQ_sub(:3, :3) / sigma_ij
                        end associate
                    end do
                end do
                contr = freq_w / pi * damp%mayer_scaling**2 &
                    * damp_alpha%sigma**2 / (3 * alpha(:, i_freq)) &
                    * dble(dQ%contract_n33_rows())
                res%dE%dalpha_dyn(:, i_freq) = res%dE%dalpha_dyn(:, i_freq) + contr
            end if
        end if
    end do
end function

#ifndef DO_COMPLEX_TYPE
#   define DO_COMPLEX_TYPE
#   include "mbd_rpa.F90"

end module

#endif
