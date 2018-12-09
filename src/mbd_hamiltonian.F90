! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MBD_TYPE
module mbd_hamiltonian

use mbd_constants
use mbd_common, only: omega_eff
use mbd_damping, only: damping_t
use mbd_dipole, only: dipole_matrix
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_t, grad_matrix_re_t, grad_matrix_cplx_t, grad_request_t
use mbd_matrix, only: matrix_re_t, matrix_cplx_t
use mbd_utils, only: result_t, tostr

implicit none

#   ifndef MODULE_UNIT_TESTS
private
public :: get_mbd_hamiltonian_energy
#   endif

interface get_mbd_hamiltonian_energy
    module procedure get_mbd_hamiltonian_energy_real
    module procedure get_mbd_hamiltonian_energy_complex
end interface

contains

#   define MBD_TYPE 0
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
type(result_t) function get_mbd_hamiltonian_energy_real( &
        geom, alpha_0, C6, damp, dene, grad) result(res)
#elif MBD_TYPE == 1
type(result_t) function get_mbd_hamiltonian_energy_complex( &
        geom, alpha_0, C6, damp, dene, grad, k_point) result(res)
#endif
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    type(grad_t), intent(out) :: dene
    type(grad_request_t), intent(in) :: grad
#if MBD_TYPE == 1
    real(dp), intent(in) :: k_point(3)
#endif

#if MBD_TYPE == 0
    type(matrix_re_t) :: relay, dQ, T, modes, c_lambda12i_c
    type(grad_matrix_re_t) :: dT
    integer :: i_xyz
#elif MBD_TYPE == 1
    type(matrix_cplx_t) :: relay, T, modes
    type(grad_matrix_cplx_t) :: dT
#endif
    real(dp), allocatable :: eigs(:), omega(:)
    type(grad_t) :: domega
    integer :: n_negative_eigs, n_atoms
    character(120) :: msg

    n_atoms = geom%siz()
#if MBD_TYPE == 0
    T = dipole_matrix(geom, damp, dT, grad)
#elif MBD_TYPE == 1
    T = dipole_matrix(geom, damp, dT, grad, k_point)
#endif
    if (geom%has_exc()) return
    if (grad%any()) then
        call relay%copy_from(T)
    else
        call relay%move_from(T)
    end if
    omega = omega_eff(C6, alpha_0, domega, grad)
    call relay%mult_cross(omega*sqrt(alpha_0))
    call relay%add_diag(omega**2)
    call geom%clock(21)
    if (geom%calc%get_modes .or. grad%any()) then
        call modes%alloc_from(relay)
        allocate (eigs(3*n_atoms))
        call modes%eigh(eigs, geom%calc%exc, src=relay)
        if (geom%calc%get_modes) then
#if MBD_TYPE == 0
            call move_alloc(modes%val, res%modes)
#elif MBD_TYPE == 1
            call move_alloc(modes%val, res%modes_k_single)
#endif
        end if
    else
        eigs = relay%eigvalsh(geom%calc%exc, destroy=.true.)
    end if
    if (geom%has_exc()) return
    call geom%clock(-21)
    if (geom%calc%get_eigs) res%mode_eigs = eigs
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        msg = "CDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues"
        if (geom%calc%param%zero_negative_eigs) then
            where (eigs < 0) eigs = 0d0
            geom%calc%info%neg_eigvals = msg
        else
            geom%calc%exc%code = MBD_EXC_NEG_EIGVALS
            geom%calc%exc%msg = msg
            return
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
#   undef MBD_TYPE
#   define MBD_TYPE 1
#   include "mbd_hamiltonian.F90"

end module

#endif
