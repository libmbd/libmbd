! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_common

use mbd_constants
use mbd_calc, only: calc_t
use mbd_gradients, only: grad_t, grad_request_t
use mbd_utils, only: tostr

implicit none

private
public :: omega_eff, sigma_selfint, scale_ts, alpha_dynamic_ts, C6_from_alpha

contains

!> \f[
!> \omega=\frac{4C_6}{3\alpha_{0}^2},\qquad
!> \partial\omega=\omega\left(
!> \frac{\partial C_6}{C_6}-\frac{2\partial\alpha_0}{\alpha_0}
!> \right)
!> \f]
function omega_eff(C6, alpha, domega, grad) result(omega)
    real(dp), intent(in) :: C6(:)
    real(dp), intent(in) :: alpha(:)
    type(grad_t), intent(out), optional :: domega
    type(grad_request_t), intent(in), optional :: grad
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
!> \bar X=X\left(\frac{\bar\alpha_0}{\alpha_0}\right)^q,\qquad
!> \partial\bar X=\bar X\left(
!> \frac{\partial X}{X}+
!> q\frac{\partial\bar\alpha_0}{\bar\alpha_0}-
!> q\frac{\partial\alpha_0}{\alpha_0}
!> \right)
!> \f]
function scale_ts(X_free, V, V_free, q, dX, grad) result(X)
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

end module
