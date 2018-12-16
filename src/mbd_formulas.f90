! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_formulas
!! Common formulas used at multiple places.

use mbd_constants
use mbd_calc, only: calc_t
use mbd_gradients, only: grad_t, grad_request_t
use mbd_utils, only: tostr

implicit none

private
public :: omega_qho, alpha_dyn_qho, C6_from_alpha, sigma_selfint, scale_with_ratio

contains

function omega_qho(C6, alpha, domega, grad) result(omega)
    !! $$
    !! \omega=\frac{4C_6}{3\alpha_{0}^2},\qquad
    !! \partial\omega=\omega\left(
    !! \frac{\partial C_6}{C_6}-\frac{2\partial\alpha_0}{\alpha_0}
    !! \right)
    !! $$
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

function alpha_dyn_qho(calc, alpha_0, omega, dalpha, grad) result(alpha)
    !! $$
    !! \alpha(\mathrm iu)=\frac{\alpha_0}{1+u^2/\omega^2},\qquad
    !! \partial\alpha(\mathrm iu)=\alpha(\mathrm iu)\left(
    !! \frac{\partial\alpha_0}{\alpha_0}+
    !! \frac2\omega\frac{\partial\omega}{1+\omega^2/u^2}
    !! \right)
    !! $$
    type(calc_t), intent(in) :: calc
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: omega(:)
    type(grad_t), allocatable, intent(out) :: dalpha(:)
    type(grad_request_t), intent(in) :: grad
    real(dp) :: alpha(size(alpha_0), 0:ubound(calc%omega_grid, 1))

    integer :: i_freq, n_atoms

    n_atoms = size(alpha_0)
    allocate (dalpha(0:ubound(alpha, 2)))
    do i_freq = 0, ubound(alpha, 2)
        associate (alpha => alpha(:, i_freq), u => calc%omega_grid(i_freq))
            alpha = alpha_0/(1+(u/omega)**2)
            if (grad%dalpha) dalpha(i_freq)%dalpha = alpha/alpha_0
            if (grad%domega) dalpha(i_freq)%domega = alpha*2d0/omega/(1d0+(omega/u)**2)
        end associate
    end do
end function

function C6_from_alpha(calc, alpha, dC6_dalpha, grad) result(C6)
    !! $$
    !! \bar C_6=\frac3\pi\int_0^\infty\mathrm du\,\bar\alpha(u)^2,\qquad
    !! \partial\bar C_6=\frac6\pi\int_0^\infty\mathrm du
    !! \bar\alpha(u)\partial\bar\alpha(u)
    !! $$
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

function sigma_selfint(alpha, dsigma_dalpha, grad) result(sigma)
    !! $$
    !! \begin{gathered}
    !! \sigma_i(u)=\left(\frac13\sqrt{\frac2\pi}\alpha_i(u)\right)^{\frac13},\qquad
    !! \partial\sigma_i=\sigma_i\frac{\partial\alpha_i}{3\alpha_i}
    !! \\ \sigma_{ij}(u)=\sqrt{\sigma_i(u)^2+\sigma_j(u)^2},\qquad
    !! \partial\sigma_{ij}=
    !! \frac{\sigma_i\partial\sigma_i+\sigma_j\partial\sigma_j}{\sigma_{ij}}
    !! \end{gathered}
    !! $$
    real(dp), intent(in) :: alpha(:)
    real(dp), allocatable, intent(out), optional :: dsigma_dalpha(:)
    logical, intent(in), optional :: grad
    real(dp) :: sigma(size(alpha))

    sigma = (sqrt(2d0/pi)*alpha/3d0)**(1d0/3)
    if (.not. present(grad)) return
    if (grad) dsigma_dalpha = sigma/(3*alpha)
end function

function scale_with_ratio(x, yp, y, q, dx, grad) result(xp)
    !! $$
    !! x'=x\left(\frac{y'}y\right)^q,\qquad
    !! \partial x'=x\left(
    !! \frac{\partial x}x+
    !! q\frac{\partial y'}{y'}-
    !! q\frac{\partial y}{y}
    !! \right)
    !! $$
    real(dp), intent(in) :: x(:), yp(:), y(:)
    real(dp), intent(in) :: q
    type(grad_t), intent(out), optional :: dx
    type(grad_request_t), intent(in), optional :: grad
    real(dp) :: xp(size(x))

    xp = x*(yp/y)**q
    if (.not. present(grad)) return
    if (grad%dX_free) dx%dX_free = xp/x
    if (grad%dV) dx%dV = xp*q/yp
    if (grad%dV_free) dx%dV_free = -xp*q/y
end function

end module
