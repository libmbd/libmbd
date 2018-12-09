! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_common

use mbd_constants
use mbd_gradients, only: grad_t, grad_request_t

implicit none

private
public :: omega_eff, sigma_selfint

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

end module
