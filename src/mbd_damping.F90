! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "defaults.h"

module mbd_damping
!! Damping functions.

use mbd_constants
use mbd_gradients, only: grad_scalar_t, grad_request_t
use mbd_utils, only: lower, exception_t

implicit none

private
public :: damping_fermi, damping_sqrtfermi, op1minus_grad

type, public :: damping_t
    !! Represents a damping function.
    character(len=20) :: version
    real(dp) :: beta = 0d0
    real(dp) :: a = MBD_DAMPING_A
    real(dp) :: ts_d = TS_DAMPING_D
    real(dp) :: ts_sr = 0d0
    real(dp) :: mayer_scaling = 1d0
    real(dp), allocatable :: r_vdw(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: damping_custom(:, :)
    real(dp), allocatable :: potential_custom(:, :, :, :)
    contains
    procedure :: set_params_from_xc => damping_set_params_from_xc
end type

contains

real(dp) function damping_fermi(r, s_vdw, d, df, grad) result(f)
    !! $$
    !! \begin{gathered}
    !! f_{(ij)}=\frac1{1+\exp\big({-}a(\eta-1)\big)},\qquad
    !! \eta=\frac{R_{(ij)}}{S_{\text{vdW}(ij)}}\equiv
    !! \frac{R_{(ij)}}{\beta R_{\text{vdW}(ij)}}
    !! \\ \frac{\mathrm df}{\mathrm dR_c}=
    !! \frac a{2+2\cosh\big(a(\eta-1)\big)}\frac{\mathrm d\eta}{\mathrm dR_c},\qquad
    !! \frac{\mathrm d\eta}{\mathrm dR_c}=
    !! \frac{R_c}{RS_\text{vdW}}-
    !! \frac{R}{S_\text{vdW}^2}\frac{\mathrm dS_\text{vdW}}{\mathrm dR_c}
    !! \end{gathered}
    !! $$
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: s_vdw
    real(dp), intent(in) :: d
    type(grad_scalar_t), intent(out), optional :: df
    type(grad_request_t), intent(in), optional :: grad

    real(dp) :: pre, eta, r_1

    r_1 = sqrt(sum(r**2))
    eta = r_1/s_vdw
    f = 1d0/(1+exp(-d*(eta-1)))
    if (.not. present(grad)) return
    pre = d/(2+2*cosh(d-d*eta))
    if (grad%dcoords) df%dr = pre*r/(r_1*s_vdw)
    if (grad%dr_vdw) df%dvdw = -pre*r_1/s_vdw**2
end function

real(dp) function damping_sqrtfermi(r, s_vdw, d) result(f)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: s_vdw
    real(dp), intent(in) :: d

    f = sqrt(damping_fermi(r, s_vdw, d))
end function

subroutine op1minus_grad(f, df)
    real(dp), intent(inout) :: f
    type(grad_scalar_t), intent(inout) :: df

    f = 1-f
    if (allocated(df%dr)) df%dr = -df%dr
    if (allocated(df%dvdw)) df%dvdw = -df%dvdw
end subroutine

type(exception_t) function damping_set_params_from_xc(this, xc, variant) result(exc)
    class(damping_t), intent(inout) :: this
    character(len=*), intent(in) :: xc
    character(len=*), intent(in) :: variant

    select case (lower(variant))
    case ('ts')
        select case (lower(xc))
        case ('pbe')
            this%ts_sr = 0.94d0
        case ('pbe0')
            this%ts_sr = 0.96d0
        case ('hse')
            this%ts_sr = 0.96d0
        case ('blyp')
            this%ts_sr = 0.62d0
        case ('b3lyp')
            this%ts_sr = 0.84d0
        case ('revpbe')
            this%ts_sr = 0.60d0
        case ('am05')
            this%ts_sr = 0.84d0
        case default
            exc%code = MBD_EXC_DAMPING
            exc%msg = 'Damping parameter S_r of method TS unknown for ' // trim(xc)
        end select
    case ('mbd-rsscs')
        select case (lower(xc))
        case ('pbe')
            this%beta = 0.83d0
        case ('pbe0')
            this%beta = 0.85d0
        case ('hse')
            this%beta = 0.85d0
        case default
            exc%code = MBD_EXC_DAMPING
            exc%msg = 'Damping parameter beta of method MBD@rsSCS unknown for ' // trim(xc)
        end select
    case ('mbd-ts')
        select case (lower(xc))
        case ('pbe')
            this%beta = 0.81d0
        case ('pbe0')
            this%beta = 0.83d0
        case ('hse')
            this%beta = 0.83d0
        case default
            exc%code = MBD_EXC_DAMPING
            exc%msg = 'Damping parameter beta of method MBD@TS unknown for ' // trim(xc)
        end select
    case ('mbd-scs')
        this%beta = 1d0
        select case (lower(xc))
        case ('pbe')
            this%a = 2.56d0
        case ('pbe0')
            this%a = 2.53d0
        case ('hse')
            this%a = 2.53d0
        case default
            exc%code = MBD_EXC_DAMPING
            exc%msg = 'Damping parameter a of method MBD@SCS unknown for ' // trim(xc)
        end select
    case default
        exc%code = MBD_EXC_DAMPING
        exc%msg = 'Method is unkonwn: ' // trim(variant)
    end select
end function

end module
