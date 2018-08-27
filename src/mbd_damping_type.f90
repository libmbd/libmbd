! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_damping_type

use mbd_constants
use mbd_common, only: lower
use mbd_gradients_type, only: mbd_grad_scalar, mbd_grad_switch

implicit none

private
public :: mbd_damping, damping_fermi, set_damping_parameters

type :: mbd_damping
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
end type mbd_damping

contains

!> \f[
!> f_{(ij)}=\frac1{1+\exp\big({-}a(\eta-1)\big)},\qquad
!> \eta=\frac{R_{(ij)}}{S_{\text{vdW}(ij)}}\equiv
!> \frac{R_{(ij)}}{\beta R_{\text{vdW}(ij)}}
!> \f]
!>
!> \f[
!> \frac{\mathrm df}{\mathrm dR_c}=
!> \frac a{2+2\cosh\big(a(\eta-1)\big)}\frac{\mathrm d\eta}{\mathrm dR_c},\qquad
!> \frac{\mathrm d\eta}{\mathrm dR_c}=
!> \frac{R_c}{RS_\text{vdW}}-
!> \frac{R}{S_\text{vdW}^2}\frac{\mathrm dS_\text{vdW}}{\mathrm dR_c}
!> \f]
real(dp) function damping_fermi(r, s_vdw, d, df, grad) result(f)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: s_vdw
    real(dp), intent(in) :: d
    type(mbd_grad_scalar), intent(out), optional :: df
    type(mbd_grad_switch), intent(in), optional :: grad

    real(dp) :: pre, eta, r_1

    r_1 = sqrt(sum(r**2))
    eta = r_1/s_vdw
    f = 1d0/(1+exp(-d*(eta-1)))
    if (.not. present(grad)) return
    pre = d/(2+2*cosh(d-d*eta))
    if (grad%dcoords) df%dr = pre*r/(r_1*s_vdw)
    if (grad%dr_vdw) df%dvdw = -pre*r_1/s_vdw**2
end function

subroutine set_damping_parameters(xc, ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, &
        mbd_ts_erf_beta, mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta)
    character(len=*), intent(in) :: xc
    real(dp), intent(out) :: &
        ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, mbd_ts_erf_beta, &
        mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta

    ts_d = 20d0
    ts_s_r = 1d0
    mbd_scs_a = 2d0
    mbd_ts_a = 6d0
    mbd_ts_erf_beta = 1d0
    mbd_ts_fermi_beta = 1d0
    mbd_rsscs_a = 6d0
    mbd_rsscs_beta = 1d0
    select case (lower(xc))
        case ("pbe")
            ts_s_r = 0.94d0
            mbd_scs_a = 2.56d0
            mbd_ts_erf_beta = 1.07d0
            mbd_ts_fermi_beta = 0.81d0
            mbd_rsscs_beta = 0.83d0
        case ("pbe0")
            ts_s_r = 0.96d0
            mbd_scs_a = 2.53d0
            mbd_ts_erf_beta = 1.08d0
            mbd_ts_fermi_beta = 0.83d0
            mbd_rsscs_beta = 0.85d0
        case ("hse")
            ts_s_r = 0.96d0
            mbd_scs_a = 2.53d0
            mbd_ts_erf_beta = 1.08d0
            mbd_ts_fermi_beta = 0.83d0
            mbd_rsscs_beta = 0.85d0
        case ("blyp")
            ts_s_r = 0.62d0
        case ("b3lyp")
            ts_s_r = 0.84d0
        case ("revpbe")
            ts_s_r = 0.60d0
        case ("am05")
            ts_s_r = 0.84d0
    endselect
end subroutine set_damping_parameters

end module
