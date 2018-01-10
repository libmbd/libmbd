! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_types

use mbd_common, only: dp

implicit none

type :: mat3n3n
    real(dp), allocatable :: re(:, :)
    complex(dp), allocatable :: cplx(:, :)
    real(dp), allocatable :: re_dr(:, :, :)
end type

type :: mat33
    real(dp) :: val(3, 3)
    ! explicit derivative, [abc] ~ dval_{ab}/dR_c
    real(dp) :: dr(3, 3, 3)
    logical :: has_vdw = .false.
    real(dp) :: dvdw(3, 3)
    logical :: has_sigma = .false.
    real(dp) :: dsigma(3, 3)
end type

type :: scalar
    real(dp) :: val
    real(dp) :: dr(3)  ! explicit derivative
    real(dp) :: dvdw
end type

end module
