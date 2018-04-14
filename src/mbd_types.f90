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
    real(dp), allocatable :: re_dvdw(:, :)
    real(dp), allocatable :: re_dsigma(:, :)
    contains
    procedure  :: siz => mat3n3n_siz
end type

type :: mat33
    real(dp) :: val(3, 3)
    ! explicit derivative, [abc] ~ dval_{ab}/dR_c
    real(dp), allocatable :: dr(:, :, :)
    real(dp), allocatable :: dvdw(:, :)
    real(dp), allocatable :: dsigma(:, :)
end type

type :: vecn
    real(dp), allocatable :: val(:)
    real(dp), allocatable :: dr(:, :, :)
    contains
    procedure  :: siz => vecn_siz
end type

interface vecn
    module procedure vecn_constructor_no_dr
end interface

type :: scalar
    real(dp) :: val
    real(dp), allocatable :: dr(:)  ! explicit derivative
    real(dp), allocatable :: dvdw
end type

contains

type(vecn) function vecn_constructor_no_dr(x) result(vec)
    real(dp), intent(in) :: x(:)

    vec%val = x
end function

integer function mat3n3n_siz(this, ndim)
    class(mat3n3n), intent(in) :: this
    integer, intent(in) :: ndim

    if (allocated(this%re)) then
        mat3n3n_siz = size(this%re, ndim)
    elseif (allocated(this%cplx)) then
        mat3n3n_siz = size(this%cplx, ndim)
    else
        mat3n3n_siz = 0
    end if
end function

integer function vecn_siz(this)
    class(vecn), intent(in) :: this

    if (allocated(this%val)) then
        vecn_siz = size(this%val)
    else
        vecn_siz = 0
    end if
end function

end module
