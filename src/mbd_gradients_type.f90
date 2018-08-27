! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_gradients_type

use mbd_constants

implicit none

private
public :: mbd_gradients, mbd_grad_switch

type :: mbd_gradients
    real(dp), allocatable :: dcoords(:, :)  ! n_atoms by 3
    real(dp), allocatable :: dalpha(:)
    real(dp), allocatable :: dalpha_dyn(:, :)  ! n_atoms by 0:n_freq
    real(dp), allocatable :: dC6(:)
    real(dp), allocatable :: dr_vdw(:)
    real(dp), allocatable :: domega(:)
    real(dp), allocatable :: dV(:)
    real(dp), allocatable :: dV_free(:)
    real(dp), allocatable :: dX_free(:)
end type

type :: mbd_grad_switch
    logical :: dcoords = .false.
    logical :: dalpha = .false.
    logical :: dalpha_dyn = .false.
    logical :: dC6 = .false.
    logical :: dr_vdw = .false.
    logical :: domega = .false.
    logical :: dV = .false.
    logical :: dV_free = .false.
    logical :: dX_free = .false.
    contains
    procedure :: any => mbd_grad_switch_any
end type

contains

logical function mbd_grad_switch_any(this) result(any)
    class(mbd_grad_switch), intent(in) :: this

    any = this%dcoords .or. this%dalpha .or. this%dC6 .or. &
        this%dr_vdw .or. this%domega
end function

end module
