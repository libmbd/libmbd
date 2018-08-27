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
    contains
    procedure :: copy_alloc => mbd_gradients_copy_alloc
    procedure :: has_grad => mbd_gradients_has_grad
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
end type

contains

subroutine mbd_gradients_copy_alloc(this, other)
    class(mbd_gradients), intent(in) :: this
    type(mbd_gradients), intent(out) :: other

    if (allocated(this%dcoords)) &
        allocate (other%dcoords(size(this%dcoords, 2), 3))
    if (allocated(this%dalpha)) &
        allocate (other%dalpha(size(this%dalpha)))
    if (allocated(this%dC6)) allocate (other%dC6(size(this%dC6)))
    if (allocated(this%dr_vdw)) allocate (other%dr_vdw(size(this%dr_vdw)))
    if (allocated(this%domega)) allocate (other%domega(size(this%domega)))
end subroutine

logical function mbd_gradients_has_grad(this) result(has_grad)
    class(mbd_gradients), intent(in) :: this

    has_grad = allocated(this%dcoords) .or. &
        allocated(this%dalpha) .or. allocated(this%dC6) .or. &
        allocated(this%dr_vdw) .or. allocated(this%domega)
end function

end module
