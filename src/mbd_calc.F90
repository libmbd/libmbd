! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_calc
!! Representing an MBD calculation.

use mbd_constants
use mbd_utils, only: tostr, exception_t, clock_t, abstract_printer, stdout_printer

implicit none

private

type, public :: calc_t
    !! Represents an MBD calculation.
    type(clock_t) :: clock
    type(exception_t) :: exc
    logical :: muted = .false.
    procedure(abstract_printer), pointer, nopass :: printer
    contains
    procedure :: init => calc_init
    procedure :: print => calc_print
    procedure :: destroy => calc_destroy
end type

contains

subroutine calc_init(this)
    class(calc_t), intent(inout) :: this

    this%printer => stdout_printer
    call this%clock%init(100)
end subroutine

subroutine calc_print(this, msg)
    class(calc_t), intent(in) :: this
    character(len=*), intent(in) :: msg

    if (this%muted) return
    call this%printer(msg)
end subroutine

subroutine calc_destroy(this)
    class(calc_t), intent(inout) :: this

    deallocate (this%clock%timestamps, this%clock%counts)
end subroutine

end module
