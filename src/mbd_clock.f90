! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_clock

implicit none

private
public :: clock_t

integer, parameter :: n_timestamps = 100

type :: clock_t
    logical :: active = .true.
    integer :: timestamps(n_timestamps), counts(n_timestamps)
    integer, private :: cnt, rate, cnt_max, aid
    contains
    procedure :: clock => clock_clock
end type clock_t

contains

subroutine clock_clock(this, id)
    class(clock_t), intent(inout) :: this
    integer, intent(in) :: id

    if (.not. this%active) return
    call system_clock(this%cnt, this%rate, this%cnt_max)
    if (id > 0) then
        this%timestamps(id) = this%timestamps(id)-this%cnt
    else
        this%aid = abs(id)
        this%timestamps(this%aid) = this%timestamps(this%aid)+this%cnt
        this%counts(this%aid) = this%counts(this%aid)+1
    end if
end subroutine

end module
