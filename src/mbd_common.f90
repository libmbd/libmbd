! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_common

implicit none

private
public :: tostr

interface tostr
    module procedure tostr_int_
    module procedure tostr_dble_
end interface

contains


character(len=50) elemental function tostr_int_(k, format)
    implicit none

    integer, intent(in) :: k
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (tostr_int_, format) k
    else
        write (tostr_int_, "(i20)") k
    end if
    tostr_int_ = adjustl(tostr_int_)
end function tostr_int_


character(len=50) elemental function tostr_dble_(x, format)
    implicit none

    double precision, intent(in) :: x
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (tostr_dble_, format) x
    else
        write (tostr_dble_, "(g50.17e3)") x
    end if
    tostr_dble_ = adjustl(tostr_dble_)
end function tostr_dble_

end module
