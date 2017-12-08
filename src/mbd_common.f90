! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_common

implicit none

private
public :: tostr, diff3, diff5, print_matrix

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


real(8) pure function diff3(x, delta)
    real(8), intent(in) :: x(-1:)
    real(8), intent(in) :: delta

    diff3 = (x(1)-x(-1))/(2*delta)
end function


real(8) pure function diff5(x, delta)
    real(8), intent(in) :: x(-2:)
    real(8), intent(in) :: delta

    diff5 = (1.d0/12*x(-2)-2.d0/3*x(-1)+2.d0/3*x(1)-1.d0/12*x(2))/delta
end function


subroutine print_matrix(label, A)
    character(len=*), intent(in) :: label
    real(8), intent(in) :: A(:, :)

    integer :: m, n, i, j

    m = size(A, 1)
    n = size(A, 2)
    write (6, '(A,":")') label
    do i = 1, m
        do j = 1, n
            write (6, "(g10.3)", advance="no") A(i, j)
        end do
        write (6, *)
    end do
end subroutine

end module
