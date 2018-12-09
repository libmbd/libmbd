! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_utils

use mbd_constants

implicit none

private
public :: tostr, diff3, diff5, print_matrix, lower, exception_t, diff7, &
    findval, printer, shift_cell, result_t, clock_t

interface tostr
    module procedure tostr_int_
    module procedure tostr_dble_
end interface

type :: exception_t
    integer :: code = 0
    character(50) :: origin = '(unknown)'
    character(150) :: msg = ''
end type

type :: result_t
    real(dp) :: energy
    real(dp), allocatable :: mode_eigs(:)
    real(dp), allocatable :: modes(:, :)
    real(dp), allocatable :: rpa_orders(:)
    real(dp), allocatable :: mode_eigs_k(:, :)
    complex(dp), allocatable :: modes_k(:, :, :)
    complex(dp), allocatable :: modes_k_single(:, :)
    real(dp), allocatable :: rpa_orders_k(:, :)
end type

integer, parameter :: n_timestamps = 100

type :: clock_t
    logical :: active = .true.
    integer :: timestamps(n_timestamps), counts(n_timestamps)
    integer, private :: cnt, rate, cnt_max, aid
    contains
    procedure :: clock => clock_clock
end type clock_t

abstract interface
    subroutine printer(msg)
        character(len=*), intent(in) :: msg
    end subroutine
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

real(dp) pure function diff3(x, delta)
    real(dp), intent(in) :: x(-1:)
    real(dp), intent(in) :: delta

    diff3 = (x(1)-x(-1))/(2*delta)
end function

real(dp) pure function diff5(x, delta)
    real(dp), intent(in) :: x(-2:)
    real(dp), intent(in) :: delta

    diff5 = (1.d0/12*x(-2)-2.d0/3*x(-1)+2.d0/3*x(1)-1.d0/12*x(2))/delta
end function

real(dp) pure function diff7(x, delta)
    real(dp), intent(in) :: x(-3:)
    real(dp), intent(in) :: delta

    diff7 = (-1.d0/60*x(-3)+3.d0/20*x(-2)-3.d0/4*x(-1)+3.d0/4*x(1)-3.d0/20*x(2)+1.d0/60*x(3))/delta
end function

subroutine print_matrix(label, A, prec)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: A(:, :)
    integer, optional, intent(in) :: prec

    integer :: m, n, i, j, prec_
    character(len=10) :: fm

    if (present(prec)) then
        prec_ = prec
    else
        prec_ = 3
    end if
    m = size(A, 1)
    n = size(A, 2)
    write (fm, '("(g",i2,".",i1,")")') prec_+8, prec_
    write (6, '(A,":")') label
    do i = 1, m
        do j = 1, n
            write (6, fm, advance="no") A(i, j)
        end do
        write (6, *)
    end do
end subroutine

pure function lower(str)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lower

    integer :: i

    do i = 1, len(str)
        select case (str(i:i))
            case ('A':'Z')
                lower(i:i) = achar(iachar(str(i:i))+32)
            case default
                lower(i:i) = str(i:i)
        end select
    end do
end function

integer pure function findval(array, val)
    integer, intent(in) :: array(:), val

    integer :: i

    findval = 0
    do i = 1, size(array)
        if (val == array(i)) then
            findval = i
            return
        end if
    end do
end function

subroutine shift_cell(ijk, first_cell, last_cell)
    integer, intent(inout) :: ijk(3)
    integer, intent(in) :: first_cell(3), last_cell(3)

    integer :: i_dim, i

    do i_dim = 3, 1, -1
        i = ijk(i_dim)+1
        if (i <= last_cell(i_dim)) then
            ijk(i_dim) = i
            return
        else
            ijk(i_dim) = first_cell(i_dim)
        end if
    end do
end subroutine

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
