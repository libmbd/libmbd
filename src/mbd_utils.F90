! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_utils
!! Utility types, interfaces, and procedures.

use mbd_constants
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

private
public :: tostr, diff3, diff5, print_matrix, lower, diff7, findval, shift_idx, &
    is_true, printer

interface tostr
    module procedure tostr_int
    module procedure tostr_real
end interface

interface findval
    module procedure findval_int
end interface

type, public :: exception_t
    !! Represents an exception.
    integer :: code = 0
    character(50) :: origin = '(unknown)'
    character(150) :: msg = ''
end type

type, public :: result_t
    !! Stores results from an MBD calculation
    real(dp) :: energy
    real(dp), allocatable :: mode_eigs(:)
    real(dp), allocatable :: modes(:, :)
    real(dp), allocatable :: rpa_orders(:)
    real(dp), allocatable :: mode_eigs_k(:, :)
    complex(dp), allocatable :: modes_k(:, :, :)
    complex(dp), allocatable :: modes_k_single(:, :)
    real(dp), allocatable :: rpa_orders_k(:, :)
end type

type, public :: atom_index_t
    !! Maps from atom indexes to positions in matrices.
    integer, allocatable :: i_atom(:)
    integer, allocatable :: j_atom(:)
    integer :: n_atoms
#   ifdef WITH_SCALAPACK
    logical :: parallel
#   endif
end type

type, public :: clock_t
    !! Used for measuring performance.
    logical :: active = .true.
    integer(kind=8), allocatable :: timestamps(:), counts(:)
    contains
    procedure :: init => clock_init
    procedure :: clock => clock_clock
    procedure :: print => clock_print
end type

type, public :: quad_pt_t
    !! Represents a 1D quadrature point
    real(dp) :: val
    real(dp) :: weight
end type

contains

character(len=50) elemental function tostr_int(k, format) result(s)
    integer, intent(in) :: k
    character(len=*), intent(in), optional :: format

    if (present(format)) then
        write (s, format) k
    else
        write (s, "(i20)") k
    end if
    s = adjustl(s)
end function

character(len=50) elemental function tostr_real(x, format) result(s)
    real(dp), intent(in) :: x
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (s, format) x
    else
        write (s, "(g50.17e3)") x
    end if
    s = adjustl(s)
end function

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

    diff7 = ( &
        -1.d0/60*x(-3) &
        + 3.d0/20*x(-2) &
        - 3.d0/4*x(-1) &
        + 3.d0/4*x(1) &
        - 3.d0/20*x(2) &
        + 1.d0/60*x(3) &
    )/delta
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

integer pure function findval_int(array, val) result(idx)
    integer, intent(in) :: array(:), val

    integer :: i

    idx = 0
    do i = 1, size(array)
        if (val == array(i)) then
            idx = i
            return
        end if
    end do
end function

subroutine shift_idx(idx, start, finish)
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: start(:), finish(:)

    integer :: i_dim, i

    do i_dim = size(idx), 1, -1
        i = idx(i_dim)+1
        if (i <= finish(i_dim)) then
            idx(i_dim) = i
            return
        else
            idx(i_dim) = start(i_dim)
        end if
    end do
end subroutine

subroutine clock_init(this, n)
    class(clock_t), intent(inout) :: this
    integer, intent(in) :: n

    allocate (this%timestamps(n), source=0_8)
    allocate (this%counts(n), source=0_8)
end subroutine

subroutine clock_clock(this, id)
    class(clock_t), intent(inout) :: this
    integer, intent(in) :: id

    integer(kind=8) :: cnt, rate, cnt_max
    integer :: absid

    if (.not. this%active) return
    call system_clock(cnt, rate, cnt_max)
    if (id > 0) then
        this%timestamps(id) = this%timestamps(id) - cnt
    else
        absid = abs(id)
        this%timestamps(absid) = this%timestamps(absid) + cnt
        this%counts(absid) = this%counts(absid)+1
    end if
end subroutine

subroutine clock_print(this)
    class(clock_t), intent(inout) :: this

    integer(kind=8) :: cnt, rate, cnt_max
    integer :: i
    character(len=20) :: label

#ifdef WITH_MPI
    if (mpi_get_rank() /= 0) return
#endif
    call system_clock(cnt, rate, cnt_max)

    print '(A20,A10,A10)', 'id', 'count', 'time (s)'
    do i = 1, size(this%counts)
        if (this%counts(i) == 0) cycle
        select case (i)
        case (11); label = 'dipole real space'
        case (12); label = 'dipole rec space'
        case (13); label = 'dipole real 3x3'
        case (21); label = 'eig MBD'
        case (23); label = 'eig RPA'
        case (24); label = 'eig RPA orders'
        case (25); label = 'MBD forces'
        case (32); label = 'inv SCS'
        case (50); label = 'SCS'
        case (51); label = 'single k-point'
        case (90); label = 'MBD@rsSCS'
        case (91); label = 'MBD@rsSCS forces'
        case default
            label = '(' // trim(tostr(i)) // ')'
        end select
        print '(A20,I10,F10.3)', label, this%counts(i), dble(this%timestamps(i))/rate
    end do
end subroutine

subroutine printer(str)
    character(len=*) :: str

#ifdef WITH_MPI
    if (mpi_get_rank() /= 0) return
#endif
    print *, str
end subroutine

logical function is_true(val) result(res)
    logical, intent(in), optional :: val

    res = .false.
    if (present(val)) res = val
end function

end module
