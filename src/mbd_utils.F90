! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_utils
!! Utility types, interfaces, and procedures.

use mbd_constants
use mbd_gradients, only: grad_t
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

private
public :: tostr, diff3, diff5, lower, diff7, findval, shift_idx, &
    is_true, printer_i, printer

integer, parameter :: i8 = selected_int_kind(18)

interface tostr
    module procedure tostr_int
    module procedure tostr_real
end interface

interface findval
    module procedure findval_int
end interface

abstract interface
    subroutine printer_i(str)
        character(len=*), intent(in) :: str
    end subroutine
end interface

type, public :: logger_t
    integer :: level = MBD_LOG_LVL_WARN
    ! TODO cannot use printer() as default because of PGI 19.4
    procedure(printer_i), nopass, pointer :: printer => null()
    contains
    procedure :: info => logger_info
    procedure :: debug => logger_debug
    procedure :: warn => logger_warn
    procedure :: error => logger_error
end type

type, public :: exception_t
    !! Represents an exception.
    integer :: code = 0
    character(50) :: origin = '(unknown)'
    character(150) :: msg = ''
end type

type, public :: result_t
    !! Stores results from an MBD calculation
    real(dp) :: energy
    type(grad_t) :: dE
    real(dp), allocatable :: mode_eigs(:)
    real(dp), allocatable :: modes(:, :)
    real(dp), allocatable :: rpa_orders(:)
    real(dp), allocatable :: mode_eigs_k(:, :)
    complex(dp), allocatable :: modes_k(:, :, :)
    complex(dp), allocatable :: modes_k_single(:, :)
    real(dp), allocatable :: rpa_orders_k(:, :)
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
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
    integer :: level = 0
    integer(i8), allocatable :: timestamps(:), counts(:)
    integer, allocatable :: levels(:)
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

    diff3 = (x(1) - x(-1)) / (2 * delta)
end function

real(dp) pure function diff5(x, delta)
    real(dp), intent(in) :: x(-2:)
    real(dp), intent(in) :: delta

    diff5 = (1.d0 / 12 * x(-2) - 2.d0 / 3 * x(-1) + 2.d0 / 3 * x(1) - 1.d0 / 12 * x(2)) / delta
end function

real(dp) pure function diff7(x, delta)
    real(dp), intent(in) :: x(-3:)
    real(dp), intent(in) :: delta

    diff7 = ( &
        -1.d0 / 60 * x(-3) &
        + 3.d0 / 20 * x(-2) &
        - 3.d0 / 4 * x(-1) &
        + 3.d0 / 4 * x(1) &
        - 3.d0 / 20 * x(2) &
        + 1.d0 / 60 * x(3) &
    ) / delta
end function

pure function lower(str)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lower

    integer :: i

    do i = 1, len(str)
        select case (str(i:i))
            case ('A':'Z')
                lower(i:i) = achar(iachar(str(i:i)) + 32)
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
        i = idx(i_dim) + 1
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

    allocate (this%timestamps(n), source=0_i8)
    allocate (this%counts(n), source=0_i8)
    allocate (this%levels(n), source=0)
end subroutine

subroutine clock_clock(this, id)
    class(clock_t), intent(inout) :: this
    integer, intent(in) :: id

    integer(i8) :: cnt, rate, cnt_max
    integer :: id_

    if (.not. this%active) return
    call system_clock(cnt, rate, cnt_max)
    id_ = abs(id)
    if (id > 0) then
        this%timestamps(id_) = this%timestamps(id_) - cnt
        this%levels(id_) = this%level
        this%level = this%level + 1
    else
        this%timestamps(id_) = this%timestamps(id_) + cnt
        this%counts(id_) = this%counts(id_) + 1
        this%level = this%level - 1
    end if
end subroutine

subroutine clock_print(this)
    class(clock_t), intent(inout) :: this

    integer(i8) :: cnt, rate, cnt_max
    integer :: i
    character(len=20) :: label

#ifdef WITH_MPI
    if (mpi_get_rank() /= 0) return
#endif
    call system_clock(cnt, rate, cnt_max)

    print '(A5,A7,A20,A10,A16)', 'id', 'level', 'label', 'count', 'time (s)'
    do i = 1, size(this%counts)
        if (this%counts(i) == 0) cycle
        select case (i)
        case (11); label = 'dipmat real'
        case (12); label = 'dipmat rec'
        case (13); label = 'P_EVR'
        case (14); label = 'mmul'
        case (15); label = 'force contractions'
        case (16); label = 'PDGETRF'
        case (17); label = 'PDGETRI'
        case (18); label = 'ELSI ev'
        case (20); label = 'MBD dipole'
        case (21); label = 'MBD eig'
        case (22); label = 'MBD forces'
        case (23); label = 'RPA eig'
        case (30); label = 'SCS dipole'
        case (32); label = 'SCS inv'
        case (33); label = 'SCS forces'
        case (50); label = 'SCS'
        case (51); label = 'MBD k-point'
        case (52); label = 'MBD'
        case (90); label = 'MBD@rsSCS'
        case (91); label = 'MBD@rsSCS forces'
        case default
            label = '('//trim(tostr(i))//')'
        end select
        print '(I5,I7,"  ",A20,I10,F16.6)', i, this%levels(i), label, this%counts(i), &
            dble(this%timestamps(i)) / rate
    end do
end subroutine

subroutine printer(str)
    character(len=*), intent(in) :: str

#ifdef WITH_MPI
    if (mpi_get_rank() /= 0) return
#endif
    print *, str
end subroutine

subroutine logger_debug(this, str)
    class(logger_t), intent(in) :: this
    character(len=*), intent(in) :: str

    if (this%level <= MBD_LOG_LVL_DEBUG) call this%printer(str)
end subroutine

subroutine logger_info(this, str)
    class(logger_t), intent(in) :: this
    character(len=*), intent(in) :: str

    if (this%level <= MBD_LOG_LVL_INFO) call this%printer(str)
end subroutine

subroutine logger_warn(this, str)
    class(logger_t), intent(in) :: this
    character(len=*), intent(in) :: str

    if (this%level <= MBD_LOG_LVL_WARN) call this%printer(str)
end subroutine

subroutine logger_error(this, str)
    class(logger_t), intent(in) :: this
    character(len=*), intent(in) :: str

    if (this%level <= MBD_LOG_LVL_ERROR) call this%printer(str)
end subroutine

logical function is_true(val) result(res)
    logical, intent(in), optional :: val

    res = .false.
    if (present(val)) res = val
end function

end module
