module mbd_interface

use mpi

implicit none

contains

subroutine sync_sum_array(array, n_array)
    implicit none

    integer, intent(in) :: n_array
    real*8, intent(inout) :: array(n_array)

    real*8 :: array_buff(n_array)
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        array, array_buff, n_array, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    array = array_buff
end subroutine sync_sum_array

subroutine broadcast_array(array, n_array)
    implicit none

    integer, intent(in) :: n_array
    real*8, intent(inout) :: array(n_array)

    integer :: mpi_err

    call MPI_BCAST(array, n_array, MPI_DOUBLE_PRECISION, 0, &
        MPI_COMM_WORLD, mpi_err)
end subroutine broadcast_array

subroutine sync_sum_array_c (array, n_array)
    implicit none

    integer, intent(in) :: n_array
    complex(kind=8), intent(inout) :: array(n_array)

    complex(kind=8) :: array_buff(n_array)
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        array, array_buff, n_array, MPI_COMPLEX16, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    array = array_buff
end subroutine sync_sum_array_c

subroutine sync_sum_number(num)
    implicit none

    real*8, intent(inout) :: num
    real*8 :: num_buff
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        num, num_buff, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    num = num_buff
end subroutine sync_sum_number

subroutine print_log(str)
    implicit none

    character(len=*), intent(in) :: str
    integer :: myid, error

    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (6, *) str
    end if
end subroutine print_log

subroutine print_warning(str)
    implicit none

    character(len=*), intent(in) :: str
    integer :: myid, error

    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (0, *) "Warning: " // str
    end if
end subroutine print_warning

subroutine print_error(str)
    implicit none

    character(len=*), intent(in) :: str
    integer :: myid, error

    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (0, *) "Error: " // str
    end if
end subroutine print_error

end module mbd_interface
