module mbd_interface

use mpi

implicit none

interface sync_sum
    module procedure sync_sum_array_dble_
    module procedure sync_sum_array_cmplx_
    module procedure sync_sum_matrix_dble_
    module procedure sync_sum_matrix_cmplx_
    module procedure sync_sum_vector_dble_
    module procedure sync_sum_vector_cmplx_
    module procedure sync_sum_dble_
    module procedure sync_sum_cmplx_
end interface

interface broadcast
    module procedure broadcast_array_dble_
end interface

contains

subroutine sync_sum_array_dble_(array, n_array)
    implicit none

    integer, intent(in) :: n_array
    real*8, intent(inout) :: array(n_array)

    real*8 :: array_buff(n_array)
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        array, array_buff, n_array, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    array = array_buff
end subroutine

subroutine sync_sum_array_cmplx_(array, n_array)
    implicit none

    integer, intent(in) :: n_array
    complex(kind=8), intent(inout) :: array(n_array)

    complex(kind=8) :: array_buff(n_array)
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        array, array_buff, n_array, MPI_COMPLEX16, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    array = array_buff
end subroutine

subroutine sync_sum_dble_(x)
    implicit none

    real*8, intent(inout) :: x
    real*8 :: x_buff
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        x, x_buff, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    x = x_buff
end subroutine

subroutine sync_sum_cmplx_(x)
    implicit none

    complex(kind=8), intent(inout) :: x
    complex(kind=8) :: x_buff
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        x, x_buff, 1, MPI_COMPLEX16, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    x = x_buff
end subroutine

subroutine sync_sum_vector_dble_(x)
    implicit none

    real*8, intent(inout) :: x(:)

    call syn_sum(x, size(x))
end subroutine

subroutine sync_sum_vector_cmplx_(x)
    implicit none

    complex(kind=8), intent(inout) :: x(:)

    call syn_sum(x, size(x))
end subroutine

subroutine sync_sum_matrix_dble_(x)
    implicit none

    real*8, intent(inout) :: x(:, :)

    call syn_sum(x, size(x))
end subroutine

subroutine sync_sum_matrix_cmplx_(x)
    implicit none

    complex(kind=8), intent(inout) :: x(:, :)

    call syn_sum(x, size(x))
end subroutine

subroutine broadcast_array_dble_(array, n_array)
    implicit none

    integer, intent(in) :: n_array
    real*8, intent(inout) :: array(n_array)

    integer :: mpi_err

    call MPI_BCAST(array, n_array, MPI_DOUBLE_PRECISION, 0, &
        MPI_COMM_WORLD, mpi_err)
end subroutine

subroutine print_log(str)
    implicit none

    character(len=*), intent(in) :: str
    integer :: myid, error

    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (6, *) str
    end if
end subroutine

subroutine print_warning(str)
    implicit none

    character(len=*), intent(in) :: str
    integer :: myid, error

    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (0, *) "Warning: " // str
    end if
end subroutine

subroutine print_error(str)
    implicit none

    character(len=*), intent(in) :: str
    integer :: myid, error

    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (0, *) "Error: " // str
    end if
end subroutine

end module mbd_interface
