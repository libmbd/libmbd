! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_mpi

use mbd_constants
#ifdef WITH_MPIFH
include 'mpif.h'
#else
use mpi
#endif

implicit none

private
public :: mpi_all_reduce, MPI_COMM_WORLD, MPI_COMM_RANK, MPI_COMM_SIZE

interface mpi_all_reduce
    module procedure mpi_all_reduce_real_0d
    module procedure mpi_all_reduce_real_1d
    module procedure mpi_all_reduce_real_2d
end interface

contains

subroutine mpi_all_reduce_real_0d(x, comm)
    real(dp), intent(inout) :: x
    integer, intent(in) :: comm

    real(dp) :: x_buffer
    integer :: ierr

    call MPI_ALLREDUCE(x, x_buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    x = x_buffer
end subroutine

subroutine mpi_all_reduce_real(x, n, comm)
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n)
    integer, intent(in) :: comm

    real(dp), allocatable :: x_buffer(:)
    integer :: ierr

    allocate (x_buffer(n))
    call MPI_ALLREDUCE(x, x_buffer, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    x = x_buffer
end subroutine

subroutine mpi_all_reduce_real_1d(x, comm)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: comm

    call mpi_all_reduce_real(x, size(x), comm)
end subroutine

subroutine mpi_all_reduce_real_2d(x, comm)
    real(dp), intent(inout) :: x(:, :)
    integer, intent(in) :: comm

    call mpi_all_reduce_real(x, size(x), comm)
end subroutine

end module
