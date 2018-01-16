! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_mpi

use mpi, only: MPI_COMM_WORLD, MPI_COMPLEX16, MPI_DOUBLE_PRECISION, MPI_SUM
use mbd_common, only: dp

implicit none

private

public :: sync_sum, broadcast, MPI_COMM_WORLD

interface sync_sum
    module procedure sync_sum_dble_
    module procedure sync_sum_vector_dble_
    module procedure sync_sum_matrix_dble_
    module procedure sync_sum_3d_dble_
    module procedure sync_sum_4d_dble_
    module procedure sync_sum_cmplx_
    module procedure sync_sum_vector_cmplx_
    module procedure sync_sum_matrix_cmplx_
    module procedure sync_sum_3d_cmplx_
    module procedure sync_sum_4d_cmplx_
end interface

interface broadcast
    module procedure broadcast_vector_dble_
    module procedure broadcast_matrix_dble_
    module procedure broadcast_3d_dble_
    module procedure broadcast_4d_dble_
    module procedure broadcast_vector_cmplx_
    module procedure broadcast_matrix_cmplx_
end interface

external :: MPI_BCAST, MPI_ALLREDUCE

contains


subroutine sync_sum_dble_(x, comm)
    real(dp), intent(inout) :: x
    integer, intent(in) :: comm

    real(dp) :: x_buff
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        x, x_buff, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, comm, mpi_err)
    x = x_buff
end subroutine


subroutine sync_sum_array_dble_(array, n_array, comm)
    integer, intent(in) :: n_array
    real(dp), intent(inout) :: array(n_array)
    integer, intent(in) :: comm

    real(dp) :: array_buff(n_array)
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        array, array_buff, n_array, MPI_DOUBLE_PRECISION, &
        MPI_SUM, comm, mpi_err)
    array = array_buff
end subroutine


subroutine sync_sum_vector_dble_(x, comm)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: comm

    call sync_sum_array_dble_(x, size(x), comm)
end subroutine


subroutine sync_sum_matrix_dble_(x, comm)
    real(dp), intent(inout) :: x(:, :)
    integer, intent(in) :: comm

    call  sync_sum_array_dble_(x, size(x), comm)
end subroutine


subroutine sync_sum_3d_dble_(x, comm)
    real(dp), intent(inout) :: x(:, :, :)
    integer, intent(in) :: comm

    call  sync_sum_array_dble_(x, size(x), comm)
end subroutine


subroutine sync_sum_4d_dble_(x, comm)
    real(dp), intent(inout) :: x(:, :, :, :)
    integer, intent(in) :: comm

    call  sync_sum_array_dble_(x, size(x), comm)
end subroutine


subroutine sync_sum_cmplx_(x, comm)
    complex(kind=8), intent(inout) :: x
    integer, intent(in) :: comm

    complex(kind=8) :: x_buff
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        x, x_buff, 1, MPI_COMPLEX16, &
        MPI_SUM, comm, mpi_err)
    x = x_buff
end subroutine


subroutine sync_sum_array_cmplx_(array, n_array, comm)
    integer, intent(in) :: n_array
    complex(kind=8), intent(inout) :: array(n_array)
    integer, intent(in) :: comm

    complex(kind=8) :: array_buff(n_array)
    integer :: mpi_err

    call MPI_ALLREDUCE( &
        array, array_buff, n_array, MPI_COMPLEX16, &
        MPI_SUM, comm, mpi_err)
    array = array_buff
end subroutine


subroutine sync_sum_vector_cmplx_(x, comm)
    complex(kind=8), intent(inout) :: x(:)
    integer, intent(in) :: comm

    call sync_sum_array_cmplx_(x, size(x), comm)
end subroutine


subroutine sync_sum_matrix_cmplx_(x, comm)
    complex(kind=8), intent(inout) :: x(:, :)
    integer, intent(in) :: comm

    call  sync_sum_array_cmplx_(x, size(x), comm)
end subroutine


subroutine sync_sum_3d_cmplx_(x, comm)
    complex(kind=8), intent(inout) :: x(:, :, :)
    integer, intent(in) :: comm

    call  sync_sum_array_cmplx_(x, size(x), comm)
end subroutine


subroutine sync_sum_4d_cmplx_(x, comm)
    complex(kind=8), intent(inout) :: x(:, :, :, :)
    integer, intent(in) :: comm

    call  sync_sum_array_cmplx_(x, size(x), comm)
end subroutine


subroutine broadcast_array_dble_(array, n_array, comm)
    integer, intent(in) :: n_array
    real(dp), intent(inout) :: array(n_array)
    integer, intent(in) :: comm

    integer :: mpi_err

    call MPI_BCAST(array, n_array, MPI_DOUBLE_PRECISION, 0, comm, mpi_err)
end subroutine


subroutine broadcast_vector_dble_(x, comm)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: comm

    call broadcast_array_dble_(x, size(x), comm)
end subroutine


subroutine broadcast_matrix_dble_(x, comm)
    real(dp), intent(inout) :: x(:, :)
    integer, intent(in) :: comm

    call broadcast_array_dble_(x, size(x), comm)
end subroutine


subroutine broadcast_3d_dble_(x, comm)
    real(dp), intent(inout) :: x(:, :, :)
    integer, intent(in) :: comm

    call broadcast_array_dble_(x, size(x), comm)
end subroutine


subroutine broadcast_4d_dble_(x, comm)
    real(dp), intent(inout) :: x(:, :, :, :)
    integer, intent(in) :: comm

    call broadcast_array_dble_(x, size(x), comm)
end subroutine


subroutine broadcast_array_cmplx_(array, n_array, comm)
    integer, intent(in) :: n_array
    complex(dp), intent(inout) :: array(n_array)
    integer, intent(in) :: comm

    integer :: mpi_err

    call MPI_BCAST(array, n_array, MPI_COMPLEX16, 0, comm, mpi_err)
end subroutine


subroutine broadcast_vector_cmplx_(x, comm)
    complex(dp), intent(inout) :: x(:)
    integer, intent(in) :: comm

    call broadcast_array_cmplx_(x, size(x), comm)
end subroutine


subroutine broadcast_matrix_cmplx_(x, comm)
    complex(dp), intent(inout) :: x(:, :)
    integer, intent(in) :: comm

    call broadcast_array_cmplx_(x, size(x), comm)
end subroutine

end module mbd_mpi
