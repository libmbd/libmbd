! Any copyright is dedicated to the Public Domain.
! http://creativecommons.org/publicdomain/zero/1.0/

module mpi
    integer, parameter :: MPI_COMM_WORLD = 0
    integer, parameter :: MPI_DOUBLE_PRECISION = 0
    integer, parameter :: MPI_COMPLEX16 = 1
    integer, parameter :: MPI_SUM = 0
    integer, parameter :: MPI_SUCCESS = 0
end module

subroutine MPI_Comm_rank(comm, rank, err)
    integer :: comm, rank, err

    rank = 0
end subroutine mpi_comm_rank

subroutine MPI_Allreduce(sendbuf, recvbuf, cnt, datatype, op, comm, err)
    use mpi

    integer :: cnt, datatype, op, comm, err
    real(8) :: sendbuf(cnt), recvbuf(cnt)

    select case (datatype)
    case (MPI_DOUBLE_PRECISION)
        recvbuf(:) = sendbuf(:)
    case (MPI_COMPLEX16)
        call MPI_Allreduce_complex16(sendbuf, recvbuf, cnt, datatype, op, comm, err)
    end select
    err = MPI_SUCCESS
end subroutine

subroutine MPI_Allreduce_complex16(sendbuf, recvbuf, cnt, datatype, op, comm, err)
    integer :: cnt, datatype, op, comm, err
    complex(8) :: sendbuf(cnt), recvbuf(cnt)

    recvbuf(:) = sendbuf(:)
end subroutine

subroutine MPI_Bcast(buffer, cnt, datatype, root, comm, err)
    use mpi

    integer :: cnt, datatype, root, comm, err
    real(8) :: buffer(cnt)

    err = MPI_SUCCESS
end subroutine
