! Any copyright is dedicated to the Public Domain.
! http://creativecommons.org/publicdomain/zero/1.0/

module mpi
    integer, parameter :: MPI_COMM_WORLD = 0
    integer, parameter :: MPI_DOUBLE_PRECISION = 0
    integer, parameter :: MPI_COMPLEX16 = 1
    integer, parameter :: MPI_SUM = 0
    integer, parameter :: MPI_SUCCESS = 0
end module

subroutine MPI_COMM_RANK(comm, rank, err)
    integer :: comm, rank, err

    rank = 0
end subroutine

subroutine MPI_ALLREDUCE(sendbuf, recvbuf, cnt, datatype, op, comm, err)
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_COMPLEX16

    external :: MPI_ALLREDUCE_COMPLEX16

    integer :: cnt, datatype, op, comm, err
    real(8) :: sendbuf(cnt), recvbuf(cnt)

    select case (datatype)
    case (MPI_DOUBLE_PRECISION)
        recvbuf(:) = sendbuf(:)
    case (MPI_COMPLEX16)
        call MPI_ALLREDUCE_COMPLEX16(sendbuf, recvbuf, cnt, datatype, op, comm, err)
    end select
    err = MPI_SUCCESS
end subroutine

subroutine MPI_ALLREDUCE_COMPLEX16(sendbuf, recvbuf, cnt, datatype, op, comm, err)
    integer :: cnt, datatype, op, comm, err
    complex(8) :: sendbuf(cnt), recvbuf(cnt)

    recvbuf(:) = sendbuf(:)
end subroutine

subroutine MPI_BCAST(buffer, cnt, datatype, root, comm, err)
    use mpi, only: MPI_SUCCESS

    integer :: cnt, datatype, root, comm, err
    real(8) :: buffer(cnt)

    err = MPI_SUCCESS
end subroutine
