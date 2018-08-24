from .fortran import with_mpi

if with_mpi:
    from mpi4py import MPI
    from functools import wraps
    import py._io.terminalwriter

    rank = MPI.COMM_WORLD.Get_rank()
    _write_out = py._io.terminalwriter.write_out

    @wraps(_write_out)
    def write_out_wrapper(*args, **kwargs):
        if rank == 0:
            _write_out(*args, **kwargs)

    py._io.terminalwriter.write_out = write_out_wrapper
