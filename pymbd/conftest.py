try:
    from mpi4py import MPI
except ImportError:
    pass
else:  # wrap py._io.terminalwriter.write_out to print only on rank 0
    rank = MPI.COMM_WORLD.Get_rank()

    from functools import wraps
    import py._io.terminalwriter

    _write_out = py._io.terminalwriter.write_out
    @wraps(_write_out)
    def write_out_wrapper(*args, **kwargs):
        if rank == 0:
            _write_out(*args, **kwargs)

    py._io.terminalwriter.write_out = write_out_wrapper
