#!/bin/bash
set -ev
if [[ "$WITH_PIP" ]]; then
    cd
    $PYTHON -m pytest --pyargs pymbd -vs --durations=3
else
    if [[ "$MPI_NODES" ]]; then
        PREFIX="mpirun -n $MPI_NODES"
    fi
    $PREFIX make -C build check
    if [[ "$CODECOV" ]]; then
        EXTRA_FLAGS="--cov=./"
    fi
    $PREFIX $PYTHON -m pytest $EXTRA_FLAGS -vs --durations=3
fi
