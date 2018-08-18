#!/bin/bash
set -ev
if [[ "$WITH_PIP" ]]; then
    cd
    $PYTHON -m pytest --pyargs pymbd -vs --durations=3
else
    if [[ "$MPI_NODES" ]]; then
        PREFIX="mpirun -n $MPI_NODES"
    fi
    make -C build mbd_tests mbd_api_tests
    $PREFIX build/mbd_tests
    $PREFIX build/mbd_api_tests
    if [[ "$CODECOV" ]]; then
        EXTRA_FLAGS="--cov=./"
    fi
    $PREFIX $PYTHON -m pytest $EXTRA_FLAGS -vs --durations=3
fi
