#!/bin/bash
set -ev
$PYTHON -m pip install --user cffi numpy pytest mpi4py
if [[ "$WITH_PIP" ]]; then
    $PYTHON -m pip install --user .
else
    CMAKE_FLAGS=()
    if [[ "$CODECOV" ]]; then
        CMAKE_FLAGS+=(-DCMAKE_Fortran_FLAGS="-fprofile-arcs -ftest-coverage")
    fi
    if [[ "$MPI_NODES" ]]; then
        CMAKE_FLAGS+=(-DENABLE_SCALAPACK_MPI=ON)
        if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
            CMAKE_FLAGS+=(-DSCALAPACK_LIBRARIES="-lscalapack-openmpi -lblacs-openmpi")
        fi
    fi
    mkdir build
    pushd build
    cmake .. "${CMAKE_FLAGS[@]}"
    make
    popd
    $PYTHON setup.py build_ext -i
    $PYTHON -m pip install --user -e .
fi
if [[ "$CODECOV" ]]; then
    $PYTHON -m pip install --user pytest-cov
fi
