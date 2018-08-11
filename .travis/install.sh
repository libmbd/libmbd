#!/bin/bash
set -ev
$PYTHON -m pip install --user cffi numpy pytest
if [[ "$MPI_NODES" ]]; then
    $PYTHON -m pip install --user mpi4py
fi
if [[ "$WITH_PIP" ]]; then
    $PYTHON -m pip install --user .
else
    CMAKE_FLAGS=()
    if [[ "$CODECOV" ]]; then
        CMAKE_FLAGS+=(-DCMAKE_Fortran_FLAGS="-fprofile-arcs -ftest-coverage")
    fi
    if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
        SCALAPACKLIB=-lscalapack
    else
        SCALAPACKLIB="-lscalapack-openmpi -lblacs-openmpi"
    fi
    if [[ "$MPI_NODES" ]]; then
        CMAKE_FLAGS+=(-DSCALAPACK="$SCALAPACKLIB")
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
