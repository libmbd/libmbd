#!/bin/bash
set -ev

if [[ $TRAVIS_OS_NAME = osx ]]; then
    rm -rf /usr/local/include/c++
    brew update
    brew install gcc open-mpi scalapack
    brew upgrade python3
fi
pip3 install tox tox-venv
if [[ $TOXENV != doc ]]; then
    CMAKE_FLAGS=()
    if [[ $TOXENV = *codecov* ]]; then
        CMAKE_FLAGS+=(-DCMAKE_Fortran_FLAGS="-fprofile-arcs -ftest-coverage")
    fi
    if [[ $TOXENV = *mpi* ]]; then
        CMAKE_FLAGS+=(-DENABLE_SCALAPACK_MPI=ON)
        if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
            CMAKE_FLAGS+=(-DSCALAPACK_LIBRARIES="-lscalapack-openmpi -lblacs-openmpi")
        fi
    fi
    mkdir build
    cd build
    cmake .. "${CMAKE_FLAGS[@]}"
fi
