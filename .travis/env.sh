CMAKE_EXTRA_FLAGS='-DCMAKE_Fortran_FLAGS="-fprofile-arcs -ftest-coverage"'
if [[ $CMAKE_FLAGS = *"ENABLE_SCALAPACK_MPI=ON"* ]]; then
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        CMAKE_EXTRA_FLAGS="$CMAKE_EXTRA_FLAGS"' -DSCALAPACK_LIBRARIES="-lscalapack-openmpi -lblacs-openmpi"'
    fi
fi
export CMAKE_FLAGS="${CMAKE_FLAGS} ${CMAKE_EXTRA_FLAGS}"
echo CMAKE_FLAGS=$CMAKE_FLAGS
export PYTEST_FLAGS="--cov=pymbd --cov-report=xml"
echo PYTEST_FLAGS=$PYTEST_FLAGS
