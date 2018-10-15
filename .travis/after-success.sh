#!/bin/bash
set -ev

if [[ $TOXENV = *codecov* ]]; then
    cp build/src/CMakeFiles/mbd.dir/*.f90.gc?? src/
    bash <(curl -s https://codecov.io/bash) -F pytest
fi
