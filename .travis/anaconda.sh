#!/bin/bash
set -ev

case $TRAVIS_OS_NAME in
    osx)
        ANACONDA_OS_TAG=MacOSX
        ;;
    linux)
        ANACONDA_OS_TAG=Linux
        ;;
esac
wget https://repo.continuum.io/miniconda/Miniconda3-latest-$ANACONDA_OS_TAG-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda config --set always_yes yes
conda update -q conda
conda install conda-build anaconda-client
conda build .conda
anaconda -t $ANACONDA_TOKEN upload -u libmbd $(conda build --output .)
