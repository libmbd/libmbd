# `libmbd` — Many-body dispersion method


[![Travis](https://img.shields.io/travis/azag0/libmbd.svg)](https://travis-ci.org/azag0/libmbd)
[![Codecov](https://img.shields.io/codecov/c/github/azag0/libmbd.svg)](https://codecov.io/gh/azag0/libmbd)
[![Python 2.7](https://img.shields.io/badge/Python-2.7-blue.svg)]()
[![Python 3.6](https://img.shields.io/badge/Python-3.6-blue.svg)]()
[![MPL 2.0 license](https://img.shields.io/github/license/azag0/libmbd.svg)](https://github.com/azag0/libmbd/blob/master/LICENSE)

Python 2/3 package for calculating [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) energies.

Most functionality is implemented in the Fortran module in `src/mbd.F90`.

## Installing

There are two basic ways how to install pymbd. The pip-based installation is easier and intended for basic usage. The cmake-based installation is best for development.

Both ways require a Fortran compiler, Lapack, cFFI, and Numpy. All these need to be installed before installing pymbd.

On Ubuntu:

```bash
apt-get install gfortran libblas-dev liblapack-dev
pip install cffi numpy
```

On macOS:

```bash
brew install gcc
pip install cffi numpy
```

### Using pip

```
pip install git+https://github.com/azag0/libmbd.git
```

If you have pytest installed, you can run tests with

```
pytest --pyargs pymbd -v --durations=3
```

### Using Cmake

This is the recommended way for developing or if installing via pip runs into problems.

```
git clone https://github.com/azag0/libmbd.git && cd libmbd
mkdir build && pushd build && cmake .. && make && popd
python setup.py build_ext -i
pip install -e .
```

Tests can be run with

```
make -C build check
pytest -v --durations=3
```

## Usage

Pymbd doesn't have any input files, it is called directly from Python scripts. 

For examples, see the [tests](https://github.com/azag0/libmbd/blob/master/pymbd/test_pymbd.py).
