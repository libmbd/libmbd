# `pymbd` — Many-body dispersion method

[![](https://travis-ci.org/azag0/pymbd.svg?branch=master)](https://travis-ci.org/azag0/pymbd)

Python 2/3 package for calculating [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) energies.

Most functionality is implemented in the Fortran module in `src/mbd.f90`.

## Installation

There are two basic ways how to install pymbd.

Both ways require a Fortran compiler, Lapack, CFFI, and Numpy. Optional MPI support requires mpi4py installed from source. All these need to be installed before installing pymbd.

On Ubuntu:

```bash
apt-get install gfortran libblas-dev liblapack-dev [mpi-default-dev mpi-default-bin]
pip install cffi numpy [mpi4py]
```

On macOS:

```bash
brew install gcc [mpich]
```

### Using pip

```
pip install git+https://github.com/azag0/pymbd.git
```

If you have pytest installed, you can run tests with

```
pytest --pyargs pymbd -v --durations=3
```

### Using Cmake

This is the recommended way for developing or if installing via pip runs into problems.

```
git clone https://github.com/azag0/pymbd.git
mkdir build && pushd build && cmake .. && make && popd
pip install -r requirements.txt
python setup.py build_ext -i
```

Tests can be run with

```
make -C build check
pytest -v --durations=3
```

## Usage

Pymbd doesn't have any input files, it is called directly from Python scripts. 

For examples, see the [tests](https://github.com/azag0/pymbd/blob/master/pymbd/test_pymbd.py).
