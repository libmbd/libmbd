# `pymbd` — Many-body dispersion method

[![](https://travis-ci.org/azag0/pymbd.svg?branch=master)](https://travis-ci.org/azag0/pymbd)

Python 2/3 package for calculating [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) energies.

Most functionality is implemented in the Fortran module in `src/mbd.f90`.

## Installation

There are two basic ways how to install pymbd.

### Using pip

Installation of pymbd with pip requires cffi and Numpy with properly configured Lapack libraries.

On macOS, this requirement is almost always satisfied with whatever Python with simple

```
$PYTHON -m pip install cffi numpy
```

This can be followed with

```
$PYTHON -m pip install git+https://github.com/azag0/pymbd.git
```

If you have pytest installed, you can run tests with

```
$PYTHON -m pytest --pyargs pymbd -v
```

On Linux, the pip-installed Numpy sometimes doesn't report the correct Lapack [1]. To avoid this issue, use the Anaconda Python with the no-MKL Numpy (if you have Anaconda already installed, [make sure](https://www.continuum.io/blog/developer-blog/anaconda-25-release-now-mkl-optimizations) that you don't use the MKL Numpy):

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
[...]
$PATH_TO_ANACONDA/bin/conda install cffi nomkl numy
```

To install pymbd, run

```
$PATH_TO_ANACONDA/bin/pip install git+https://github.com/azag0/pymbd.git
```

And finally test with

```
$PATH_TO_ANACONDA/bin/python -m pytest --pyargs pymbd -v
```

[1]: What really matters is whether `numpy.distutils.system_info.get_info('lapack_opt')` returns the correct information. On macOS, it always returns the Accelerate framework. With Anaconda and no-MKL Numpy, it always returns the internal openBLAS libraries. But as long as it returns a valid Lapack library, you can use whatever Python installation.

If you have mpi4py installed, pymbd will use its MPI Fortran compiler to build an MPI-enabled pymbd.

### Using Cmake

This is the recommended way for developing or if installing via pip runs into problems. This installation is also tested on [Travis](https://travis-ci.org/azag0/pymbd) both on Ubuntu (Python 2) and macOS (Python 3).

The prerequisities can be installed on Ubuntu with

```
sudo apt-get install gfortran liblapack3 mpi-default-dev python-mpi4py
```

On macOS with

```
brew install python3 gcc mpich
```

Then run 

```
git clone https://github.com/azag0/pymbd.git
mkdir build && cd build && cmake .. && make && cd ..
$PYTHON -m pip install --user -r requirements.txt
$PYTHON setup.py build_ext -i
```

Tests can be run with

```
make -C build check
$PYTHON -m pytest -v --durations=3
```

## Usage

Pymbd doesn't have any input files, it is called directly from Python scripts. 

For examples, see the [tests](https://github.com/azag0/pymbd/blob/master/pymbd/test_pymbd.py)
