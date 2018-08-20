# Libmbd — many-body dispersion library


[![Travis](https://img.shields.io/travis/azag0/libmbd.svg)](https://travis-ci.org/azag0/libmbd)
[![Codecov](https://img.shields.io/codecov/c/github/azag0/libmbd.svg)](https://codecov.io/gh/azag0/libmbd)
[![Python 2.7](https://img.shields.io/badge/Python-2.7-blue.svg)]()
[![Python 3.7](https://img.shields.io/badge/Python-3.7-blue.svg)]()
[![MPL 2.0 license](https://img.shields.io/github/license/azag0/libmbd.svg)](https://github.com/azag0/libmbd/blob/master/LICENSE)

This project contains implementations of the [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) (MBD) method in several programming languages and frameworks:

- The Fortran implementation is the reference, most advanced implementation, with support for analytical gradients and distributed parallelism, and additional functionality beyond the MBD method itself. It provides a low-level and a high-level Fortran API, and a C API. Furthermore, Python bindings to the C API are provided.
- The Python/Numpy implementation is intended for prototyping, and as a high-level language reference.
- The Python/Tensorflow implemntation is an experiment that should enable rapid prototyping of machine learning applications with MBD.

Since all implementations provide Python bindings, the project is structured as a Python package called Pymbd, however, the implementations in languages other than Python can be used as standalone libraries in their respective languages.

## Installing Pymbd

There are two basic ways how to install Pymbd. The Pip installation is easier and intended for basic usage. The Cmake installation is best for development, and supports distributed parallelism with Scalapack.

Both approaches require a Fortran compiler, Lapack, cFFI, and Numpy.  The parallel version also requires MPI and Scalapack. All these need to be installed before installing Pymbd.

On Ubuntu:

```bash
apt-get install gfortran libblas-dev liblapack-dev [mpi-default-dev mpi-default-bin libscalapack-mpi-dev]
pip install cffi numpy
```

On macOS:

```bash
brew install gcc [open-mpi scalapack]
pip install cffi numpy
```

### Using Pip

```
pip install git+https://github.com/azag0/libmbd.git
```

If you have Pytest installed, you can run tests with

```
pytest --pyargs pymbd -v --durations=3
```

### Using Cmake

This is the recommended way for developing or if installing via Pip runs into problems.

```
git clone https://github.com/azag0/libmbd.git && cd libmbd
mkdir build && pushd build && cmake .. [-DSCALAPACK=<Scalapack linking>] && make && popd
python setup.py build_ext -i
pip install -e .[MPI,tensorflow]
```

If one wants to build the serial version even though MPI is present on the system, the Cmake flag `-DMPI_Fortran_COMPILER=xxxxx` disables the Scalapack compilation.

Tests can be run with

```
make -C build check
pytest -v --durations=3
```

## Example

```python
from pymbd import mbd_energy_species, ang
from pymbd.fortran import MBDCalc

ene_py = mbd_energy_species(  # pure Python implementation
    [(0, 0, 0), (0, 0, 4*ang)], ['Ar', 'Ar'], [1, 1], 0.83
)
with MBDCalc() as calc:
    ene_f = calc.mbd_energy_species(  # Fortran implementation
        [(0, 0, 0), (0, 0, 4*ang)], ['Ar', 'Ar'], [1, 1], 0.83
    )
assert abs(ene_f-ene_py) < 1e-15
```
