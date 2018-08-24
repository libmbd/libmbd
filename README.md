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

The Python-based implementations as well as Python bindings to the Libmbd C API are accessible from the Python package called Pymbd.

## Installing Libmbd

Libmbd uses CMake for building and installation, and requires a Fortran compiler, Lapack, and optionally ScaLAPACK/MPI.

On Ubuntu:

```bash
apt-get install gfortran libblas-dev liblapack-dev [mpi-default-dev mpi-default-bin libscalapack-mpi-dev]
```

On macOS:

```bash
brew install gcc [open-mpi scalapack]
```

The building and installation can then proceed with

```
git clone https://github.com/azag0/libmbd.git && cd libmbd
mkdir build && cd build
cmake .. [-DENABLE_SCALAPACK_MPI=ON]
make
make install
```

This installs the Libmbd shared library, C API header file, and high-level Fortran API module file.

Tests can be run with

```
make check
```

## Installing Pymbd

Pymbd links against Libmbd, which can be either installed on the system, or built on the fly by Pip/Setuptools. The linking requires the cFFI Python package installed prior to installing Libmbd. If the installed Libmbd is built with ScaLAPACK/MPI, Mpi4py package is required. For the Pip/Setuptools build, Fortran compiler must be available on the system (ScaLAPACK/MPI is not supported by the Setuptools build), and Numpy must be installed prior to installing Pymbd.

```
pip install cffi [numpy] [mpi4py]
pip install git+https://github.com/azag0/libmbd.git
```

If you have Pytest installed, you can run tests with

```
pytest --pyargs pymbd -v --durations=3
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

## Developing

For development, Libmbd doesn't have to be installed on the system, and Pymbd can be linked against Libmbd in `./build`

```
git clone https://github.com/azag0/libmbd.git && cd libmbd
mkdir build && cd build
cmake .. -DENABLE_SCALAPACK_MPI=ON
make
make check
cd ..
pip install cffi mpi4py
python setup.py build_ext -i
pytest -v --durations=3
```
