# Libmbd — many-body dispersion library

[![build](https://img.shields.io/travis/azag0/libmbd/master.svg)](https://travis-ci.org/azag0/libmbd)
[![coverage](https://img.shields.io/codecov/c/github/azag0/libmbd.svg)](https://codecov.io/gh/azag0/libmbd)
![python](https://img.shields.io/pypi/pyversions/pymbd.svg)
[![pypi](https://img.shields.io/pypi/v/pymbd.svg)](https://pypi.org/project/pymbd/)
[![commits since](https://img.shields.io/github/commits-since/azag0/libmbd/latest.svg)](https://github.com/azag0/libmbd/releases)
[![last commit](https://img.shields.io/github/last-commit/azag0/libmbd.svg)](https://github.com/azag0/libmbd/commits/master)
[![license](https://img.shields.io/github/license/azag0/libmbd.svg)](https://github.com/azag0/libmbd/blob/master/LICENSE)

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

Pymbd links against Libmbd, which can be either installed on the system, or built on the fly by Pip/Setuptools. If the installed Libmbd is built with ScaLAPACK/MPI, the `MPI` extras is required. For the Pip/Setuptools build, Fortran compiler must be available on the system (ScaLAPACK/MPI is not supported by the Setuptools build), and Numpy must be installed prior to installing Pymbd.

```
pip install numpy
pip install pymbd  # or pymbd[MPI]
```

If you have Pytest installed, you can run tests with

```
pytest --pyargs pymbd -v --durations=3
```

## Examples

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

```fortran
use mbd, only: mbd_input, mbd_calculation

type(mbd_input) :: inp
type(mbd_calculation) :: calc
real(dp) :: energy, gradients(3, 2)
integer :: code
character(200) :: origin, msg

inp%atom_types = ['Ar', 'Ar']
inp%coords = reshape([0d0, 0d0, 0d0, 0d0, 0d0, 7.5d0], [3, 2])
inp%xc = 'pbe'
call calc%init(inp)
call calc%get_exception(code, origin, msg)
if (code > 0) then
    print *, msg
    stop
end if
call calc%update_vdw_params_from_ratios([0.98d0, 0.98d0])
call calc%get_energy(energy)
call calc%get_gradients(gradients)
call calc%destroy()
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
python3 -m venv venv && source venv/bin/activate
pip install cffi numpy scipy mpi4py
python setup.py build_ext -i
pytest -v --durations=3
```
