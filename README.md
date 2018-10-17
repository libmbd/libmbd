# Libmbd — many-body dispersion library

[![build](https://img.shields.io/travis/azag0/libmbd/master.svg)](https://travis-ci.org/azag0/libmbd)
[![coverage](https://img.shields.io/codecov/c/github/azag0/libmbd.svg)](https://codecov.io/gh/azag0/libmbd)
![python](https://img.shields.io/pypi/pyversions/pymbd.svg)
[![release](https://img.shields.io/github/release/azag0/libmbd.svg)](https://github.com/azag0/libmbd/releases)
[![conda](https://img.shields.io/conda/v/azag0/pymbd.svg)](https://anaconda.org/azag0/pymbd)
[![pypi](https://img.shields.io/pypi/v/pymbd.svg)](https://pypi.org/project/pymbd/)
[![commits since](https://img.shields.io/github/commits-since/azag0/libmbd/latest.svg)](https://github.com/azag0/libmbd/releases)
[![last commit](https://img.shields.io/github/last-commit/azag0/libmbd.svg)](https://github.com/azag0/libmbd/commits/master)
[![license](https://img.shields.io/github/license/azag0/libmbd.svg)](https://github.com/azag0/libmbd/blob/master/LICENSE)

This project contains implementations of the [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) (MBD) method in several programming languages and frameworks:

- The Fortran implementation is the reference, most advanced implementation, with support for analytical gradients and distributed parallelism, and additional functionality beyond the MBD method itself. It provides a low-level and a high-level Fortran API, and a C API. Furthermore, Python bindings to the C API are provided.
- The Python/Numpy implementation is intended for prototyping, and as a high-level language reference.
- The Python/Tensorflow implementation is an experiment that should enable rapid prototyping of machine learning applications with MBD.

The Python-based implementations as well as Python bindings to the Libmbd C API are accessible from the Python package called Pymbd.

## Installing Pymbd

The easiest way to get Pymbd is to install the Pymbd [Conda](https://conda.io/docs/) package, which ships with pre-built Libmbd.

```
conda install -c azag0 pymbd
```

Alternatively, if you have Libmbd installed on your system (see below), you can install Pymbd via Pip, in which case it links against the installed Libmbd. To support Libmbd built with ScaLAPACK/MPI, the `MPI` extras is required.

```
pip install pymbd  # or pymbd[MPI]
```

In both cases, tests can be run with Pytest.

```
pytest -v --durations=3 --pyargs pymbd
```

If you don’t need the Fortran bindings in Pymbd, you can install it without the C extension, in which case `pymbd.fortran` becomes unimportable

```
pip install pymbd --install-option="--no-ext"
```

## Installing Libmbd

Libmbd uses CMake for building and installation, and requires a Fortran compiler, LAPACK, and optionally ScaLAPACK/MPI.

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

For development, Libmbd doesn't have to be installed on the system, and Pymbd can be linked against Libmbd in `./build`.

```
git clone https://github.com/azag0/libmbd.git && cd libmbd
mkdir build
(cd build && cmake .. -DENABLE_SCALAPACK_MPI=ON)
make -C build
make -C build check
python3 -m venv venv && source venv/bin/activate
pip install cffi numpy scipy mpi4py
python setup.py build_ext -i -Isrc -Lbuild/src -Rbuild/src
pytest -v --durations=3
```
