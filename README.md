# Libmbd

![checks](https://img.shields.io/github/checks-status/libmbd/libmbd/master.svg)
[![coverage](https://img.shields.io/codecov/c/github/libmbd/libmbd.svg)](https://codecov.io/gh/libmbd/libmbd)
![python](https://img.shields.io/pypi/pyversions/pymbd.svg)
[![conda](https://img.shields.io/conda/vn/conda-forge/libmbd.svg)](https://anaconda.org/conda-forge/libmbd)
[![pypi](https://img.shields.io/pypi/v/pymbd.svg)](https://pypi.org/project/pymbd/)
[![commits since](https://img.shields.io/github/commits-since/libmbd/libmbd/latest.svg)](https://github.com/libmbd/libmbd/releases)
[![last commit](https://img.shields.io/github/last-commit/libmbd/libmbd.svg)](https://github.com/libmbd/libmbd/commits/master)
[![license](https://img.shields.io/github/license/libmbd/libmbd.svg)](https://github.com/libmbd/libmbd/blob/master/LICENSE)
[![code style](https://img.shields.io/badge/code%20style-black-202020.svg)](https://github.com/ambv/black)
[![chat](https://img.shields.io/gitter/room/libmbd/community)](https://gitter.im/libmbd/community)
[![doi](https://img.shields.io/badge/doi-10%2Fg67f-blue)](http://doi.org/g67f)

Libmbd implements the [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) (MBD) method in several programming languages and frameworks:

- The Fortran implementation is the reference, most advanced implementation, with support for analytical gradients and distributed parallelism, and additional functionality beyond the MBD method itself. It provides a low-level and a high-level Fortran API, as well as a C API. Furthermore, Python bindings to the C API are provided.
- The Python/Numpy implementation is intended for prototyping, and as a high-level language reference.
- The Python/Tensorflow implementation is an experiment that should enable rapid prototyping of machine learning applications with MBD.

The Python-based implementations as well as Python bindings to the Libmbd C API are accessible from the Python package called Pymbd.

Libmbd is included in [FHI-aims](https://aimsclub.fhi-berlin.mpg.de), [Quantum Espresso](https://www.quantum-espresso.org), [DFTB+](https://dftbplus.org), and [ESL Bundle](https://esl.cecam.org/bundle/).

## Installing

**TL;DR** Install prebuilt Libmbd binaries via [Conda-forge](https://conda-forge.org) and Pymbd with [Pip](https://pip.pypa.io/en/stable/quickstart/).

```
conda install -c conda-forge libmbd
pip install pymbd
```

One can also install the ScaLAPACK/MPI version.

```
conda install -c conda-forge 'libmbd=*=mpi_*' mpi4py
pip install pymbd[mpi]
```

Verify installation with

```
$ python -m pymbd
Expected energy:   -0.0002462647623815428
Calculated energy: -0.0002462647623817456
```

###  Libmbd

Libmbd uses CMake for compiling and installing, and requires a Fortran compiler, LAPACK, and optionally ScaLAPACK/MPI.

On Ubuntu:

```bash
apt-get install gfortran libblas-dev liblapack-dev [mpi-default-dev mpi-default-bin libscalapack-mpi-dev]
```

On macOS:

```bash
brew install gcc [open-mpi scalapack]
```

The compiling and installation can then proceed with

```
cmake -B build [-DENABLE_SCALAPACK_MPI=ON]
make -C build install
[ctest --test-dir build]
```

This installs the Libmbd shared library, C API header file,  high-level Fortran API module file, and Cmake package files, and optionally runs tests.

### Pymbd

Pymbd can be installed and updated using [Pip](https://pip.pypa.io/en/stable/quickstart/), but requires installed Libmbd as a dependency (see above).

```
pip install pymbd
```

To support Libmbd built with ScaLAPACK/MPI, the `mpi` extras is required, which installs `mpi4py` as an extra dependency. In this case one has to make sure that `mpi4py` is linked against the same MPI library as Libmbd (for instance by compiling both manually, or installing both via Conda-forge).

```
pip install pymbd[mpi]
```

If Libmbd is installed in a non-standard location, you can point Pymbd to it with

```
env LIBMBD_PREFIX=<path to Libmbd install prefix> pip install pymbd
```

If you don’t need the Fortran bindings in Pymbd, you can install it without the C extension, in which case `pymbd.fortran` becomes unimportable:

```
env LIBMBD_PREFIX= pip install pymbd
```


## Examples

```python
from pymbd import mbd_energy_species
from pymbd.fortran import MBDGeom

# pure Python implementation
energy = mbd_energy_species([(0, 0, 0), (0, 0, 7.5)], ['Ar', 'Ar'], [1, 1], 0.83)
# Fortran implementation
energy = MBDGeom([(0, 0, 0), (0, 0, 7.5)]).mbd_energy_species(
    ['Ar', 'Ar'], [1, 1], 0.83
)
```

```fortran
use mbd, only: mbd_input_t, mbd_calc_t

type(mbd_input_t) :: inp
type(mbd_calc_t) :: calc
real(8) :: energy, gradients(3, 2)
integer :: code
character(200) :: origin, msg

inp%atom_types = ['Ar', 'Ar']
inp%coords = reshape([0d0, 0d0, 0d0, 0d0, 0d0, 7.5d0], [3, 2])
inp%xc = 'pbe'
call calc%init(inp)
call calc%get_exception(code, origin, msg)
if (code > 0) then
    print *, msg
    stop 1
end if
call calc%update_vdw_params_from_ratios([0.98d0, 0.98d0])
call calc%evaluate_vdw_method(energy)
call calc%get_gradients(gradients)
call calc%destroy()
```

## Links

- Libmbd documentation: https://libmbd.github.io
- Pymbd documentation: https://libmbd.github.io/pymbd

## Developing

For development, a top-level `Makefile` is included, which configures and compiles Libmbd, compiles the Pymbd C extension, and runs both Libmbd and Pymbd tests.

```
git clone https://github.com/libmbd/libmbd.git && cd libmbd
python3 -m venv venv && source venv/bin/activate
make
# development work...
make
```
