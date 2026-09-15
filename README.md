<h1 align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="doc/_static/libmbd-lockup-dark.svg">
  <img alt="libMBD" src="doc/_static/libmbd-lockup.svg" width="260">
</picture>
</h1>

[![tests](https://img.shields.io/github/actions/workflow/status/libmbd/libmbd/tests.yaml?branch=master&label=tests)](https://github.com/libmbd/libmbd/actions/workflows/tests.yaml)
[![coverage](https://img.shields.io/codecov/c/github/libmbd/libmbd.svg)](https://codecov.io/gh/libmbd/libmbd)
![python](https://img.shields.io/pypi/pyversions/pymbd.svg)
[![conda](https://img.shields.io/conda/vn/conda-forge/libmbd.svg)](https://anaconda.org/conda-forge/libmbd)
[![pypi](https://img.shields.io/pypi/v/pymbd.svg)](https://pypi.org/project/pymbd/)
[![commits since](https://img.shields.io/github/commits-since/libmbd/libmbd/latest.svg)](https://github.com/libmbd/libmbd/releases)
[![last commit](https://img.shields.io/github/last-commit/libmbd/libmbd.svg)](https://github.com/libmbd/libmbd/commits/master)
[![license](https://img.shields.io/github/license/libmbd/libmbd.svg)](https://github.com/libmbd/libmbd/blob/master/LICENSE)
[![ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![doi](https://img.shields.io/badge/doi-10%2Fk4bm-blue)](http://doi.org/k4bm)

> libMBD: A general-purpose package for scalable quantum many-body dispersion calculations. [J. Hermann](https://github.com/jhrmnn), [M. Stöhr](https://github.com/martin-stoehr), S. Góger, [S. Chaudhuri](https://github.com/shaychaudhuri), [B. Aradi](https://github.com/aradi), [R. J. Maurer](https://github.com/reinimaurer1) & A. Tkatchenko. [*J. Chem. Phys.* **159**, 174802](http://doi.org/k4bm) (2023)

libMBD implements the [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) (MBD) method in several programming languages and frameworks:

- The Fortran implementation is the reference, most advanced implementation, with support for analytical gradients and distributed parallelism, and additional functionality beyond the MBD method itself. It provides a low-level and a high-level Fortran API, as well as a C API. Furthermore, Python bindings to the C API are provided.
- The Python/Numpy implementation is intended for prototyping, and as a high-level language reference.
- The Python/Tensorflow implementation is an experiment that should enable rapid prototyping of machine learning applications with MBD.

The Python-based implementations as well as Python bindings to the libMBD C API are accessible from the Python package called pyMBD.

libMBD is included in [FHI-aims](https://aimsclub.fhi-berlin.mpg.de), [Quantum Espresso](https://www.quantum-espresso.org), [DFTB+](https://dftbplus.org), and [ESL Bundle](https://esl.cecam.org/en/bundle/).

## Installing

**TL;DR** Install prebuilt libMBD binaries via [Conda-forge](https://conda-forge.org) and pyMBD with [Pip](https://pip.pypa.io/en/stable/quickstart/).

```
conda install -c conda-forge libmbd
pip install pymbd
```

One can also install the ScaLAPACK/MPI version.

```
conda install -c conda-forge 'libmbd=*=mpi_*' mpi4py
pip install pymbd[mpi]
```

If an installation with support for MBD-ML is desired, do

```
conda install -c conda-forge libmbd
pip install .[mbd-ml]
```

To retrieve the exact dependencies used for the MBD-ML publication, 
a `requirements.txt` file has been deposited in the `mbd_ml_paper/` directory.

Verify installation with

```
$ python -m pymbd
Expected energy:   -0.0002462647623815428
Calculated energy: -0.0002462647623817456
```

### libMBD

libMBD uses CMake (≥ 3.22) for compiling and installing, and requires a Fortran compiler, LAPACK, and optionally ScaLAPACK/MPI.

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

This installs the libMBD shared library, C API header file,  high-level Fortran API module file, and Cmake package files, and optionally runs tests.

### pyMBD

pyMBD requires Python ≥ 3.10 and can be installed and updated using [Pip](https://pip.pypa.io/en/stable/quickstart/), but requires installed libMBD as a dependency (see above).

```
pip install pymbd
```

To support libMBD built with ScaLAPACK/MPI, the `mpi` extras is required, which installs `mpi4py` as an extra dependency. In this case one has to make sure that `mpi4py` is linked against the same MPI library as libMBD (for instance by compiling both manually, or installing both via Conda-forge).

```
pip install pymbd[mpi]
```

If libMBD is installed in a non-standard location, you can point pyMBD to it with

```
env LIBMBD_PREFIX=<path to libMBD install prefix> pip install pymbd
```

If you don’t need the Fortran bindings in pyMBD, you can install it without the C extension, in which case `pymbd.fortran` becomes unimportable:

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
use iso_fortran_env, only: real64

type(mbd_input_t) :: inp
type(mbd_calc_t) :: calc
real(real64) :: energy, gradients(3, 2)
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
### MBD-ML

To make use of the MBD-ML model, you only need to provide
your molecule or material structure as a ASE Atoms object
and the damping parameter beta corresponding to the
desired exchange-correlation functional. For example, 
to obtain the MBD energy and force correction for a PBE-level
calculation (beta=0.81), do

```python
import pymbd.mbd_ml
from ase import Atoms

mol = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])

mbd_props = pymbd.mbd_ml.mbd_properties_from_structure(atoms, beta=0.81)

E_MBD = mbd_props['E']
F_MBD = mbd_props['F']
```

for a periodic calculation, e.g a molecular crystal, you need to additionally
provide a k-mesh. For example for a urea crystal, with the following structure (.extxyz file)

`urea.extxyz`
```
16
Lattice="5.565 0.0 0.0 0.0 5.565 0.0 0.0 0.0 4.684" Properties=species:S:1:pos:R:3 pbc="T T T"
C        0.00000000       2.78250000       1.52698400
C        2.78250000       0.00000000       3.15701600
O        0.00000000       2.78250000       2.78838520
O        2.78250000       0.00000000       1.89561480
N        0.81193350       3.59443350       0.82719440
N        1.97056650       0.81193350       3.85680560
N        4.75306650       1.97056650       0.82719440
N        3.59443350       4.75306650       3.85680560
H        1.43298750       4.21548750       1.32416680
H        1.34951250       1.43298750       3.35983320
H        4.13201250       1.34951250       1.32416680
H        4.21548750       4.13201250       3.35983320
H        0.80191650       3.58441650       4.50600800
H        1.98058350       0.80191650       0.17799200
H        4.76308350       1.98058350       4.50600800
H        3.58441650       4.76308350       0.17799200
```

the same function call will output the stress tensor in addition to the MBD energy
and atomic forces.

```python
import pymbd.mbd_ml
import ase.io

urea = ase.io.read('urea.extxyz')
k_grid = (4, 4, 4)

mbd_props = pymbd.mbd_ml.mbd_properties_from_structure(atoms, beta=0.81, k_grid=k_grid)

E_MBD = mbd_props['E']
F_MBD = mbd_props['F']
S_MBD = mbd_props['S'] #stress tensor for periodic systems
```

If only the MBD ratios (a0/C6) are desired, use instead the `ratios_from_mbdml()` function (no k-mesh needed) 

```python
import pymbd.mbd_ml
import ase.io

urea = ase.io.read('urea.extxyz')

mbd_ratios = pymbd.mbd_ml.ratios_from_mbdml(urea)

a0 = mbd_ratios['a0']
c6 = mbd_ratios['c6']
```

To use the MBD-ML framework in combination with an ASE calculator, see the examples/ folder.

## Links

- libMBD documentation: https://libmbd.nablapsi.science
- pyMBD documentation: https://libmbd.nablapsi.science/pymbd

## Developing

For development, a top-level `Makefile` is included, which configures and compiles libMBD, compiles the pyMBD C extension, and runs both libMBD and pyMBD tests.

```
git clone https://github.com/libmbd/libmbd.git && cd libmbd
python3 -m venv venv && source venv/bin/activate
make
# development work...
make
```
