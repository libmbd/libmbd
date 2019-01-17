---
project: Libmbd
summary: Many-body dispersion library
license: by
src_dir: ../src
css: tweaks.css
hide_undoc: true
preprocessor: gfortran -cpp -E -P -DWITH_MPI -DWITH_SCALAPACK
exclude:
    mbd_blacs.f90
    mbd_c_api.F90
    mbd_coulomb.f90
    mbd_density.f90
    mbd_lapack.f90
    mbd_linalg.F90
    mbd_matrix.F90
    mbd_mpi.F90
    mbd_rpa.F90
    mbd_scalapack.f90
    mbd_vdw_param.f90
---

At the moment the documentation consists of an automatically generated API reference and a miniature of a user guide in the following paragraph. All mathematical formulas used in the code are documented directly in the source code and rendered in [Procedures](lists/procedures.html). Installation instructions can be found in the [Readme](https://github.com/jhrmnn/libmbd/blob/master/README.md).

The user-facing Fortran API of Libmbd is contained in the [[mbd]] module and consists of the [[mbd_input_t]] and [[mbd_calc_t]] derived types. A [[mbd_input_t]] object serves to set various options for the calculation and is used to initialize a [[mbd_calc_t]] object, which is then used to actually perform the MBD calculation.

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
    stop
end if
call calc%update_vdw_params_from_ratios([0.98d0, 0.98d0])
call calc%evaluate_vdw_method(energy)
call calc%get_gradients(gradients)
call calc%destroy()
```
