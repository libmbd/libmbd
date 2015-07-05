**!! functional but experimental !!**

Implementation of the MBD method in Fortran, ported to Python via `f2py`.

To compile, first adapt `system.example.mk` and save it as `system.mk` (which is sourced by `Makefile`) to suit your system, and make sure you have `f2py` installed. Then, run `make`. The build process also runs a test at the end.
