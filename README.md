**!! functional but experimental !!**

Implementation of the MBD method in Fortran, ported to Python via `f2py`.

To compile, first adapt `system.example.mk` and save it as `system.mk` (which is sourced by `Makefile`) to suit your system, and make sure you have `f2py` installed. Then, run `make`. The build process also runs a test at the end.

```
>>> import pymbd
>>> for key, val in pymbd.main('mbddata.json').items():
...     print key, val
... 
Running on 1 nodes...
Evaluating TS...
Evaluating MBD@TS...
Evaluating SCS...
Evaluating TS@SCS...
Evaluating MBD@SCS...
Evaluating rsSCS...
Evaluating MBD@rsSCS...
MBD(TS)@rsSCS~fermi@rsSCS,dip -0.000379972055702
MBD(RPA)@rsSCS~fermi@rsSCS,dip [ -3.80232137e-04  -3.79972056e-04  -3.07640882e-23  -2.59715992e-07
   1.32553274e-25  -3.64538951e-10   0.00000000e+00  -6.63334333e-13
   0.00000000e+00  -1.36246806e-15]
MBD@TS~fermi@TS,dip -0.00044416943627
MBD(nbody)@rsSCS~fermi@rsSCS,dip [-0.00038023 -0.00038023  0.        ]
MBD@TS~erf@TS,dip -0.000222862762834
TS@SCS~fermi@SCS -8.86147033383e-05
TS@TS~fermi@TS -9.13517851601e-05
MBD@SCS~dip,1mexp@SCS -0.000246993634594
MBD@rsSCS~fermi@rsSCS,dip -0.000380232170903

```
