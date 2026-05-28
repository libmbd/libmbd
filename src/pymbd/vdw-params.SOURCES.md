# Provenance of `vdw-params.csv`

The file `vdw-params.csv` is a compilation of free-atom static dipole
polarizabilities (`alpha_0`, in a.u.), dipole–dipole dispersion coefficients
(`C6`, in a.u.) and van der Waals radii (`R_vdw`, in Å) used by the
Tkatchenko–Scheffler (TS) and many-body dispersion (MBD) methods.

The compiled table itself is dedicated to the public domain via CC0-1.0
(see `vdw-params.csv.license`). The underlying numerical values originate
from the published sources below.

## `alpha_0(BG)`, `C6(BG)` and all ion rows (symbols containing `+` or `-`)

Tim Gould and Tomáš Bučko,
*"C6 coefficients and dipole polarizabilities for all atoms and many ions
in rows 1–6 of the periodic table"*,
**J. Chem. Theory Comput.** 12, 3603–3613 (2016).
DOI: [10.1021/acs.jctc.6b00361](https://doi.org/10.1021/acs.jctc.6b00361)

These values were added to libMBD in commit
[`471020e`](https://github.com/libmbd/libmbd/commit/471020e882a775d5164315cce62ef70cb7d53bd5)
("add vdw parameters from gould-bucko").

## `alpha_0(TS)`, `C6(TS)`, `R_vdw(TS)` for Z = 1–86 (excluding lanthanides)

Original TS reference values as established in

Alexandre Tkatchenko and Matthias Scheffler,
*"Accurate Molecular Van Der Waals Interactions from Ground-State Electron
Density and Free-Atom Reference Data"*,
**Phys. Rev. Lett.** 102, 073005 (2009).
DOI: [10.1103/PhysRevLett.102.073005](https://doi.org/10.1103/PhysRevLett.102.073005)

These values were inherited from the original `pymbd/vdw_param.py` table
(public-domain via CC0-1.0) when the file was reorganised into CSV form in
commit
[`8c31758`](https://github.com/libmbd/libmbd/commit/8c31758fe79cdd60756aea0c9614eb200e31a67e).

## `alpha_0(TS)`, `C6(TS)`, `R_vdw(TS)` for La–Lu and Fr–No

Added in commit
[`3c6900c`](https://github.com/libmbd/libmbd/commit/3c6900c058748cd7145a219dbb8e766fcddac2ca).
Sourced from a comprehensive periodic-table-wide TS-style tabulation
(Tomáš Bučko, habilitation thesis, TU Berlin, Table A.1;
<https://depositonce.tu-berlin.de/items/40bff578-b71a-4ed6-b7c6-667226e7af5c>).

## `alpha_0(TSsurf)`, `C6(TSsurf)`, `R_vdw(TSsurf)`

Ni, Cu, Pd, Ag, Pt, Au:

Victor G. Ruiz, Wei Liu, Egbert Zojer, Matthias Scheffler and Alexandre
Tkatchenko,
*"Density-Functional Theory with Screened van der Waals Interactions for
the Modeling of Hybrid Inorganic–Organic Systems"*,
**Phys. Rev. Lett.** 108, 146103 (2012).
DOI: [10.1103/PhysRevLett.108.146103](https://doi.org/10.1103/PhysRevLett.108.146103)

Zn (and consistent extensions):

Reinhard J. Maurer, Victor G. Ruiz, Javier Camarillo-Cisneros, Wei Liu,
Nicola Ferri, Karsten Reuter and Alexandre Tkatchenko,
*"Adsorption structures and energetics of molecules on metal surfaces:
Benchmarking DFT-based vdW methods"*,
**Phys. Rev. B** 93, 035118 (2016).
DOI: [10.1103/PhysRevB.93.035118](https://doi.org/10.1103/PhysRevB.93.035118)

Added in commit
[`c195111`](https://github.com/libmbd/libmbd/commit/c19511127e86404b8cebcc10b06d664c54812b16).
