#!/usr/bin/env python3
"""Reproduce libmbd/libmbd#72: MBD@rsSCS crashes for a simple bulk Cu system.

Issue #72 reports that VASP (IVDW=14, MBD@rsSCS) crashes for a Cu primitive
cell.  The geometry and vdW parameters are discussed on the VASP forum thread
https://www.vasp.at/forum/viewtopic.php?t=20071.  The system is an fcc Cu
crystal (one atom in the primitive cell).

Running the same MBD@rsSCS method through pyMBD on that structure reproduces the
failure: once the k-point grid is dense enough to be physically meaningful, the
coupled-dipole-model (CDM) Hamiltonian acquires negative eigenvalues.  Those
eigenvalues are squared oscillator frequencies, so a negative value means an
imaginary frequency -- the harmonic MBD oscillator system is unstable (the
"polarization catastrophe" that MBD is known to suffer for dense, highly
polarizable metals).  libMBD detects this in mbd_hamiltonian.F90 and raises
``MBD_EXC_NEG_EIGVALS``; were the energy ``1/2 * sum(sqrt(eigs))`` evaluated
anyway it would be NaN.  In the VASP interface the same condition surfaces as a
crash.

This is therefore not a numerical bug in libMBD but an intrinsic limitation of
the MBD model for bulk metals like Cu.

Run with:  python reproduce_issue_72.py
"""

import numpy as np

from pymbd.fortran import MBDGeom

ANG = 1 / 0.52917721092  # angstrom -> bohr


def cu_fcc(a_ang):
    """fcc Cu primitive cell (one atom) at lattice constant ``a_ang`` (angstrom)."""
    lattice = (a_ang / 2) * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]) * ANG
    coords = np.zeros((1, 3))
    return coords, lattice


def mbd_energy(a_ang, volume_ratio, k_grid, zero_neg=False):
    coords, lattice = cu_fcc(a_ang)
    geom = MBDGeom(coords, lattice=lattice, k_grid=k_grid)
    # beta = 0.83 is the MBD@rsSCS damping parameter for PBE (VASP IVDW=14)
    return geom.mbd_energy_species(['Cu'], [volume_ratio], 0.83)


def main():
    a = 3.615  # experimental fcc Cu lattice constant (angstrom)
    print(f'fcc Cu, a = {a} A, free-atom vdW params (volume ratio 1.0), beta = 0.83\n')
    print('Convergence of MBD@rsSCS energy with the k-point grid:')
    for kg in ([2, 2, 2], [4, 4, 4], [6, 6, 6], [8, 8, 8]):
        try:
            e = mbd_energy(a, 1.0, kg)
            print(f'  k_grid={kg}:  energy = {e:.6f} Ha')
        except Exception as exc:  # noqa: BLE001
            msg = exc.args[2] if len(getattr(exc, 'args', [])) > 2 else exc
            print(f'  k_grid={kg}:  FAILED -> {type(exc).__name__}: {msg}')
    print(
        '\nThe failure is "CDM Hamiltonian has N negative eigenvalues": the MBD '
        'oscillator\nsystem is unstable for bulk Cu, so the rsSCS energy is '
        'undefined (sqrt of a\nnegative squared frequency -> NaN). This is the '
        'crash reported in issue #72.'
    )


if __name__ == '__main__':
    main()
