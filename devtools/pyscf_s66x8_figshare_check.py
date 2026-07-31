#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""Validate the pymbd.pyscf bridge against published all-electron MBD data.

For a subset of S66x8 dimers this recomputes the MBD@rsSCS interaction energy
with the PySCF -> pyMBD bridge (:mod:`pymbd.pyscf`) and compares it, point by
point, with the FHI-aims MBD energies of

    Hermann & Tkatchenko, "Electronic exchange and correlation in van der Waals
    systems", J. Chem. Theory Comput. (2018); doi:10.1021/acs.jctc.7b01172,

whose processed results are published on figshare. Because FHI-aims shares no
code, basis, or Hirshfeld implementation with PySCF/pyMBD, agreement on the MBD
term is an independent cross-code validation of the bridge.

Use a basis of at least triple-zeta quality. Taking the free-atom reference in
the molecular basis makes the Hirshfeld ratios far less basis-sensitive than the
density itself, but not insensitive: over the default subset def2-SVP gives an
MBD MAE of 0.037 kcal/mol against def2-TZVP's 0.018, and the gap is concentrated
in the systems with pi systems (0.11 for pyridine-ethene, 0.09 for
uracil-ethyne, against 0.005-0.009 for the small hydrogen-bonded ones).

The integration grid matters much less than it does for the DFT, because the
same grid integrates the in-molecule and the free-atom volume and its error
largely cancels in their ratio, so this script defaults to level 1 rather than
PySCF's 3. It does not survive level 0, though, and how far it survives is
system-dependent: on an aromatic dimer levels 0 to 4 span 7e-4 kcal/mol, but on
pentane-pentane, the most hydrogen-rich system in S66, level 0 is 0.029 off
while level 1 is within 3e-4 of a converged grid. Hydrogen carries the coarsest
atomic grid, so the aliphatic systems set the requirement. The DFT (PBE) term is
reported alongside, but at this grid its residual against the all-electron
reference has a grid component, on top of the basis set and (no counterpoise)
BSSE -- raise --grid before reading anything into it. None of it is a bridge
effect.

The SCF convergence threshold can be loosened a long way for the same reason, so
it defaults to 1e-4 rather than PySCF's 1e-9. Measured over the default subset
against a tight run, that shifts the MBD energies by 0.0002 kcal/mol on average
and 0.002 at worst, and takes the run from 220 s to 139 s. Going on to 1e-2
reaches 123 s but the scatter grows to 0.0035 and 0.012, which is an appreciable
fraction of the deviation being measured, and being scatter rather than bias it
does not average out.

Successive points of a curve start from the previous point's converged density,
which combines with the loose threshold without the error compounding along the
curve: the deviation from a tight run is flat from the first point, which has no
guess, to the last.

Data (downloaded/located, not shipped):
  * figshare all-data.h5  https://ndownloader.figshare.com/files/9775933
    (article 5117167, "dft-vdw-range Data")
  * geometries + labels from the vdwsets data files (read directly, no import)
    https://github.com/azag0/vdwsets  (git clone it and pass --vdwsets)

Usage:
    OPENBLAS_NUM_THREADS=1 python pyscf_s66x8_figshare_check.py \\
        --vdwsets /path/to/vdwsets/clone \\
        [--h5 all-data.h5] [--idx 1 2 8 32 33 51 59 60] [--basis def2-tzvp]

Almost all the run time is the DFT itself, so it is worth giving OpenBLAS a
single thread and leaving the rest to PySCF's OpenMP. Where the two thread pools
both try to use every core they contend badly: a water dimer took 1.28 s with
both at four threads against 0.15 s with OPENBLAS_NUM_THREADS=1, and a 17-atom
dimer 13.4 s against 9.9 s. Restricting OpenMP instead is not a substitute --
two threads was slower than four at both sizes.
"""

import argparse
import csv
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np

AU2KCAL = 627.509474
H5_URL = 'https://ndownloader.figshare.com/files/9775933'
# MBD@rsSCS for PBE (matches the reference calculations)
MBD_A, MBD_BETA = 6.0, 0.83


def _norm(label):
    return re.sub(r'\s+', ' ', label).strip().lower()


def read_xyz(path):
    species, coords = [], []
    with open(path) as f:
        n = int(f.readline())
        f.readline()  # comment
        for _ in range(n):
            el, x, y, z = f.readline().split()[:4]
            species.append(el)
            coords.append([float(x), float(y), float(z)])
    return species, coords


def load_geometries(vdwsets_root, idx_filter):
    """Yield (idx, label, scale, fragments) from a vdwsets checkout's data files.

    Reproduces the S66x8 filename convention of vdwsets.get_s66x8 without
    importing the package (which needs pkg_resources). ``fragments`` maps
    'complex'/'fragment-1'/'fragment-2' to (species, coords).
    """
    data = Path(vdwsets_root) / 'vdwsets' / 'data' / 's66x8'
    with open(data / 'energies.csv') as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        idx, label, scale = int(row['idx']), row['label'], float(row['scale'])
        if idx_filter is not None and idx not in idx_filter:
            continue
        file_lbl = label.replace(' ', '-').lower()
        frags = {}
        for j, frag in enumerate(['complex', 'fragment', 'fragment']):
            scale_lbl = scale if frag == 'complex' else 1.0
            path = (
                data / 'geoms' / f'{idx:02}-{scale_lbl:.2f}_{file_lbl}_{frag}_{j}.xyz'
            )
            key = 'complex' if frag == 'complex' else f'fragment-{j}'
            frags[key] = read_xyz(path)
        yield idx, label, scale, frags


def ensure_h5(path):
    if not path.exists():
        print(f'downloading all-data.h5 -> {path} ...', flush=True)
        urllib.request.urlretrieve(H5_URL, path)
    return path


def load_paper(h5_path):
    import pandas as pd

    with pd.HDFStore(str(h5_path), 'r') as store:
        mbd = store['mbd'].xs(('S66x8', MBD_A, MBD_BETA), level=('name', 'a', 'beta'))[
            'ene'
        ]
        scf = store['scf'].xs(('S66x8', 'pbe', False), level=('name', 'xc', 'cp'))
    pbe = {(_norm(s), float(d)): v for (s, d), v in scf['ene'].items()}
    ref = {(_norm(s), float(d)): v for (s, d), v in scf['ref'].items()}
    mbd = {(_norm(s), float(d)): v for (s, d), v in mbd.items()}
    return pbe, mbd, ref


_ENERGIES = {}


def our_energies(species, coords, basis, xc, grid, conv_tol, dm0=None):
    """Return (SCF energy, MBD energy, density matrix) for one fragment.

    Memoized on the geometry: the two monomers of a system are taken at its
    equilibrium separation whatever the separation of the complex, so they are the
    same fragment at all eight points of a dissociation curve and only need
    computing once. The density matrix is returned for chaining (below) and is
    None on a cache hit, since it is not worth keeping every fragment's.

    ``dm0`` is an initial guess, meant to be the converged density of the previous
    point of a dissociation curve. The monomers are rigid along one, so successive
    complexes differ only in their separation and the previous density is a good
    starting point. At a tight threshold it cannot change the converged result; at
    a loose one it can shift where the SCF stops, but only within that threshold's
    own noise, and measurably without compounding along the curve. Either way it
    is not part of the cache key, which identifies a geometry rather than a path
    taken to it.
    """
    key = (basis, xc, grid, conv_tol, tuple(species), tuple(map(tuple, coords)))
    if key in _ENERGIES:
        return (*_ENERGIES[key], None)

    from pyscf import dft, gto

    from pymbd.pyscf import mbd_energy

    mol = gto.M(
        atom=[[s, tuple(c)] for s, c in zip(species, coords)],
        basis=basis,
        unit='Angstrom',
        verbose=0,
    )
    mf = dft.RKS(mol, xc=xc).density_fit()
    mf.grids.level = grid
    mf.conv_tol = conv_tol
    mf.kernel(dm0=dm0)
    _ENERGIES[key] = (mf.e_tot, mbd_energy(mf, beta=MBD_BETA))
    return (*_ENERGIES[key], mf.make_rdm1())


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vdwsets', required=True, help='path to a vdwsets checkout')
    p.add_argument('--h5', default='all-data.h5', help='path to all-data.h5')
    p.add_argument('--basis', default='def2-tzvp')
    p.add_argument('--xc', default='PBE')
    p.add_argument(
        '--grid',
        type=int,
        default=1,
        help='PySCF grid level (default 1: enough for MBD, see the module '
        'docstring, but too coarse for the DFT energies)',
    )
    p.add_argument(
        '--conv-tol',
        type=float,
        default=1e-4,
        help='SCF convergence threshold (default 1e-4: see the module docstring; '
        "PySCF's own default is 1e-9)",
    )
    p.add_argument(
        '--idx',
        type=int,
        nargs='+',
        default=[1, 2, 8, 32, 33, 51, 59, 60],
        help='S66 system indices to run (default: a small E/D/M-spanning subset)',
    )
    p.add_argument('--all', action='store_true', help='run all 66 systems')
    args = p.parse_args(argv)

    pbe_ref, mbd_ref, ccsdt = load_paper(ensure_h5(Path(args.h5)))
    idx_filter = None if args.all else set(args.idx)

    rows = []
    seen = set()
    # the converged density of the previous point of the current dissociation
    # curve, fed to the next one as its initial guess; only one is ever held,
    # since the points arrive grouped by system
    curve_idx, curve_dm = None, None
    for idx, label, scale, frags in load_geometries(args.vdwsets, idx_filter):
        key = (_norm(label), float(scale))
        if key not in mbd_ref:
            continue
        if idx != curve_idx:
            curve_idx, curve_dm = idx, None
        e = {}
        for name, g in frags.items():
            guess = curve_dm if name == 'complex' else None
            *e[name], dm = our_energies(
                *g, args.basis, args.xc, args.grid, args.conv_tol, guess
            )
            if name == 'complex':
                curve_dm = dm
        pbe_int = (e['complex'][0] - e['fragment-1'][0] - e['fragment-2'][0]) * AU2KCAL
        mbd_int = (e['complex'][1] - e['fragment-1'][1] - e['fragment-2'][1]) * AU2KCAL
        rows.append(
            {
                'idx': idx,
                'label': label,
                'scale': float(scale),
                'mbd_our': mbd_int,
                'mbd_paper': mbd_ref[key],
                'pbe_our': pbe_int,
                'pbe_paper': pbe_ref[key],
                'ref': ccsdt[key],
            }
        )
        seen.add(idx)
        r = rows[-1]
        print(
            f'{idx:2d} {label:28.28s} {scale:4.2f} | '
            f'MBD our {r["mbd_our"]:7.3f}  paper {r["mbd_paper"]:7.3f}  '
            f'Δ {r["mbd_our"] - r["mbd_paper"]:+.3f}',
            flush=True,
        )

    if not rows:
        print('no matching systems found', file=sys.stderr)
        return 1
    dmbd = np.array([r['mbd_our'] - r['mbd_paper'] for r in rows])
    dpbe = np.array([r['pbe_our'] - r['pbe_paper'] for r in rows])
    print(f'\n{len(rows)} points over systems {sorted(seen)}')
    print(
        f'MBD  (ours vs FHI-aims):  MAE {np.abs(dmbd).mean():.4f}  '
        f'max|Δ| {np.abs(dmbd).max():.4f}  bias {dmbd.mean():+.4f} kcal/mol'
    )
    print(
        f'PBE  (our {args.basis}, grid {args.grid}, conv {args.conv_tol:g}, no CP, '
        f'vs all-electron):  '
        f'MAE {np.abs(dpbe).mean():.4f}  bias {dpbe.mean():+.4f} kcal/mol  '
        f'(basis set + BSSE + grid, not the bridge)'
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
