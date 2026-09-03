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

Besides the MBD term itself the summary reports what that deviation is worth on
the benchmark the method is judged by, as a MARE against the CCSD(T)/CBS S66x8
reference energies. Both the reconstructed and the reference total take the DFT
part from the reference's own all-electron PBE interaction energies, so the DFT
is identical in the two and the whole gap between them is the dispersion term.
Using our PBE instead would drown it, being uncorrected for BSSE and in a finite
basis. Both are broken down by separation alongside the relative MBD deviation,
since that deviation is very nearly a pure systematic offset and so cancels
against the reference's own error at some separations and compounds at others --
which means a basis can rank well on the compressed geometries for reasons that
reverse at long range.

The defaults are all chosen for the MBD term rather than for the DFT:

  * ``--basis aug-cc-pvdz``. What matters is diffuse functions rather than zeta
    level, since the r^3 weight of the volume integrals probes the density tail;
    an augmented double-zeta set beats an unaugmented triple-zeta one. The
    augmentation has to be even-handed across elements, because the ratio is
    taken against a promolecule -- the ma-def2 sets augment the heavy atoms only,
    which shifts weight off hydrogen and reallocates the error between the
    hydrogen-bonded and pi-stacked systems rather than removing it.
  * ``--grid 1``, against PySCF's 3. The same grid integrates the in-molecule and
    the free-atom volume, so its error largely cancels in the ratio -- but not
    unconditionally: hydrogen carries the coarsest atomic grid, so the aliphatic
    systems set the requirement, and level 0 does not meet it.
  * ``--conv-tol 1e-4``, against PySCF's 1e-9, resting on the same cancellation.
    Loosening it further turns a negligible shift into scatter comparable with
    the deviation being measured, and scatter, unlike a bias, does not average
    out.

Successive points of a curve start from the previous point's converged density,
which combines with the loose threshold without the error compounding along the
curve.

The DFT (PBE) term is reported alongside but is not a bridge diagnostic: with no
counterpoise correction its residual against the all-electron reference is basis
set incompleteness plus BSSE plus grid. The augmented default makes it worse
rather than better, since diffuse functions inflate BSSE. Raise --grid and
--basis before reading anything into it.

Data (downloaded/located, not shipped):
  * figshare all-data.h5  https://ndownloader.figshare.com/files/9775933
    (article 5117167, "dft-vdw-range Data")
  * geometries + labels from the vdwsets data files (read directly, no import)
    https://github.com/azag0/vdwsets  (git clone it and pass --vdwsets)

Requires ``pyscf``, ``pandas``, ``tables`` and a pyMBD with the compiled
extension.

Usage:
    OPENBLAS_NUM_THREADS=1 python pyscf_s66x8_figshare_check.py \\
        --vdwsets /path/to/vdwsets/clone \\
        [--h5 all-data.h5] [--idx 1 2 8 32 33 51 59 60] \\
        [--basis aug-cc-pvdz] [--out-json results/s66x8-b0.json]

Almost all the run time is the DFT itself, so it is worth giving OpenBLAS a
single thread and leaving the rest to PySCF's OpenMP: where the two thread pools
both try to use every core they contend badly. Restricting OpenMP instead is not
a substitute.
"""

import argparse
import csv
import json
import os
import re
import sys
import time
import traceback
import urllib.request
from importlib.metadata import PackageNotFoundError, version
from itertools import groupby
from operator import itemgetter
from pathlib import Path

import numpy as np

AU2KCAL = 627.509474
H5_URL = 'https://ndownloader.figshare.com/files/9775933'
# MBD@rsSCS for PBE (matches the reference calculations)
MBD_A, MBD_BETA = 6.0, 0.83
# separations per system in S66x8; a curve with fewer points ran short
N_SCALES = 8
# version of the --out-json payload, checked by the aggregator
SCHEMA = 1


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
        # renamed on success, so a half-finished download is not mistaken for a
        # good one by the next run
        part = path.with_name(path.name + '.part')
        urllib.request.urlretrieve(H5_URL, part)
        os.replace(part, path)
    return path


def load_paper(h5_path):
    """Return the reference values as {'pbe'|'mbd'|'ref': {(label, scale): value}}."""
    import pandas as pd

    with pd.HDFStore(str(h5_path), 'r') as store:
        mbd = store['mbd'].xs(('S66x8', MBD_A, MBD_BETA), level=('name', 'a', 'beta'))[
            'ene'
        ]
        scf = store['scf'].xs(('S66x8', 'pbe', False), level=('name', 'xc', 'cp'))
    pbe = {(_norm(s), float(d)): v for (s, d), v in scf['ene'].items()}
    ref = {(_norm(s), float(d)): v for (s, d), v in scf['ref'].items()}
    mbd = {(_norm(s), float(d)): v for (s, d), v in mbd.items()}
    return {'pbe': pbe, 'mbd': mbd, 'ref': ref}


_ENERGIES = {}


def our_energies(species, coords, basis, xc, grid, conv_tol, dm0=None):
    """Return (SCF energy, MBD energy, density matrix) for one fragment.

    Memoized on the geometry: the two monomers of a system are taken at its
    equilibrium separation whatever the separation of the complex, so they are the
    same fragment at all eight points of a dissociation curve and only need
    computing once. The density matrix is returned for chaining (below) and is
    None on a cache hit, since it is not worth keeping every fragment's -- which
    doubles as how the caller counts the SCFs it actually ran.

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


def run_curve(points, refs, args):
    """Return the rows for one system's dissociation curve.

    The curve is the unit of work because the density chaining is curve-local,
    which is also what lets a failing system be dropped without poisoning the
    next one.
    """
    rows, dm = [], None
    for idx, label, scale, frags in points:
        key = (_norm(label), float(scale))
        if key not in refs['mbd']:
            continue
        t0 = time.perf_counter()
        e, n_scf = {}, 0
        for name, g in frags.items():
            guess = dm if name == 'complex' else None
            *e[name], new_dm = our_energies(
                *g, args.basis, args.xc, args.grid, args.conv_tol, guess
            )
            n_scf += new_dm is not None  # None means served from _ENERGIES
            if name == 'complex':
                dm = new_dm
        pbe_int = (e['complex'][0] - e['fragment-1'][0] - e['fragment-2'][0]) * AU2KCAL
        mbd_int = (e['complex'][1] - e['fragment-1'][1] - e['fragment-2'][1]) * AU2KCAL
        rows.append(
            {
                'idx': idx,
                'label': label,
                'scale': float(scale),
                'mbd_our': mbd_int,
                'mbd_paper': refs['mbd'][key],
                'pbe_our': pbe_int,
                'pbe_paper': refs['pbe'][key],
                'ref': refs['ref'][key],
                'wall': time.perf_counter() - t0,
                'n_scf': n_scf,
            }
        )
        r = rows[-1]
        print(
            f'{idx:2d} {label:28.28s} {scale:4.2f} | '
            f'MBD our {r["mbd_our"]:7.3f}  paper {r["mbd_paper"]:7.3f}  '
            f'Δ {r["mbd_our"] - r["mbd_paper"]:+.3f}',
            flush=True,
        )
    return rows


def run_systems(refs, args):
    """Run every requested system, one dissociation curve at a time.

    A system that raises is recorded and skipped rather than taking the rest of
    the shard down with it.
    """
    idx_filter = None if args.all else set(args.idx)
    systems, done = [], set()
    for idx, group in groupby(load_geometries(args.vdwsets, idx_filter), itemgetter(0)):
        if idx in done:
            raise SystemExit(f'system {idx} is not contiguous in energies.csv')
        done.add(idx)
        group = list(group)
        t0 = time.perf_counter()
        try:
            rows, error = run_curve(group, refs, args), None
        except Exception as exc:
            rows, error = [], f'{type(exc).__name__}: {exc}'
            traceback.print_exc()
        systems.append(
            {
                'idx': idx,
                'label': group[0][1],
                'error': error,
                'wall': time.perf_counter() - t0,
                'n_scf': sum(r['n_scf'] for r in rows),
                'points': rows,
            }
        )
    return systems


def flatten_rows(systems):
    """Every system's point rows, in the shape summarize/format_report want."""
    return [p for s in systems for p in s['points']]


def _version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def meta_from_args(args):
    """The settings a result was produced under, for the JSON and the report header."""
    return {
        'basis': args.basis,
        'xc': args.xc,
        'grid': args.grid,
        'conv_tol': args.conv_tol,
        'mbd_a': MBD_A,
        'mbd_beta': MBD_BETA,
        'versions': {n: _version(n) for n in ('pymbd', 'pyscf', 'numpy')},
    }


def summarize(rows):
    """Aggregate statistics over point rows.

    Shared by the report, the CI aggregator (which calls it per system too, to
    gate on ``mbd_rel_dev``) and the plot script.
    """
    mbd_our = np.array([r['mbd_our'] for r in rows])
    mbd_paper = np.array([r['mbd_paper'] for r in rows])
    pbe_paper = np.array([r['pbe_paper'] for r in rows])
    dpbe = np.array([r['pbe_our'] - r['pbe_paper'] for r in rows])
    ref = np.array([r['ref'] for r in rows])
    scales = np.array([r['scale'] for r in rows])
    dmbd = mbd_our - mbd_paper
    argmax = int(np.abs(dmbd).argmax())

    def mare(total, sel=slice(None)):
        return 100 * np.mean(np.abs((total[sel] - ref[sel]) / ref[sel]))

    def rel_dev(sel=slice(None)):
        return 100 * np.mean(mbd_our[sel] / mbd_paper[sel] - 1)

    return {
        'n': len(rows),
        'systems': sorted({r['idx'] for r in rows}),
        'wall': sum(r['wall'] for r in rows),
        'mbd_mae': float(np.abs(dmbd).mean()),
        'mbd_absmax': float(np.abs(dmbd).max()),
        'mbd_argmax': (rows[argmax]['idx'], rows[argmax]['scale']),
        'mbd_bias': float(dmbd.mean()),
        'mbd_rel_dev': float(rel_dev()),
        'pbe_mae': float(np.abs(dpbe).mean()),
        'pbe_bias': float(dpbe.mean()),
        'mare_recon': float(mare(pbe_paper + mbd_our)),
        'mare_ref': float(mare(pbe_paper + mbd_paper)),
        'mare_nodisp': float(mare(pbe_paper)),
        'by_scale': [
            {
                'scale': float(s),
                'n': int((scales == s).sum()),
                'mbd_rel_dev': float(rel_dev(scales == s)),
                'mare_recon': float(mare(pbe_paper + mbd_our, scales == s)),
                'mare_ref': float(mare(pbe_paper + mbd_paper, scales == s)),
                'mare_nodisp': float(mare(pbe_paper, scales == s)),
                'abs_ref': float(np.mean(np.abs(ref[scales == s]))),
            }
            for s in sorted(set(scales))
        ],
    }


def format_report(rows, seen, meta, failed=()):
    """The summary block: MBD deviation, reconstructed MARE, and both by separation.

    The reconstructed total energy is our MBD term on top of the *reference's* own
    all-electron PBE interaction energy. That keeps the DFT part identical to the
    reference, so the whole difference between the reconstructed and the reference
    MARE is the dispersion term, i.e. the bridge. Adding our own PBE instead would
    swamp it, that term being uncorrected for BSSE and in a finite basis.
    """
    s = summarize(rows)
    out = [
        f'\n{s["n"]} points over systems {sorted(seen)}',
        f'MBD  (ours vs FHI-aims):  MAE {s["mbd_mae"]:.4f}  '
        f'max|Δ| {s["mbd_absmax"]:.4f}  bias {s["mbd_bias"]:+.4f} kcal/mol  '
        f'({s["mbd_rel_dev"]:+.2f}% mean relative)',
        f'PBE  (our {meta["basis"]}, grid {meta["grid"]}, '
        f'conv {meta["conv_tol"]:g}, no CP, vs all-electron):  '
        f'MAE {s["pbe_mae"]:.4f}  bias {s["pbe_bias"]:+.4f} kcal/mol  '
        f'(basis set + BSSE + grid, not the bridge)',
        # against the CCSD(T)/CBS benchmark, with the reference's PBE in both totals
        f'\nvs CCSD(T)/CBS, reference PBE + MBD:  '
        f'reconstructed (our MBD) MARE {s["mare_recon"]:.1f}%   '
        f'reference MARE {s["mare_ref"]:.1f}%   '
        f'no dispersion {s["mare_nodisp"]:.1f}%',
        '\nby separation:',
        f'{"scale":>7} {"n":>5} {"MBD rel dev":>13} '
        f'{"MARE recon":>12} {"MARE ref":>10} {"|E| ref":>9}',
    ]
    out += [
        f'{b["scale"]:>7.2f} {b["n"]:>5d} '
        f'{b["mbd_rel_dev"]:>+12.2f}% '
        f'{b["mare_recon"]:>11.1f}% '
        f'{b["mare_ref"]:>9.1f}% '
        f'{b["abs_ref"]:>9.3f}'
        for b in s['by_scale']
    ]
    if failed:
        out.append('\nfailed systems:')
        out += [f'{f["idx"]:2d} {f["label"]:28.28s} {f["error"]}' for f in failed]
    return '\n'.join(out) + '\n'


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vdwsets', required=True, help='path to a vdwsets checkout')
    p.add_argument('--h5', default='all-data.h5', help='path to all-data.h5')
    p.add_argument('--basis', default='aug-cc-pvdz')
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
    p.add_argument(
        '--out-json',
        metavar='PATH',
        help='write the per-point results here, for pyscf_s66x8_aggregate.py',
    )
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    systems = run_systems(load_paper(ensure_h5(Path(args.h5))), args)
    # written before the exit code is decided, so a shard that lost a system
    # still hands the rest over to the aggregator
    if args.out_json:
        Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_json, 'w') as f:
            json.dump(
                {'schema': SCHEMA, 'meta': meta_from_args(args), 'systems': systems},
                f,
                indent=1,
            )

    rows = flatten_rows(systems)
    failed = [s for s in systems if s['error']]
    if not rows:
        print('no matching systems found', file=sys.stderr)
        for s in failed:
            print(f'{s["idx"]:2d} {s["label"]}: {s["error"]}', file=sys.stderr)
        return 1
    print(
        format_report(
            rows,
            {s['idx'] for s in systems if s['points']},
            meta_from_args(args),
            failed,
        ),
        end='',
    )
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
