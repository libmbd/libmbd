#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""Plot the S66x8 deviation, and what it is worth, against separation.

The relative deviation of our MBD term from the FHI-aims one per system, over
the MARE against CCSD(T)/CBS with our MBD term, with the reference's, and with
no dispersion. Being an offset rather than scatter, the deviation cancels
against the reference's own error at some separations and compounds at others,
which the summary tables can state but only this shows.

Usage:
    python pyscf_s66x8_plot.py --results-dir all-results [--out PATH]
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt  # noqa: E402  -- after the backend is fixed

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pyscf_s66x8_aggregate import load_shards
from pyscf_s66x8_figshare_check import flatten_rows, summarize


def plot(systems, meta, out_path):
    """Draw both panels and write the figure."""
    rows = flatten_rows(systems.values())
    overall = summarize(rows)
    scales = [b['scale'] for b in overall['by_scale']]

    fig, (top, bottom) = plt.subplots(
        2, 1, sharex=True, figsize=(7, 7), constrained_layout=True
    )
    for system in systems.values():
        if len(system['points']) < 2:
            continue
        by_scale = summarize(system['points'])['by_scale']
        top.plot(
            [b['scale'] for b in by_scale],
            [b['mbd_rel_dev'] for b in by_scale],
            color='0.7',
            lw=0.8,
            zorder=1,
        )
    top.plot(
        scales,
        [b['mbd_rel_dev'] for b in overall['by_scale']],
        'o-',
        color='C0',
        lw=2,
        zorder=3,
        label=f'mean over {len(overall["systems"])} systems',
    )
    top.axhline(0, color='k', lw=0.8, zorder=2)
    top.set_ylabel('MBD deviation from FHI-aims (%)')
    top.set_title(
        f'S66x8, {meta["basis"]} / {meta["xc"]} / grid {meta["grid"]} / '
        f'conv {meta["conv_tol"]:g}'
    )
    top.legend(loc='best', fontsize='small')

    for key, label, style in (
        ('mare_nodisp', 'no dispersion', ':'),
        ('mare_ref', 'reference MBD', '--'),
        ('mare_recon', 'our MBD', '-'),
    ):
        bottom.plot(
            scales,
            [b[key] for b in overall['by_scale']],
            style,
            marker='o',
            ms=3,
            label=label,
        )
    bottom.set_xlabel('separation (fraction of the equilibrium distance)')
    bottom.set_ylabel('MARE vs CCSD(T)/CBS (%)')
    bottom.legend(loc='best', fontsize='small')

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f'wrote {out_path}')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-dir', required=True, help='directory of shard JSONs')
    p.add_argument('--out', default='results/convergence-s66x8.png')
    args = p.parse_args(argv)

    # the workflow runs this unconditionally, so an empty or failed run is a
    # no-op rather than a second red step on top of the aggregator's
    try:
        systems, meta = load_shards(args.results_dir)
    except SystemExit as exc:
        print(f'nothing to plot: {exc}')
        return 0
    if not flatten_rows(systems.values()):
        print('nothing to plot: every system failed')
        return 0
    plot(systems, meta, args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
