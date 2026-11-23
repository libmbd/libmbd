#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""Merge the S66x8 benchmark shards, gate on the committed baseline, and report.

The gate compares each system against what a previous run of the same setup
measured, in pyscf_s66x8_reference.json, not against the FHI-aims values
themselves: the deviation from those is basis-set incompleteness in the
Hirshfeld ratios and is not meant to be zero. The baseline therefore has to be
seeded from a full run before it gates anything; every run writes the one it
measured, ready to commit.

Needs numpy and matplotlib and nothing else: the runner's pyscf, pymbd and
pandas imports are lazy, and importing from it here must not change that.

Usage:
    python pyscf_s66x8_aggregate.py --results-dir all-results [--report-only]
"""

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pyscf_s66x8_figshare_check import N_SCALES, SCHEMA, flatten_rows, summarize

REFERENCE = Path(__file__).resolve().parent / 'pyscf_s66x8_reference.json'
REFERENCE_COMMENT = (
    'Per-system baseline deviation of the PySCF -> pyMBD MBD term from the '
    'FHI-aims values, not zero by construction. Regenerate from a full run, or '
    'null a system out to stop gating it.'
)
# what a result must share with the baseline for the comparison to mean anything
SETTINGS = ('basis', 'xc', 'grid', 'conv_tol')
# gated quantity -> unit, absolute floor of the tolerance, fraction of |baseline|
GATED = (
    ('mbd_rel_dev', '%', 0.5, 0.25),
    ('mbd_mae', ' kcal/mol', 0.010, 0.30),
)
NO_FILE = 'no results (the shard produced no file)'


def load_shards(results_dir):
    """Return (systems by index, shared meta) from every shard file in a directory."""
    systems, source, meta = {}, {}, None
    paths = sorted(Path(results_dir).glob('s66x8-*.json'))
    if not paths:
        raise SystemExit(f'no s66x8-*.json in {results_dir}')
    for path in paths:
        shard = json.loads(path.read_text())
        if shard.get('schema') != SCHEMA:
            raise SystemExit(f'{path}: unsupported schema {shard.get("schema")}')
        if meta is None:
            meta = shard['meta']
        elif [meta[k] for k in SETTINGS] != [shard['meta'][k] for k in SETTINGS]:
            raise SystemExit(f'{path}: settings differ from the other shards')
        for system in shard['systems']:
            idx = system['idx']
            if idx in systems:
                raise SystemExit(
                    f'system {idx} is in both {source[idx]} and {path.name}; '
                    'the shard layout is not a partition'
                )
            systems[idx], source[idx] = system, path.name
    return systems, meta


def system_stats(system):
    """Per-system statistics in the shape the gate and the table want."""
    derived = dict.fromkeys(('mbd_rel_dev', 'mbd_mae', 'mbd_absmax'))
    if system['points']:
        derived.update({k: summarize(system['points'])[k] for k in derived})
    return {
        **{k: system[k] for k in ('idx', 'label', 'error', 'wall')},
        'n': len(system['points']),
        **derived,
    }


def regression_reason(stats, ref):
    """Return None if a system passes the gate, else why it did not.

    ``ref`` of None disables the numeric comparison, which is how a new or
    deliberately unpinned system is carried. A system that errored or ran short
    still fails: that is a broken run rather than an unset baseline.
    """
    if stats['error']:
        return stats['error']
    if stats['n'] != N_SCALES:
        return f'{stats["n"]}/{N_SCALES} points'
    for key, unit, floor, frac in GATED:
        want = (ref or {}).get(key)
        if want is None:
            continue
        got, tol = stats[key], max(floor, frac * abs(want))
        if abs(got - want) > tol:
            return (
                f'{key} {got:+.3f}{unit} vs baseline {want:+.3f}{unit} '
                f'({got - want:+.3f}, tol {tol:.3f})'
            )
    return None


def load_reference(meta):
    """Return (per-system baselines, why the gate is off) for this run's settings."""
    reference = json.loads(REFERENCE.read_text())
    if reference.get('schema') != SCHEMA:
        raise SystemExit(f'{REFERENCE}: unsupported schema {reference.get("schema")}')

    def fmt(settings):
        return ', '.join(f'{k} {settings[k]}' for k in SETTINGS)

    # the settings are the runner's own defaults, so a mismatch means someone
    # changed them: turn the gate off loudly rather than fail on a known change
    if [reference['settings'][k] for k in SETTINGS] != [meta[k] for k in SETTINGS]:
        return (
            {},
            f'this run is {fmt(meta)}, the baseline is {fmt(reference["settings"])}',
        )
    return {int(k): v for k, v in reference['systems'].items()}, None


def _rnd(value, digits):
    return None if value is None else round(value, digits)


def write_reference(path, stats, meta, vdwsets_ref):
    """Write the baseline this run measured, ready to commit as the reference file."""
    head = json.dumps(
        {
            'schema': SCHEMA,
            'comment': REFERENCE_COMMENT,
            'settings': {k: meta[k] for k in SETTINGS},
            'vdwsets_ref': vdwsets_ref,
        },
        indent=1,
    )
    # one line per system: indent=1 throughout would spend five on each, and 66
    # systems of that is most of the file
    systems = ',\n'.join(
        f'  "{s["idx"]}": '
        + json.dumps(
            {
                'label': s['label'],
                'mbd_rel_dev': _rnd(s['mbd_rel_dev'], 3),
                'mbd_mae': _rnd(s['mbd_mae'], 4),
            }
        )
        for s in stats
    )
    # head[:-2] drops its closing "\n}"
    path.write_text(f'{head[:-2]},\n "systems": {{\n{systems}\n }}\n}}\n')


def write_plot(path, systems, meta):
    """The deviation, and what it is worth, against separation.

    Being very nearly a pure offset rather than scatter, the deviation cancels
    against the reference's own error at some separations and compounds at
    others, which the table can state but only this shows.
    """
    import matplotlib

    matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    overall = summarize(flatten_rows(systems.values()))
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
    fig.savefig(path, dpi=150)


def _f(value, spec, unit=''):
    return '-' if value is None else format(value, spec) + unit


def _duration(seconds):
    return f'{seconds / 60:.0f} min' if seconds < 7200 else f'{seconds / 3600:.1f} h'


def format_markdown(stats, rows, meta, reference, gate_off, verdicts):
    """The job summary: headline numbers, the per-system table, the gate's verdict."""
    s = summarize(rows)
    versions = ' - '.join(
        f'{n} {v}' for n, v in meta['versions'].items() if v and n != 'numpy'
    )
    out = [
        '# S66x8: the PySCF -> pyMBD bridge against FHI-aims',
        '',
        f'{meta["basis"]} / {meta["xc"]} / grid {meta["grid"]} / '
        f'conv {meta["conv_tol"]:g} - MBD@rsSCS (a={meta["mbd_a"]}, '
        f'beta={meta["mbd_beta"]}) - {versions}',
        '',
        f'{len(stats)} systems, {s["n"]} points, '
        f'{_duration(sum(x["wall"] for x in stats))} of compute. '
        f'MBD vs FHI-aims: MAE {s["mbd_mae"]:.4f} kcal/mol, '
        f'{s["mbd_rel_dev"]:+.2f} % mean relative. '
        f'MARE vs CCSD(T)/CBS {s["mare_recon"]:.1f} % against the '
        f'reference MBD\'s {s["mare_ref"]:.1f} %.',
        '',
        '## Per system',
        '',
        '| # | System | n | rel dev | baseline | MAE | max abs dev | Status |',
        '| ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |',
    ]
    for x in stats:
        ref = reference.get(x['idx']) or {}
        if x['idx'] in verdicts:
            status = '**error**' if x['error'] else '**regressed**'
        elif gate_off or ref.get('mbd_rel_dev') is None:
            status = 'ungated'
        else:
            status = 'ok'
        out.append(
            f'| {x["idx"]} | {x["label"]} | {x["n"]} | '
            f'{_f(x["mbd_rel_dev"], "+.2f", " %")} | '
            f'{_f(ref.get("mbd_rel_dev"), "+.2f", " %")} | '
            f'{_f(x["mbd_mae"], ".4f")} | {_f(x["mbd_absmax"], ".4f")} | {status} |'
        )
    if gate_off:
        out += ['', f'> The gate is off: {gate_off}.']
    if verdicts:
        out += ['', '## Regressions and failures', '']
        out += [
            f'- **{x["idx"]} {x["label"]}**: {verdicts[x["idx"]]}'
            for x in stats
            if x['idx'] in verdicts
        ]
    out += [
        '',
        'The figure and the baseline this run measured are in the '
        '**s66x8-summary** artifact, the baseline ready to commit as '
        '`devtools/pyscf_s66x8_reference.json`.',
        '',
    ]
    return '\n'.join(out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-dir', required=True, help='directory of shard JSONs')
    p.add_argument('--out-dir', default='results', help='where the report goes')
    p.add_argument('--vdwsets-ref', default='', help='stamped into the baseline')
    p.add_argument('--report-only', action='store_true', help='render without gating')
    args = p.parse_args(argv)

    systems, meta = load_shards(args.results_dir)
    reference, gate_off = load_reference(meta)
    for idx, ref in reference.items():
        # a shard that died outright leaves no file; a placeholder keeps the
        # system in the table and fails the gate instead of dropping it silently
        systems.setdefault(
            idx,
            {
                'idx': idx,
                'label': ref['label'],
                'error': NO_FILE,
                'wall': 0.0,
                'points': [],
            },
        )
    stats = [system_stats(systems[i]) for i in sorted(systems)]
    rows = flatten_rows(systems[i] for i in sorted(systems))
    if not rows:
        raise SystemExit('every system failed; nothing to report')

    verdicts = {
        x['idx']: regression_reason(x, None if gate_off else reference.get(x['idx']))
        for x in stats
    }
    verdicts = {k: v for k, v in verdicts.items() if v}

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_reference(out / 'reference-s66x8.json', stats, meta, args.vdwsets_ref)
    # before the gate's exit code, so a red run still leaves the figure
    write_plot(out / 'convergence-s66x8.png', systems, meta)
    summary = format_markdown(stats, rows, meta, reference, gate_off, verdicts)
    (out / 'summary-s66x8.md').write_text(summary, encoding='utf-8')
    if os.environ.get('GITHUB_STEP_SUMMARY'):
        with open(os.environ['GITHUB_STEP_SUMMARY'], 'a', encoding='utf-8') as f:
            f.write(summary)
    print(summary)
    for idx in sorted(verdicts):
        print(f'{idx}: {verdicts[idx]}', file=sys.stderr)
    return 1 if verdicts and not args.report_only else 0


if __name__ == '__main__':
    raise SystemExit(main())
