#!/usr/bin/env python3
"""Summarise local coverage the way Codecov does.

Reads the .gcov files produced by gen-fortran-coverage.sh and the Python
coverage.xml, and reports per-file / total line coverage using Codecov's
semantics: coverage is tracked per *source line*, and a line counts as hit if
*any* of its gcov entries executed. gcov emits some source lines more than once
(e.g. `end module`, or template-instantiated lines), so a naive count of raw
.gcov rows under-reports coverage; de-duplicating by line number is what makes
the local number match Codecov.

Usage: coverage-summary.py [--gcov-dir coverage-fortran] [--xml coverage.xml]
"""

import argparse
import glob
import os
import xml.etree.ElementTree as ET


def gcov_file_coverage(path):
    """Return (source_path, {lineno: hit_bool}) merged per source line."""
    lines = {}
    src = None
    with open(path, encoding='utf-8', errors='replace') as f:
        for row in f:
            # Format: '<count>:<lineno>:<source>'; headers use lineno 0.
            parts = row.split(':', 2)
            if len(parts) < 3:
                continue
            count, lineno = parts[0].strip(), parts[1].strip()
            if lineno == '0':
                if count == '-' and parts[2].startswith('Source:'):
                    src = parts[2][len('Source:') :].strip()
                continue
            if count == '-':  # non-executable line
                continue
            ln = int(lineno)
            hit = count not in ('#####', '=====')  # any digit count => executed
            lines[ln] = lines.get(ln, False) or hit  # covered-if-any
    return src, lines


def repo_relative(src):
    """Map an absolute source path to its repo-relative src/... form."""
    if '/src/' in src:
        return 'src/' + src.split('/src/', 1)[1]
    return os.path.basename(src)


def python_coverage(xml_path):
    """Return {filename: (hits, misses)} from a coverage.xml report."""
    out = {}
    root = ET.parse(xml_path).getroot()
    for cls in root.iter('class'):
        fn = cls.get('filename')
        hits = total = 0
        for line in cls.iter('line'):
            total += 1
            hits += 1 if int(line.get('hits')) > 0 else 0
        out[fn] = (hits, total - hits)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gcov-dir', default='coverage-fortran')
    p.add_argument('--xml', default='coverage.xml')
    args = p.parse_args()

    files = {}  # name -> (hits, misses)
    for g in glob.glob(os.path.join(args.gcov_dir, '*.gcov')):
        src, lines = gcov_file_coverage(g)
        name = repo_relative(src) if src else os.path.basename(g)
        hits = sum(1 for v in lines.values() if v)
        misses = sum(1 for v in lines.values() if not v)
        files[name] = (hits, misses)  # one .gcov per source with gcov -p

    if os.path.exists(args.xml):
        files.update(python_coverage(args.xml))

    total_hits = total_lines = 0
    print(f'{"file":32} {"cov%":>7} {"hit":>5} {"lines":>6}')
    for name in sorted(files):
        hits, misses = files[name]
        tot = hits + misses
        total_hits += hits
        total_lines += tot
        cov = 100 * hits / tot if tot else 0.0
        print(f'{name.replace("src/", ""):32} {cov:7.2f} {hits:5d} {tot:6d}')
    cov = 100 * total_hits / total_lines if total_lines else 0.0
    print(f'{"TOTAL":32} {cov:7.2f} {total_hits:5d} {total_lines:6d}')


if __name__ == '__main__':
    main()
