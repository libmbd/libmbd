#!/usr/bin/env python3
import re
from pathlib import Path
from string import Template
import subprocess


def parse_modules(lines):
    defined = []
    used = set()
    for line in lines:
        line = line.lstrip()
        if not line:
            continue
        if line[0] == '!':
            continue
        word = line.split(' ', 1)[0].lower()
        if word == 'module':
            module = re.match(
                r'module\s+(\w+)\s*', line, re.IGNORECASE
            ).group(1).lower()
            if module != 'procedure':
                defined.append(module)
        elif word == 'use':
            module = re.match(
                r'use\s+(\w+)\s*', line, re.IGNORECASE
            ).group(1).lower()
            used.add(module)
    used.difference_update(defined)
    return defined, used


def topsort(tree):
    idxs = {node: n for n, node in enumerate(tree)}
    outgoing = [
        [idxs[child] for child in children]
        for node, children in tree.items()
    ]
    N = len(tree)
    nincoming = N*[0]
    for edges in outgoing:
        for n in edges:
            nincoming[n] += 1
    L = []
    S = [n for n in range(N) if nincoming[n] == 0]
    while S:
        n = S.pop()
        L.append(n)
        for m in outgoing[n]:
            nincoming[m] -= 1
            if not nincoming[m]:
                S.append(m)
    if sum(nincoming):
        raise RuntimeError('Graph with cycles')
    iidxs = {n: node for node, n in idxs.items()}
    return [iidxs[n] for n in L]


def dep_tree(paths, macros=None):
    defs = {}
    used = {}
    for path in paths:
        if path.suffix == '.F90':
            args = ['gfortran', '-E', '-P', str(path)]
            if macros:
                args += [f'-D{macro}' for macro in macros]
            lines = subprocess.run(args, capture_output=True) \
                .stdout.decode().split('\n')
        else:
            lines = path.read_text().split('\n')
        modules, used[path] = parse_modules(lines)
        for mod in modules:
            defs[mod] = path
    deps = {name: sorted(set(
        defs[mod] for mod in mods if mod not in {'iso_c_binding', 'mpi'}
    )) for name, mods in used.items()}
    ordered = list(reversed(topsort(deps)))
    return ordered, deps


if __name__ == '__main__':
    root = Path(__file__).parents[1]
    sources = sorted(
        path for path in (root/'src').glob('*.?90')
        if 'tests' not in path.name and path.name != 'mbd_density.f90'
    )
    configs = [
        (
            {'mbd_mpi.F90', 'mbd_blacs.f90', 'mbd_scalapack.f90', 'mbd_elsi.F90'},
            [],
            'serial'
        ),
        (
            {'mbd_elsi.F90'},
            ['WITH_SCALAPACK', 'WITH_MPI'],
            'scalapack-mpi'
        ),
    ]
    config_file = root/'setup.cfg'
    for filtered, macros, label in configs:
        srcs = [src for src in sources if src.name not in filtered]
        ordered, deps = dep_tree(srcs, macros)
        template = Template((root/'utils/Makefile.in').read_text())
        objs = {path: re.sub(r'\.[fF]90$', '.o', path.name) for path in srcs}
        (root/'src'/f'{label}.make').write_text(template.substitute(
            objs=' '.join(objs[path] for path in srcs),
            deps='\n'.join(
                f'{objs[path]}: {" ".join(objs[p] for p in sorted(deps[path]))}'
                for path in srcs
            ),
            macros=' '.join(f'-D{macro}' for macro in macros),
        ))
