#!/usr/bin/env python3
import json
import numpy as np
import sys


colors = {
    'red': '\x1b[01;31m',
    'green': '\x1b[32m',
    'yellow': '\x1b[33m',
    'normal': '\x1b[0m',
}


def colstr(s, color):
    return colors[color] + s + colors['normal']


results = json.load(sys.stdin)
with open(sys.argv[1]) as f:
    refs = json.load(f)
errors = []
for key in results:
    my = results[key]
    ref = refs[key]
    if np.linalg.norm(np.array(my)-np.array(ref)) > 1e-10:
        errors.append((key, ref, my))
if not errors:
    print(colstr('OK: all energies match', 'green'))
else:
    print(colstr('FAIL: there were some erors:', 'red'))
    for key, ref, my in errors:
        print(key, ref, my)
