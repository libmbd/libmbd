from __future__ import print_function
import json
import numpy as np
import sys

results = json.load(sys.stdin)
with open('argon-dimer.json') as f:
    data = json.load(f)
errors = []
for key in results:
    my = results[key]
    for ene in data['energy']:
        if ene['name'] == key:
            ref = ene['value']
            break
    if np.linalg.norm(np.array(my)-np.array(ref)) > 1e-10:
        errors.append((key, ref, my))
if not errors:
    print('Success: all energies match')
else:
    print('There were some erors:')
    for key, ref, my in errors:
        print(key, ref, my)
