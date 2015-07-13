from __future__ import print_function
import pymbd
import json
import numpy as np

energy = pymbd.main('mbddata.json')
data = json.load(open('mbddata.json'))
errors = []
for key in energy:
    my = energy[key]
    ref = [e['value'] for e in data['energy'] if e['name'] == key][0]
    if np.linalg.norm(my-ref) > 1e-10:
        errors.append((key, ref, my))
if not errors:
    print('Success: all energies match')
else:
    print('There were some erors:')
    for key, ref, my in errors:
        print(key, ref, my)
