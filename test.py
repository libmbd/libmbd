import pymbd
import json
import numpy as np

data = json.load(open('results.json'))
energy = pymbd.main('results.json')
for key in energy:
    my = energy[key]
    ref = [e['value'] for e in data['energy'] if e['name'] == key][0]
    assert np.linalg.norm(my-ref) < 1e-10
print('Success: all energies match')
