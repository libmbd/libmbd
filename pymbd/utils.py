# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import print_function, division

import numpy as np


def numerical_gradients(f, coords, *args, **kwargs):
    delta = kwargs.pop('delta', 1e-3)  # support python 2
    coords = np.array(coords)
    gradients = np.zeros(coords.shape)
    for i_atom in range(coords.shape[0]):
        for i_xyz in range(3):
            ene = {}
            for step in [-2, -1, 1, 2]:
                coords_diff = coords.copy()
                coords_diff[i_atom, i_xyz] += step*delta
                ene[step] = f(coords_diff, *args, **kwargs)
            gradients[i_atom, i_xyz] = _diff5(ene, delta)
    return gradients


def _diff5(x, delta):
    return (1./12*x[-2]-2./3*x[-1]+2./3*x[1]-1./12*x[2])/delta
