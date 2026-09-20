# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import division, print_function

import numpy as np

__all__ = ()


def numerical_gradients(geom, func, *args, delta=1e-3, **kwargs):
    steps, diff = finite_diff_gen(kwargs.pop('npts', 5))
    coords_0 = geom.coords
    gradients = np.zeros(coords_0.shape)
    for i_atom in range(coords_0.shape[0]):
        for i_xyz in range(3):
            ene = {}
            for step in steps:
                coords = coords_0.copy()
                coords[i_atom, i_xyz] += step * delta
                geom.coords = coords
                ene[step] = getattr(geom, func)(*args, **kwargs)
            gradients[i_atom, i_xyz] = diff(ene, delta)
    return gradients


def numerical_latt_gradients(geom, func, *args, **kwargs):
    delta = kwargs.pop('delta', 1e-3)
    steps, diff = finite_diff_gen(kwargs.pop('npts', 5))
    lattice_0 = geom.lattice
    gradients = np.zeros((3, 3))
    for i_vec in range(3):
        for i_xyz in range(3):
            ene = {}
            for step in steps:
                lattice = lattice_0.copy()
                lattice[i_vec, i_xyz] += step * delta
                geom.lattice = lattice
                ene[step] = getattr(geom, func)(*args, **kwargs)
            gradients[i_vec, i_xyz] = diff(ene, delta)
    return gradients


def numerical_vdw_params_gradients(func, alpha_0, C6, R_vdw, delta=1e-3, npts=5):
    """Numerically differentiate ``func(alpha_0, C6, R_vdw)`` in its arguments.

    Unlike the coordinate steps above, the step is relative, because the three
    kinds of vdW parameter differ in magnitude by an order of magnitude or two.
    Returns the three gradients in the order of the arguments.
    """
    steps, diff = finite_diff_gen(npts)
    params_0 = [np.array(param, dtype=float) for param in (alpha_0, C6, R_vdw)]
    gradients = []
    for i_param, param_0 in enumerate(params_0):
        gradient = np.zeros_like(param_0)
        for i_atom in range(len(param_0)):
            step_size = delta * param_0[i_atom]
            ene = {}
            for step in steps:
                params = [param.copy() for param in params_0]
                params[i_param][i_atom] += step * step_size
                ene[step] = func(*params)
            gradient[i_atom] = diff(ene, step_size)
        gradients.append(gradient)
    return gradients


def _diff3(x, delta):
    return (-1.0 / 2 * x[-1] + 1.0 / 2 * x[1]) / delta


def _diff5(x, delta):
    return (
        1.0 / 12 * x[-2] - 2.0 / 3 * x[-1] + 2.0 / 3 * x[1] - 1.0 / 12 * x[2]
    ) / delta


def finite_diff_gen(npts):
    return {3: ([-1, 1], _diff3), 5: ([-2, -1, 1, 2], _diff5)}[npts]
