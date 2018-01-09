# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import print_function
import numpy as np
import pkg_resources
import sys
import csv
try:
    from mpi4py import MPI
    MPI.COMM_WORLD
except ImportError:
    pass

from ._libmbd import ffi as _ffi, lib as _lib

ang = 1/0.529177249


class MBDCalc(object):
    def __init__(self):
        self._calc = None

    def __enter__(self):
        self._calc = _lib.mbd_init_calc()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        _lib.mbd_destroy_calc(self._calc)
        self._calc = None

    def mbd_energy(self, coords, alpha_0, C6, R_vdw, beta,
                   lattice=None, k_grid=None,
                   a=6., func='calc_mbd_rsscs_energy', force=False,
                   damping='fermi,dip'):
        if not self._calc:
            raise RuntimeError('MBDCalc must be used as a context manager')
        coords = np.array(coords, dtype=float, order='F')
        alpha_0 = np.array(alpha_0, dtype=float)
        C6 = np.array(C6, dtype=float)
        R_vdw = np.array(R_vdw, dtype=float)
        periodic = lattice is not None
        if periodic:
            lattice = np.array(lattice, dtype=float, order='F')
            k_grid = np.array(k_grid, dtype='i4')
        n_atoms = len(coords)
        system = _lib.mbd_init_system(
            self._calc,
            n_atoms,
            _ffi.cast('double*', coords.ctypes.data),
            periodic,
            _ffi.cast('double*', lattice.ctypes.data) if periodic else _ffi.NULL,
            _ffi.cast('int*', k_grid.ctypes.data) if periodic else _ffi.NULL,
        )
        if force:
            system.do_force[0] = True
        damping = _lib.mbd_init_damping(
            n_atoms, damping.encode(), _ffi.cast('double*', R_vdw.ctypes.data), beta, a,
        )
        ene = getattr(_lib, func)(
            system,
            n_atoms,
            _ffi.cast('double*', alpha_0.ctypes.data),
            _ffi.cast('double*', C6.ctypes.data),
            damping
        )
        if force:
            ret_val = ene, _ndarray(system.forces, (n_atoms, 3)).copy()
        else:
            ret_val = ene
        _lib.mbd_destroy_damping(damping)
        _lib.mbd_destroy_system(system)
        return ret_val

    def mbd_energy_species(self, coords, species, volumes, beta, **kwargs):
        alpha_0, C6, R_vdw = (
            np.array([vdw_params[sp][param] for sp in species])
            for param in ['alpha_0', 'C6', 'R_vdw']
        )
        volumes = np.array(volumes)
        alpha_0 *= volumes
        C6 *= volumes**2
        R_vdw *= volumes**(1./3)
        return self.mbd_energy(coords, alpha_0, C6, R_vdw, beta, **kwargs)


def numerical_forces(f, coords, *args, **kwargs):
    delta = kwargs.pop('delta', 1e-3)  # support python 2
    coords = np.array(coords)
    forces = np.zeros(coords.shape)
    for i_atom in range(coords.shape[0]):
        for i_xyz in range(3):
            ene = {}
            for step in [-2, -1, 1, 2]:
                coords_diff = coords.copy()
                coords_diff[i_atom, i_xyz] += step*delta
                ene[step] = f(coords_diff, *args, **kwargs)
            forces[i_atom, i_xyz] = _diff5(ene, delta)
    return forces


def _diff5(x, delta):
    return (1./12*x[-2]-2./3*x[-1]+2./3*x[1]-1./12*x[2])/delta


def _get_vdw_params():
    csv_lines = pkg_resources.resource_string(__name__, 'vdw-params.csv').split(b'\n')
    if sys.version_info[0] > 2:
        csv_lines = [l.decode() for l in csv_lines]
    reader = csv.DictReader(csv_lines, quoting=csv.QUOTE_NONNUMERIC)
    vdw_params = {}
    for row in reader:
        vdw_params[row.pop('species')] = row
    return vdw_params


vdw_params = _get_vdw_params()


def _ndarray(ptr, shape=None, dtype='float'):
    return np.ndarray(
        buffer=_ffi.buffer(
            ptr,
            (np.prod(shape) if shape else 1)*np.dtype(dtype).itemsize
        ),
        shape=shape,
        dtype=dtype,
        order='F'
    )
