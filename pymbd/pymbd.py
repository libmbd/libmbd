# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import print_function
import numpy as np
import pkg_resources
import sys
import csv
from itertools import product
try:
    from mpi4py import MPI
    MPI.COMM_WORLD
except ImportError:
    pass

from ._libmbd import ffi as _ffi, lib as _lib

ang = 1/0.529177249


def _array(obj, *args, **kwargs):
    if obj is not None:
        return np.array(obj, *args, **kwargs)


def _cast(ctype, array):
    return _ffi.NULL if array is None else _ffi.cast(ctype, array.ctypes.data)


class MBDCalc(object):
    def __init__(self, n_freq=15):
        self.__calc = None
        self.n_freq = n_freq

    def __enter__(self):
        self.__calc = _lib.mbd_init_calc(self.n_freq)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        _lib.mbd_destroy_calc(self._calc)
        self.__calc = None

    @property
    def _calc(self):
        if not self.__calc:
            raise RuntimeError('MBDCalc must be used as a context manager')
        return self.__calc

    @property
    def omega_grid(self):
        return (
            _ndarray(self._calc.omega_grid, (self._calc.n_freq+1,)).copy(),
            _ndarray(self._calc.omega_grid_w, (self._calc.n_freq+1,)).copy(),
        )

    def mbd_energy(self, coords, alpha_0, C6, R_vdw, beta,
                   lattice=None, k_grid=None,
                   a=6., func='calc_mbd_rsscs_energy', force=False,
                   damping='fermi,dip'):
        coords = _array(coords, dtype=float, order='F')
        alpha_0 = _array(alpha_0, dtype=float)
        C6 = _array(C6, dtype=float)
        R_vdw = _array(R_vdw, dtype=float)
        lattice = _array(lattice, dtype=float, order='F')
        k_grid = _array(k_grid, dtype='i4')
        n_atoms = len(coords)
        system = _lib.mbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _cast('double*', lattice),
            _cast('int*', k_grid),
        )
        if force:
            system.do_force[0] = True
        damping = _lib.mbd_init_damping(
            n_atoms, damping.encode(), _cast('double*', R_vdw), _ffi.NULL, beta, a,
        )
        ene = getattr(_lib, func)(
            system,
            n_atoms,
            _cast('double*', alpha_0),
            _cast('double*', C6),
            damping
        )
        if force:
            ret_val = ene, _ndarray(system.forces, (n_atoms, 3)).copy()
        else:
            ret_val = ene
        _lib.mbd_destroy_damping(damping)
        _lib.mbd_destroy_system(system)
        return ret_val

    def dipole_matrix(self, coords, damping, beta=0., lattice=None, k_point=None,
                      R_vdw=None, sigma=None, a=6.):
        coords = _array(coords, dtype=float, order='F')
        R_vdw = _array(R_vdw, dtype=float)
        sigma = _array(sigma, dtype=float)
        lattice = _array(lattice, dtype=float, order='F')
        k_point = _array(k_point, dtype=float)
        n_atoms = len(coords)
        system = _lib.mbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _cast('double*', lattice),
            _ffi.NULL,
        )
        damping = _lib.mbd_init_damping(
            n_atoms, damping.encode(),
            _cast('double*', R_vdw),
            _cast('double*', sigma),
            beta, a,
        )
        dipmat = np.empty(
            (3*n_atoms, 3*n_atoms),
            dtype=float if k_point is None else complex,
            order='F'
        )
        _lib.calc_dipole_matrix(
            system,
            damping,
            _cast('double*', k_point),
            _cast('double*', dipmat),
        )
        _lib.mbd_destroy_damping(damping)
        _lib.mbd_destroy_system(system)
        return dipmat

    def pymbd_energy(self, coords, alpha_0, C6, R_vdw, beta,
                     lattice=None, k_grid=None):
        coords = _array(coords, dtype=float, order='F')
        alpha_0 = _array(alpha_0, dtype=float)
        C6 = _array(C6, dtype=float)
        R_vdw = _array(R_vdw, dtype=float)
        freq, freq_w = self.omega_grid
        omega = 4./3*C6/alpha_0**2
        alpha_dyn = alpha_0/(1+(freq[:, None]/omega)**2)
        alpha_dyn_rsscs = np.empty_like(alpha_dyn)
        for a, a_scr in zip(alpha_dyn, alpha_dyn_rsscs):
            sigma = (np.sqrt(2./np.pi)*a/3)**(1./3)
            a_nlc = np.linalg.inv(
                np.diag(np.repeat(1./a, 3)) +
                self.dipole_matrix(
                    coords, 'fermi,dip,gg', sigma=sigma, R_vdw=R_vdw, beta=beta,
                    lattice=lattice
                )
            )
            a_scr[:] = np.array(
                [a_nlc[i::3, i::3].sum(1) for i in range(3)]
            ).sum(0)/3
        C6_rsscs = 3./np.pi*np.sum(freq_w[:, None]*alpha_dyn_rsscs**2, 0)
        R_vdw_rsscs = R_vdw*(alpha_dyn_rsscs[0, :]/alpha_0)**(1./3)
        omega_rsscs = 4./3*C6_rsscs/alpha_dyn_rsscs[0, :]**2
        pre = np.repeat(omega_rsscs*np.sqrt(alpha_dyn_rsscs[0, :]), 3)
        if lattice is None:
            k_grid = [None]
        else:
            assert k_grid is not None
            k_grid = get_kgrid(lattice, k_grid)
        ene = 0
        for k_point in k_grid:
            eigs = np.linalg.eigvalsh(
                np.diag(np.repeat(omega_rsscs**2, 3)) +
                np.outer(pre, pre)*self.dipole_matrix(
                    coords, 'fermi,dip', R_vdw=R_vdw_rsscs, beta=beta,
                    lattice=lattice, k_point=k_point
                )
            )
            ene += np.sum(np.sqrt(eigs))/2-3*np.sum(omega_rsscs)/2
        ene /= len(k_grid)
        return ene

    def mbd_energy_species(self, coords, species, volumes, beta, **kwargs):
        alpha_0, C6, R_vdw = from_volumes(species, volumes)
        return self.mbd_energy(coords, alpha_0, C6, R_vdw, beta, **kwargs)

    def pymbd_energy_species(self, coords, species, volumes, beta, **kwargs):
        alpha_0, C6, R_vdw = from_volumes(species, volumes)
        return self.pymbd_energy(coords, alpha_0, C6, R_vdw, beta, **kwargs)


def from_volumes(species, volumes):
    alpha_0, C6, R_vdw = (
        np.array([vdw_params[sp][param] for sp in species])
        for param in ['alpha_0', 'C6', 'R_vdw']
    )
    volumes = np.array(volumes)
    alpha_0 *= volumes
    C6 *= volumes**2
    R_vdw *= volumes**(1./3)
    return alpha_0, C6, R_vdw


def get_kgrid(lattice, k_grid, shift=0.5):
    k_grid = np.array(k_grid)
    lattice = np.array(lattice)
    idx_grid = np.array(list(product(*map(range, k_grid))))+shift
    idx_grid /= k_grid
    idx_grid = np.where(idx_grid > 0.5, idx_grid-1, idx_grid)
    recp_lattice = 2*np.pi*np.linalg.inv(lattice.T)
    k_grid = idx_grid.dot(recp_lattice)
    return k_grid


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
