# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import division, print_function

import numpy as np

from .pymbd import _array, from_volumes
from ._libmbd import ffi as _ffi
from ._libmbd import lib as _lib

with_scalapack = _lib.with_scalapack
if with_scalapack:
    from mpi4py import MPI  # noqa


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

    def ts_energy(self, coords, alpha_0, C6, R_vdw, sR,
                  lattice=None, d=20., damping='fermi'):
        coords = _array(coords, dtype=float)
        alpha_0 = _array(alpha_0, dtype=float)
        C6 = _array(C6, dtype=float)
        R_vdw = _array(R_vdw, dtype=float)
        lattice = _array(lattice, dtype=float)
        n_atoms = len(coords)
        system = _lib.mbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _cast('double*', lattice),
            _ffi.NULL,
        )
        damping = _lib.mbd_init_damping(
            n_atoms, damping.encode(), _cast('double*', R_vdw), _ffi.NULL, sR, d,
        )
        ene = _lib.calc_ts_energy(
            system,
            n_atoms,
            _cast('double*', alpha_0),
            _cast('double*', C6),
            damping,
            _ffi.NULL,
        )
        _lib.mbd_destroy_damping(damping)
        _lib.mbd_destroy_system(system)
        return ene

    def mbd_energy(self, coords, alpha_0, C6, R_vdw, beta,
                   lattice=None, k_grid=None,
                   a=6., func='calc_mbd_rsscs_energy', force=False,
                   damping='fermi,dip'):
        coords = _array(coords, dtype=float)
        alpha_0 = _array(alpha_0, dtype=float)
        C6 = _array(C6, dtype=float)
        R_vdw = _array(R_vdw, dtype=float)
        lattice = _array(lattice, dtype=float)
        k_grid = _array(k_grid, dtype='i4')
        n_atoms = len(coords)
        system = _lib.mbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _cast('double*', lattice),
            _cast('int*', k_grid),
        )
        damping = _lib.mbd_init_damping(
            n_atoms, damping.encode(), _cast('double*', R_vdw), _ffi.NULL, beta, a,
        )
        gradients = np.zeros((n_atoms, 3)) if force else None
        ene = getattr(_lib, func)(
            system,
            n_atoms,
            _cast('double*', alpha_0),
            _cast('double*', C6),
            damping,
            _cast('double*', gradients),
        )
        _lib.mbd_destroy_damping(damping)
        _lib.mbd_destroy_system(system)
        if force:
            return ene, gradients
        return ene

    def dipole_matrix(self, coords, damping, beta=0., lattice=None, k_point=None,
                      R_vdw=None, sigma=None, a=6.):
        coords = _array(coords, dtype=float)
        R_vdw = _array(R_vdw, dtype=float)
        sigma = _array(sigma, dtype=float)
        lattice = _array(lattice, dtype=float)
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

    def mbd_energy_species(self, coords, species, vols, beta, **kwargs):
        alpha_0, C6, R_vdw = from_volumes(species, vols)
        return self.mbd_energy(coords, alpha_0, C6, R_vdw, beta, **kwargs)

    def ts_energy_species(self, coords, species, vols, beta, **kwargs):
        alpha_0, C6, R_vdw = from_volumes(species, vols)
        return self.ts_energy(coords, alpha_0, C6, R_vdw, beta, **kwargs)


def _ndarray(ptr, shape=None, dtype='float'):
    return np.ndarray(
        buffer=_ffi.buffer(
            ptr,
            (np.prod(shape) if shape else 1)*np.dtype(dtype).itemsize
        ),
        shape=shape,
        dtype=dtype,
    )


def _cast(ctype, array):
    return _ffi.NULL if array is None else _ffi.cast(ctype, array.ctypes.data)
