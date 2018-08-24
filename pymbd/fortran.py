# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import division, print_function

import numpy as np

from .pymbd import _array, from_volumes
from ._libmbd import ffi as _ffi, lib as _lib

with_mpi = _lib.cmbd_with_mpi
with_scalapack = _lib.cmbd_with_scalapack

if with_mpi:
    from mpi4py import MPI  # noqa


class MBDFortranException(Exception):
    def __init__(self, code, origin, msg):
        super(MBDFortranException, self).__init__(msg)
        self.code = code
        self.origin = origin


class MBDCalc(object):
    def __init__(self, n_freq=15):
        self._calc_obj = None
        self.n_freq = n_freq

    def __enter__(self):
        self._calc_obj = _lib.cmbd_init_calc(self.n_freq)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        _lib.cmbd_destroy_calc(self._calc)
        self._calc_obj = None

    def _check_exc(self):
        code = _array(0, dtype=int)
        origin = _ffi.new('char[50]')
        msg = _ffi.new('char[150]')
        _lib.cmbd_get_exception(self._calc, _cast('int*', code), origin, msg)
        if code != 0:
            raise MBDFortranException(
                int(code), _ffi.string(origin).decode(), _ffi.string(msg).decode()
            )

    @property
    def _calc(self):
        if not self._calc_obj:
            raise RuntimeError('MBDCalc must be used as a context manager')
        return self._calc_obj

    @property
    def omega_grid(self):
        return (
            _ndarray(self._calc.omega_grid, (self.n_freq+1,)).copy(),
            _ndarray(self._calc.omega_grid_w, (self.n_freq+1,)).copy(),
        )

    def ts_energy(self, coords, alpha_0, C6, R_vdw, sR,
                  lattice=None, d=20., damping='fermi'):
        coords, alpha_0, C6, R_vdw, lattice = \
            map(_array, (coords, alpha_0, C6, R_vdw, lattice))
        n_atoms = len(coords)
        system = _lib.cmbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _cast('double*', lattice),
            _ffi.NULL,
        )
        damping = _lib.cmbd_init_damping(
            n_atoms, damping.encode(), _cast('double*', R_vdw), _ffi.NULL, sR, d
        )
        ene = _lib.cmbd_ts_energy(
            system,
            n_atoms,
            _cast('double*', alpha_0),
            _cast('double*', C6),
            damping,
            _ffi.NULL,
        )
        _lib.cmbd_destroy_damping(damping)
        _lib.cmbd_destroy_system(system)
        self._check_exc()
        return ene

    def mbd_energy(self, coords, alpha_0, C6, R_vdw, beta,
                   lattice=None, k_grid=None,
                   a=6., func='mbd_rsscs_energy', force=False,
                   damping='fermi,dip', spectrum=False):
        coords, alpha_0, C6, R_vdw, lattice = \
            map(_array, (coords, alpha_0, C6, R_vdw, lattice))
        k_grid = _array(k_grid, dtype='i4')
        n_atoms = len(coords)
        system = _lib.cmbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _cast('double*', lattice),
            _cast('int*', k_grid),
        )
        damping = _lib.cmbd_init_damping(
            n_atoms, damping.encode(), _cast('double*', R_vdw), _ffi.NULL, beta, a,
        )
        gradients = np.zeros((n_atoms, 3)) if force else None
        eigs, modes = None, None
        args = (
            system,
            n_atoms,
            _cast('double*', alpha_0),
            _cast('double*', C6),
            damping,
            _cast('double*', gradients),
        )
        if func == 'mbd_rsscs_energy':
            if spectrum:
                eigs = np.zeros(3*n_atoms)
                modes = np.zeros((3*n_atoms, 3*n_atoms), order='F')
            args += (_cast('double*', eigs), _cast('double*', modes))
        ene = getattr(_lib, 'cmbd_' + func)(*args)
        _lib.cmbd_destroy_damping(damping)
        _lib.cmbd_destroy_system(system)
        self._check_exc()
        if spectrum:
            ene = ene, eigs, modes
        if force:
            return ene, gradients
        return ene

    def dipole_matrix(self, coords, damping, beta=0., lattice=None, k_point=None,
                      R_vdw=None, sigma=None, a=6.):
        coords, R_vdw, sigma, lattice, k_point = \
            map(_array, (coords, R_vdw, sigma, lattice, k_point))
        n_atoms = len(coords)
        system = _lib.cmbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _cast('double*', lattice),
            _ffi.NULL,
        )
        damping = _lib.cmbd_init_damping(
            n_atoms, damping.encode(),
            _cast('double*', R_vdw),
            _cast('double*', sigma),
            beta, a,
        )
        dipmat = np.empty(
            (3*n_atoms, 3*n_atoms),
            dtype=float if k_point is None else complex,
        )
        _lib.cmbd_dipole_matrix(
            system,
            damping,
            _cast('double*', k_point),
            _cast('double*', dipmat),
        )
        _lib.cmbd_destroy_damping(damping)
        _lib.cmbd_destroy_system(system)
        self._check_exc()
        return dipmat

    def mbd_energy_species(self, coords, species, vols, beta, **kwargs):
        alpha_0, C6, R_vdw = from_volumes(species, vols)
        return self.mbd_energy(coords, alpha_0, C6, R_vdw, beta, **kwargs)

    def ts_energy_species(self, coords, species, vols, beta, **kwargs):
        alpha_0, C6, R_vdw = from_volumes(species, vols)
        return self.ts_energy(coords, alpha_0, C6, R_vdw, beta, **kwargs)

    def dipole_energy(self, coords, a0, w, w_t, version, r_vdw, beta, a, C):
        n_atoms = len(coords)
        system = _lib.cmbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _ffi.NULL,
            _ffi.NULL,
        )
        res = _lib.cmbd_dipole_energy(
            system,
            n_atoms,
            _cast('double*', a0),
            _cast('double*', w),
            _cast('double*', w_t),
            version.encode(),
            _cast('double*', r_vdw),
            beta,
            a,
            _cast('double*', C),
        )
        self._check_exc()
        return res

    def coulomb_energy(self, coords, q, m, w_t, version, r_vdw, beta, a, C):
        n_atoms = len(coords)
        system = _lib.cmbd_init_system(
            self._calc,
            n_atoms,
            _cast('double*', coords),
            _ffi.NULL,
            _ffi.NULL,
        )
        res = _lib.cmbd_coulomb_energy(
            system,
            n_atoms,
            _cast('double*', q),
            _cast('double*', m),
            _cast('double*', w_t),
            version.encode(),
            _cast('double*', r_vdw),
            beta,
            a,
            _cast('double*', C),
        )
        self._check_exc()
        return res


def _ndarray(ptr, shape=None, dtype='float'):
    buffer_size = (np.prod(shape) if shape else 1)*np.dtype(dtype).itemsize
    return np.ndarray(
        buffer=_ffi.buffer(ptr, buffer_size), shape=shape, dtype=dtype,
    )


def _cast(ctype, array):
    return _ffi.NULL if array is None else _ffi.cast(ctype, array.ctypes.data)
