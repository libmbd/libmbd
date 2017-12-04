# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import print_function
from ._libmbd import ffi as _ffi, lib as _lib
import numpy as np
from mpi4py import MPI
MPI.COMM_WORLD

bohr = 0.529177249


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


def calculate(coords, alpha_0, omega, R_vdw, beta, a):
    coords = np.array(coords, dtype=float, order='F')/bohr
    alpha_0 = np.array(alpha_0, dtype=float)
    omega = np.array(omega, dtype=float)
    R_vdw = np.array(R_vdw, dtype=float)/bohr
    energy = np.zeros(())
    n_atoms = len(coords)
    calc = _lib.mbd_init_calc(17)
    damping = _lib.mbd_init_damping(
        n_atoms,
        _ffi.cast('double *', R_vdw.ctypes.data),
        beta,
        a)
    _lib.mbd_calculate(
        calc,
        n_atoms,
        _ffi.cast('double *', coords.ctypes.data),
        _ffi.cast('double *', alpha_0.ctypes.data),
        _ffi.cast('double *', omega.ctypes.data),
        damping,
        _ffi.cast('double *', energy.ctypes.data),
    )
    _lib.mbd_destroy_damping(damping)
    _lib.mbd_destroy_calc(calc)
    return float(energy)


# class Settings(object):
#     _fields = {
#         'econv_thr': {},
#         'n_quad_pts': {'dtype': 'intc'},
#         'verbosity': {'dtype': 'intc'},
#         'ewald': {'dtype': 'bool'},
#         'timing': {'dtype': 'bool'},
#         'low_dim': {'dtype': 'bool'},
#         'vacuum': {'shape': (3,), 'dtype': 'bool'}
#     }
#
#     def __init__(self):
#         self._c_sett = ffi.new('struct Settings *')
#         mbdvdw.c_init_settings(self._c_sett)
#         self._sett = {key: get_ndarray(
#             getattr(self._c_sett, key),
#             **self._fields[key]
#         ) for key in self.keys()}
#
#     def keys(self):
#         return dir(self._c_sett)
#
#     def __getitem__(self, key):
#         return self._sett[key]
#
#     def __setitem__(self, key, value):
#         if 'shape' in self._fields[key]:
#             self._sett[key][:] = value
#         else:
#             self._sett[key][()] = value
#
#     def __repr__(self):
#         return pformat(self._sett)
#
#
# settings = Settings()
#
#
# if __name__ == '__main__':
#     print(calculate(
#         [[0, 0, 0], [4/bohr, 0, 0]],
#         ['Ar', 'Ar'],
#         [1, 1],
#         0.85,
#         get_forces=True
#     ))
#     print(calculate(
#         [[0, 0, 0]],
#         ['Ar'],
#         [1],
#         0.85,
#         lattice=4/bohr*np.eye(3),
#         k_grid=[4, 4, 4],
#         get_forces=True
#     ))
