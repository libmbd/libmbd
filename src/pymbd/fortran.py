# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import division, print_function

from functools import wraps

import numpy as np

from .__init__ import __version__ as PYMBD_VERSION
from .pymbd import _array, from_volumes

try:
    from ._libmbd import ffi as _ffi, lib as _lib
except ImportError:
    raise Exception('Pymbd C extension unimportable, cannot use Fortran')

__all__ = ['MBDGeom', 'with_mpi', 'with_scalapack']

with_mpi = _lib.cmbd_with_mpi
"""Whether Libmbd was compiled with MPI"""

with_scalapack = _lib.cmbd_with_scalapack
"""Whether Libmbd was compiled with Scalapack"""

LIBMBD_VERSION = (
    _lib.cmbd_version_major,
    _lib.cmbd_version_minor,
    _lib.cmbd_version_patch,
)

# do not test versions when running autodoc
if PYMBD_VERSION and isinstance(LIBMBD_VERSION[0], int):
    assert PYMBD_VERSION[0] == LIBMBD_VERSION[0]
    if PYMBD_VERSION[0] == 0:
        assert PYMBD_VERSION[1] == LIBMBD_VERSION[1]
    else:
        assert PYMBD_VERSION[1] <= LIBMBD_VERSION[1]

if with_mpi:
    from mpi4py import MPI  # noqa


class MBDFortranException(Exception):
    def __init__(self, code, origin, msg):
        super(MBDFortranException, self).__init__(msg)
        self.code = code
        self.origin = origin


def _auto_context(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        if self._geom_f:
            return method(self, *args, **kwargs)
        with self:
            return method(self, *args, **kwargs)

    return wrapper


class MBDGeom(object):
    """Represents an initialized Libmbd `geom_t <../type/geom_t.html>`_ object.

    :param array-like coords: (a.u.) atomic coordinates as rows
    :param array-like lattice: (a.u.) lattice vectors as rows
    :param array-like k_grid: number of :math:`k`-points per reciprocal lattice vector
    :param array-like custom_k_pts: (a.u.) custom :math:`k`-points as rows
    :param int n_freq: number of quadrature points for frequency integration
    :param bool do_rpa: whether to calculate MBD energy via frequency integration
    :param bool get_spectrum: whether to return eigenvalues and eigenvectors
    :param bool get_rpa_orders: whether to return RPA order decomposition
    :param bool rpa_rescale_eigs: whether to rescale RPA eigenvalues
    """

    def __init__(
        self,
        coords,
        lattice=None,
        k_grid=None,
        custom_k_pts=None,
        n_freq=None,
        do_rpa=False,
        get_spectrum=False,
        get_rpa_orders=False,
        rpa_rescale_eigs=False,
        max_atoms_per_block=None,
    ):
        self._geom_f = None
        self._coords, self._lattice = map(_array, (coords, lattice))
        self._k_grid = _array(k_grid, dtype='i4')
        self._custom_k_pts = _array(custom_k_pts)
        self._n_freq = n_freq
        self._do_rpa = do_rpa
        self._get_spectrum = get_spectrum
        self._get_rpa_orders = get_rpa_orders
        self._rpa_rescale_eigs = rpa_rescale_eigs
        self._max_atoms_per_block = max_atoms_per_block

    def __len__(self):
        return len(self._coords)

    def __enter__(self):
        self._geom_f = _lib.cmbd_init_geom(
            len(self),
            _cast('double*', self._coords),
            _cast('double*', self._lattice),
            _cast('int*', self._k_grid),
            len(self._custom_k_pts) if self._custom_k_pts is not None else 0,
            _cast('double*', self._custom_k_pts),
            self._n_freq or 0,
            self._do_rpa,
            self._get_spectrum,
            self._get_rpa_orders,
            self._rpa_rescale_eigs,
            self._max_atoms_per_block or 0,
        )
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        _lib.cmbd_destroy_geom(self._geom_f)
        self._geom_f = None

    def _check_exc(self):
        code = _array(0, dtype=int)
        origin = _ffi.new('char[50]')
        msg = _ffi.new('char[150]')
        _lib.cmbd_get_exception(self._geom_f, _cast('int*', code), origin, msg)
        if code != 0:
            raise MBDFortranException(
                int(code), _ffi.string(origin).decode(), _ffi.string(msg).decode()
            )

    @property
    def coords(self):
        """(a.u.) Atom coordinates in rows."""
        return self._coords.copy()

    @coords.setter
    def coords(self, coords):
        _lib.cmbd_update_coords(self._geom_f, _cast('double*', _array(coords)))

    @property
    def lattice(self):
        """(a.u.) Lattice vectors in rows."""
        return self._lattice.copy()

    @lattice.setter
    def lattice(self, lattice):
        _lib.cmbd_update_lattice(self._geom_f, _cast('double*', _array(lattice)))

    def has_lattice(self):
        """Whether structure is a crystal."""
        return self._lattice is not None

    @_auto_context
    def ts_energy(self, alpha_0, C6, R_vdw, sR, d=20.0, damping='fermi', force=False):
        """Calculate a TS energy.

        :param array-like alpha_0: (a.u.) atomic polarizabilities
        :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
        :param array-like R_vdw: (a.u.) atomic vdW radii
        :param float sR: TS damping parameter :math:`s_R`
        :param float d: TS damping parameter :math:`d`
        :param damping str: type of damping
        :param force bool: if True, calculate energy gradients
        """
        alpha_0, C6, R_vdw = map(_array, (alpha_0, C6, R_vdw))
        n_atoms = len(self)
        damping_f = _lib.cmbd_init_damping(
            len(self), damping.encode(), _cast('double*', R_vdw), _ffi.NULL, sR, d
        )
        res_f = _lib.cmbd_ts_energy(
            self._geom_f,
            _cast('double*', alpha_0),
            _cast('double*', C6),
            damping_f,
            force,
        )
        _lib.cmbd_destroy_damping(damping_f)
        self._check_exc()
        ene = np.empty(1)  # for some reason np.array(0) doesn't work
        gradients, lattice_gradients = 2 * [None]
        if force:
            gradients = np.zeros((n_atoms, 3))
            if self.has_lattice():
                lattice_gradients = np.zeros((3, 3))
        results = ene, gradients, lattice_gradients, None, None, None, None, None
        _lib.cmbd_get_results(res_f, *(_cast('double*', x) for x in results))
        _lib.cmbd_destroy_result(res_f)
        ene = ene.item()
        if force:
            ene = (ene, gradients)
            if self.has_lattice():
                ene += (lattice_gradients,)
        return ene

    @_auto_context
    def mbd_energy(
        self,
        alpha_0,
        C6,
        R_vdw=None,
        beta=0.0,
        a=6.0,
        sigma=None,
        damping='fermi,dip',
        variant='rsscs',
        force=False,
    ):
        r"""Calculate an MBD energy.

        :param array-like alpha_0: (a.u.) atomic polarizabilities
        :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
        :param array-like R_vdw: (a.u.) atomic vdW radii
        :param array-like sigma: (a.u.) oscillator widths
        :param float beta: MBD damping parameter :math:`\beta`
        :param float a: MBD damping parameter :math:`a`
        :param damping str: type of damping
        :param variant str: one of 'plain', 'scs', 'rsscs'
        :param force bool: if True, calculate energy gradients
        """
        alpha_0, C6, R_vdw, sigma = map(_array, (alpha_0, C6, R_vdw, sigma))
        n_atoms = len(self)
        damping_f = _lib.cmbd_init_damping(
            n_atoms,
            damping.encode(),
            _cast('double*', R_vdw),
            _cast('double*', sigma),
            beta,
            a,
        )
        args = (
            self._geom_f,
            _cast('double*', alpha_0),
            _cast('double*', C6),
            damping_f,
            force,
        )
        if variant == 'plain':
            res_f = _lib.cmbd_mbd_energy(*args)
        else:
            args = args[:1] + (variant.encode(),) + args[1:]
            res_f = _lib.cmbd_mbd_scs_energy(*args)
        _lib.cmbd_destroy_damping(damping_f)
        self._check_exc()
        ene = np.empty(1)  # for some reason np.array(0) doesn't work
        gradients, lattice_gradients = 2 * [None]
        eigs, modes, rpa_orders, eigs_k, modes_k = 5 * [None]
        if force:
            gradients = np.zeros((n_atoms, 3))
            if self.has_lattice():
                lattice_gradients = np.zeros((3, 3))
        if self._get_spectrum:
            if self.has_lattice():
                n_kpts = (
                    len(self._custom_k_pts)
                    if self._custom_k_pts is not None
                    else self._k_grid.prod()
                )
                eigs_k = np.zeros((3 * n_atoms, n_kpts), order='F')
                modes_k = np.zeros(
                    (3 * n_atoms, 3 * n_atoms, n_kpts), dtype=complex, order='F'
                )
            else:
                eigs = np.zeros(3 * n_atoms)
                modes = np.zeros((3 * n_atoms, 3 * n_atoms), order='F')
        elif self._get_rpa_orders:
            rpa_orders = np.zeros(10)
        results = (
            ene,
            gradients,
            lattice_gradients,
            eigs,
            modes,
            rpa_orders,
            eigs_k,
            modes_k,
        )
        _lib.cmbd_get_results(res_f, *(_cast('double*', x) for x in results))
        _lib.cmbd_destroy_result(res_f)
        ene = ene.item()
        if self.has_lattice():
            eigs, modes = eigs_k, modes_k
        if self._get_spectrum:
            ene = ene, eigs, modes
        elif self._get_rpa_orders:
            ene = ene, rpa_orders
        if force:
            ene = (ene, gradients)
            if self.has_lattice():
                ene += (lattice_gradients,)
        return ene

    @_auto_context
    def dipole_matrix(
        self, damping, beta=0.0, k_point=None, R_vdw=None, sigma=None, a=6.0
    ):  # noqa: D102
        R_vdw, sigma, k_point = map(_array, (R_vdw, sigma, k_point))
        n_atoms = len(self)
        damping_f = _lib.cmbd_init_damping(
            n_atoms,
            damping.encode(),
            _cast('double*', R_vdw),
            _cast('double*', sigma),
            beta,
            a,
        )
        dipmat = np.empty(
            (3 * n_atoms, 3 * n_atoms), dtype=float if k_point is None else complex
        )
        _lib.cmbd_dipole_matrix(
            self._geom_f, damping_f, _cast('double*', k_point), _cast('double*', dipmat)
        )
        _lib.cmbd_destroy_damping(damping_f)
        self._check_exc()
        return dipmat

    def mbd_energy_species(self, species, volume_ratios, beta, **kwargs):
        r"""Calculate an MBD energy from atom types and Hirshfed-volume ratios.

        :param array-like species: atom types (elements)
        :param array-like volume_ratios: ratios of Hirshfeld volumes in
            molecule and vacuum
        :param float beta: MBD damping parameter :math:`\beta`
        :param kwargs: see :meth:`mbd_energy`
        """
        alpha_0, C6, R_vdw = from_volumes(species, volume_ratios)
        return self.mbd_energy(alpha_0, C6, R_vdw, beta, **kwargs)

    def ts_energy_species(self, species, volume_ratios, beta, **kwargs):
        """Calculate a TS energy from atom types and Hirshfed-volume ratios.

        :param array-like species: atom types (elements)
        :param array-like volume_ratios: ratios of Hirshfeld volumes in
            molecule and vacuum
        :param float sR: TS damping parameter :math:`s_R`
        :param kwargs: see :meth:`ts_energy`
        """
        alpha_0, C6, R_vdw = from_volumes(species, volume_ratios)
        return self.ts_energy(alpha_0, C6, R_vdw, beta, **kwargs)

    @_auto_context
    def dipole_energy(self, a0, w, w_t, version, r_vdw, beta, a, C):  # noqa: D102
        n_atoms = len(self)
        res = _lib.cmbd_dipole_energy(
            self._geom_f,
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

    @_auto_context
    def coulomb_energy(self, q, m, w_t, version, r_vdw, beta, a, C):  # noqa: D102
        n_atoms = len(self)
        res = _lib.cmbd_coulomb_energy(
            self._geom_f,
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
    buffer_size = (np.prod(shape) if shape else 1) * np.dtype(dtype).itemsize
    return np.ndarray(buffer=_ffi.buffer(ptr, buffer_size), shape=shape, dtype=dtype)


def _cast(ctype, array):
    return _ffi.NULL if array is None else _ffi.cast(ctype, array.ctypes.data)
