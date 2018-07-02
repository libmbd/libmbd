# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import print_function
from __future__ import division
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erf, erfc
import pkg_resources
import sys
import csv
from itertools import product
from ._libmbd import ffi as _ffi, lib as _lib
if _lib.with_scalapack:
    from mpi4py import MPI  # noqa

with_scalapack = _lib.with_scalapack
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

    def pymbd_energy(self, coords, alpha_0, C6, R_vdw, beta,
                     lattice=None, k_grid=None, nfreq=15):
        coords = _array(coords, dtype=float)
        alpha_0 = _array(alpha_0, dtype=float)
        C6 = _array(C6, dtype=float)
        R_vdw = _array(R_vdw, dtype=float)
        freq, freq_w = freq_grid(nfreq)
        omega = 4/3*C6/alpha_0**2
        alpha_dyn = [alpha_0/(1+(u/omega)**2) for u in freq]
        alpha_dyn_rsscs = []
        for a in alpha_dyn:
            sigma = (np.sqrt(2/np.pi)*a/3)**(1/3)
            dipmat = dipole_matrix(
                coords, 'fermi,dip,gg', sigma=sigma, R_vdw=R_vdw, beta=beta,
                lattice=lattice
            )
            a_nlc = np.linalg.inv(np.diag(np.repeat(1/a, 3))+dipmat)
            a_contr = sum(np.sum(a_nlc[i::3, i::3], 1) for i in range(3))/3
            alpha_dyn_rsscs.append(a_contr)
        alpha_dyn_rsscs = np.stack(alpha_dyn_rsscs)
        C6_rsscs = 3/np.pi*np.sum(freq_w[:, None]*alpha_dyn_rsscs**2, 0)
        R_vdw_rsscs = R_vdw*(alpha_dyn_rsscs[0, :]/alpha_0)**(1/3)
        omega_rsscs = 4/3*C6_rsscs/alpha_dyn_rsscs[0, :]**2
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
                np.outer(pre, pre)*dipole_matrix(
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

    def ts_energy_species(self, coords, species, volumes, sR, **kwargs):
        alpha_0, C6, R_vdw = from_volumes(species, volumes)
        return self.ts_energy(coords, alpha_0, C6, R_vdw, sR, **kwargs)


def freq_grid(n, L=0.6):
    x, w = leggauss(n)
    w = 2*L/(1-x)**2*w
    x = L*(1+x)/(1-x)
    return np.hstack(([0], x[::-1])), np.hstack(([0], w[::-1]))


def dipole_matrix(coords, damping, beta=0., lattice=None, k_point=None,
                  R_vdw=None, sigma=None, a=6.):
    if lattice is not None:
        volume = max(np.abs(np.product(np.linalg.eigvals(lattice))), 0.2)
        ewald_alpha = 2.5/volume**(1/3)
        real_space_cutoff = 6/ewald_alpha
        range_cell = supercell_circum(lattice, real_space_cutoff)
    else:
        range_cell = (0, 0, 0)
    do_ewald = lattice is not None and damping in {'fermi,dip'}
    n = len(coords)
    dtype = float if k_point is None else complex
    dipmat = np.zeros((n, n, 3, 3), dtype=dtype)
    for idx_cell in product(*(range(-i, i+1) for i in range_cell)):
        R_cell = lattice.T.dot(idx_cell) if lattice is not None else np.zeros(3)
        Rs = coords[:, None, :]-coords[None, :, :]+R_cell
        dists = np.sqrt(np.sum(Rs**2, -1))
        if R_vdw is not None:
            S_vdw = beta*(R_vdw[:, None]+R_vdw[None, :])
        if sigma is not None:
            sigma_ij = np.sqrt(sigma[:, None]**2+sigma[None, :]**2)
        if damping == 'fermi,dip':
            T = damping_fermi(dists, S_vdw, a)[:, :, None, None]*T_bare(Rs)
        elif damping == 'fermi,dip,gg':
            T = (1-damping_fermi(dists, S_vdw, a)[:, :, None, None]) * \
                T_erf_coulomb(Rs, sigma_ij)
        else:
            raise ValueError('Unsupported damping: {0}'.format(damping))
        if do_ewald:
            T += T_erfc(Rs, ewald_alpha)-T_bare(Rs)
        if k_point is not None:
            k_pref = np.exp(-1j*np.sum(k_point*Rs, -1))
            T = k_pref[:, :, None, None]*T
        dipmat += T
    if do_ewald:
        dipmat += dipole_matrix_ewald(coords, lattice, ewald_alpha, k_point)
    n_atoms = np.shape(coords)[0]
    return np.reshape(
        np.transpose(dipmat, (0, 2, 1, 3)),
        (3*n_atoms, 3*n_atoms)
    )


def dipole_matrix_ewald(coords, lattice, alpha, k_point=None):
    Rs = coords[:, None, :]-coords[None, :, :]
    rlattice = 2*np.pi*np.linalg.inv(lattice.T)
    volume = abs(np.product(np.linalg.eigvals(lattice)))
    rec_space_cutoff = 10*alpha
    range_G_vector = supercell_circum(rlattice, rec_space_cutoff)
    dtype = float if k_point is None else complex
    dipmat = np.zeros((len(Rs), len(Rs), 3, 3), dtype=dtype)
    fourier_factor = \
        (lambda x: np.cos(x)) if k_point is None else (lambda x: np.exp(1j*x))
    for idx_gvec in product(*(range(-i, i+1) for i in range_G_vector)):
        if idx_gvec == (0, 0, 0):
            continue
        gvec = rlattice.T.dot(idx_gvec)
        k_total = k_point+gvec if k_point is not None else gvec
        k_sq = sum(k_total**2)
        if np.sqrt(k_sq) > rec_space_cutoff:
            continue
        k_prefactor = 4*np.pi/volume*np.exp(-k_sq/(4*alpha**2)) * \
            np.outer(k_total, k_total)/k_sq
        dipmat += \
            k_prefactor*fourier_factor(np.sum(gvec*Rs, -1))[:, :, None, None]
    dipmat += -np.eye(len(Rs))[:, :, None, None] * \
        np.diag(np.repeat(4*alpha**3/(3*np.sqrt(np.pi)), 3))
    k_sq = np.sum(k_point**2) if k_point is not None else 0
    if np.sqrt(k_sq) > 1e-15:
        dipmat += 4*np.pi/volume*np.outer(k_point, k_point)/k_sq * \
            np.exp(-k_sq/(4*alpha**2))
    else:
        dipmat += np.diag(np.repeat(4*np.pi/(3*volume), 3))
    return dipmat


def supercell_circum(latt, radius):
    rlatt = 2*np.pi*np.linalg.inv(latt.T)
    layer_sep = np.sum(latt*rlatt/np.sqrt(np.sum(rlatt**2, 1))[None, :], 0)
    return np.ceil(radius/layer_sep+0.5).astype(int)


def damping_fermi(R, S_vdw, d):
    return 1/(1+np.exp(-d*(R/S_vdw-1)))


def T_bare(R):
    R_2 = np.sum(R**2, -1)
    R_5 = np.where(R_2 > 0, np.sqrt(R_2)**5, np.inf)
    return (
        -3*R[:, :, :, None]*R[:, :, None, :] +
        R_2[:, :, None, None]*np.eye(3)[None, None, :, :]
    )/R_5[:, :, None, None]


def T_erf_coulomb(R, sigma):
    bare = T_bare(R)
    R_1 = np.sqrt(np.sum(R**2, -1))
    R_5 = np.where(R_1 > 0, R_1**5, np.inf)
    RR_R5 = R[:, :, :, None]*R[:, :, None, :]/R_5[:, :, None, None]
    zeta = R_1/sigma
    theta = 2*zeta/np.sqrt(np.pi)*np.exp(-zeta**2)
    erf_theta = erf(zeta)-theta
    return erf_theta[:, :, None, None]*bare + \
        (2*(zeta**2)*theta)[:, :, None, None]*RR_R5


def T_erfc(R, a):
    R_2 = np.sum(R**2, -1)
    R_1 = np.sqrt(R_2)
    R_3 = np.where(R_1 > 0, R_1**3, np.inf)
    R_5 = np.where(R_1 > 0, R_1**5, np.inf)
    B = (erfc(a*R_1)+(2*a*R_1/np.sqrt(np.pi))*np.exp(-(a*R_1)**2))/R_3
    C = (3*erfc(a*R_1)+(2*a*R_1/np.sqrt(np.pi))*(3+2*(a*R_1)**2)*np.exp(-(a*R_1)**2))/R_5
    return -C[:, :, None, None]*R[:, :, :, None]*R[:, :, None, :] + \
        B[:, :, None, None]*np.eye(3)


def from_volumes(species, volumes, kind='TS'):
    alpha_0, C6, R_vdw = (
        np.array([vdw_params[sp][param] for sp in species])
        for param in ['alpha_0({})'.format(kind), 'C6({})'.format(kind), 'R_vdw']
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


def _get_vdw_params():
    csv_lines = pkg_resources.resource_string(__name__, 'vdw-params.csv').split(b'\n')
    if sys.version_info[0] > 2:
        csv_lines = [l.decode() for l in csv_lines]
    reader = csv.DictReader(csv_lines, quoting=csv.QUOTE_NONNUMERIC)
    vdw_params = {}
    for row in reader:
        vdw_params[row.pop('symbol')] = row
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
    )
