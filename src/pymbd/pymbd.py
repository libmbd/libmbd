# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import division, print_function

import csv
import sys
from itertools import product

import numpy as np
from numpy.polynomial.legendre import leggauss
from pkg_resources import resource_string
from scipy.special import erf, erfc

__all__ = ['mbd_energy', 'mbd_energy_species', 'screening', 'ang']

ang = 1 / 0.529177249
"""(a.u.) angstrom"""


def screening(coords, alpha_0, C6, R_vdw, beta, lattice=None, nfreq=15):
    r"""Screen atomic polarizabilities.

    :param array-like coords: (a.u.) atom coordinates in rows
    :param array-like alpha_0: (a.u.) atomic polarizabilities
    :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
    :param array-like R_vdw: (a.u.) atomic vdW radii
    :param float beta: MBD damping parameter :math:`\beta`
    :param array-like lattice: (a.u.) lattice vectors in rows
    :param int nfreq: number of grid points for frequency quadrature

    Returns static polarizabilities, :math:`C_6` coefficients, and
    :math:`R_\mathrm{vdw}` coefficients (a.u.).
    """
    freq, freq_w = freq_grid(nfreq)
    omega = 4 / 3 * C6 / alpha_0 ** 2
    alpha_dyn = [alpha_0 / (1 + (u / omega) ** 2) for u in freq]
    alpha_dyn_rsscs = []
    for a in alpha_dyn:
        sigma = (np.sqrt(2 / np.pi) * a / 3) ** (1 / 3)
        dipmat = dipole_matrix(
            coords, 'fermi,dip,gg', sigma=sigma, R_vdw=R_vdw, beta=beta, lattice=lattice
        )
        a_nlc = np.linalg.inv(np.diag(np.repeat(1 / a, 3)) + dipmat)
        a_contr = sum(np.sum(a_nlc[i::3, i::3], 1) for i in range(3)) / 3
        alpha_dyn_rsscs.append(a_contr)
    alpha_dyn_rsscs = np.stack(alpha_dyn_rsscs)
    C6_rsscs = 3 / np.pi * np.sum(freq_w[:, None] * alpha_dyn_rsscs ** 2, 0)
    R_vdw_rsscs = R_vdw * (alpha_dyn_rsscs[0, :] / alpha_0) ** (1 / 3)
    return alpha_dyn_rsscs[0], C6_rsscs, R_vdw_rsscs


def mbd_energy(coords, alpha_0, C6, R_vdw, beta, lattice=None, k_grid=None, nfreq=15):
    r"""Calculate an MBD energy.

    :param array-like coords: (a.u.) atom coordinates in rows
    :param array-like alpha_0: (a.u.) atomic polarizabilities
    :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
    :param array-like R_vdw: (a.u.) atomic vdW radii
    :param float beta: MBD damping parameter :math:`\beta`
    :param array-like lattice: (a.u.) lattice vectors in rows
    :param array-like k_grid: number of :math:`k`-points along reciprocal axes
    :param int nfreq: number of grid points for frequency quadrature
    """
    coords, alpha_0, C6, R_vdw, lattice = map(
        _array, (coords, alpha_0, C6, R_vdw, lattice)
    )
    alpha_0_rsscs, C6_rsscs, R_vdw_rsscs = screening(
        coords, alpha_0, C6, R_vdw, beta, lattice=lattice, nfreq=15
    )
    omega_rsscs = 4 / 3 * C6_rsscs / alpha_0_rsscs ** 2
    pre = np.repeat(omega_rsscs * np.sqrt(alpha_0_rsscs), 3)
    if lattice is None:
        k_points = [None]
    else:
        assert k_grid is not None
        k_points = get_kpts(lattice, k_grid)
    ene = 0.0
    for k_point in k_points:
        eigs = np.linalg.eigvalsh(
            np.diag(np.repeat(omega_rsscs ** 2, 3))
            + np.outer(pre, pre)
            * dipole_matrix(
                coords,
                'fermi,dip',
                R_vdw=R_vdw_rsscs,
                beta=beta,
                lattice=lattice,
                k_point=k_point,
            )
        )
        ene += np.sum(np.sqrt(eigs)) / 2
    ene = ene / len(k_points) - 3 * np.sum(omega_rsscs) / 2
    return ene


def mbd_energy_species(coords, species, volume_ratios, beta, **kwargs):
    r"""Calculate an MBD energy from atom types and Hirshfed-volume ratios.

    :param array-like coords: (a.u.) atom coordinates in rows
    :param array-like species: atom types (elements)
    :param array-like volume_ratios: ratios of Hirshfeld volumes in molecule and vacuum
    :param float beta: MBD damping parameter :math:`\beta`
    :param kwargs: see :func:`mbd_energy`
    """
    alpha_0, C6, R_vdw = from_volumes(species, volume_ratios)
    return mbd_energy(coords, alpha_0, C6, R_vdw, beta, **kwargs)


def dipole_matrix(
    coords, damping, beta=0.0, lattice=None, k_point=None, R_vdw=None, sigma=None, a=6.0
):
    if lattice is not None:
        volume = max(np.abs(np.product(np.linalg.eigvals(lattice))), 0.2)
        ewald_alpha = 2.5 / volume ** (1 / 3)
        real_space_cutoff = 6 / ewald_alpha
        range_cell = supercell_circum(lattice, real_space_cutoff)
    else:
        range_cell = (0, 0, 0)
    do_ewald = lattice is not None and damping in {'fermi,dip'}
    n = len(coords)
    dtype = float if k_point is None else complex
    dipmat = np.zeros((n, n, 3, 3), dtype=dtype)
    if R_vdw is not None:
        S_vdw = beta * (R_vdw[:, None] + R_vdw[None, :])
    if sigma is not None:
        sigma_ij = np.sqrt(sigma[:, None] ** 2 + sigma[None, :] ** 2)
    for idx_cell in product(*(range(-i, i + 1) for i in range_cell)):
        R_cell = lattice.T.dot(idx_cell) if lattice is not None else np.zeros(3)
        Rs = coords[:, None, :] - coords[None, :, :] + R_cell
        dists = np.sqrt(np.sum(Rs ** 2, -1))
        if damping == 'fermi,dip':
            T = damping_fermi(dists, S_vdw, a)[:, :, None, None] * T_bare(Rs)
        elif damping == 'fermi,dip,gg':
            T = (1 - damping_fermi(dists, S_vdw, a)[:, :, None, None]) * T_erf_coulomb(
                Rs, sigma_ij
            )
        else:
            raise ValueError(f'Unsupported damping: {damping}')
        if do_ewald:
            T += T_erfc(Rs, ewald_alpha) - T_bare(Rs)
        if k_point is not None:
            k_pref = np.exp(-1j * np.sum(k_point * Rs, -1))
            T = k_pref[:, :, None, None] * T
        dipmat += T
    if do_ewald:
        dipmat += dipole_matrix_ewald(coords, lattice, ewald_alpha, k_point)
    n = len(coords)
    return np.reshape(np.transpose(dipmat, (0, 2, 1, 3)), (3 * n, 3 * n))


def dipole_matrix_ewald(coords, lattice, alpha, k_point=None):
    Rs = coords[:, None, :] - coords[None, :, :]
    rlattice = 2 * np.pi * np.linalg.inv(lattice.T)
    volume = abs(np.product(np.linalg.eigvals(lattice)))
    rec_space_cutoff = 10 * alpha
    range_G_vector = supercell_circum(rlattice, rec_space_cutoff)
    dtype = float if k_point is None else complex
    dipmat = np.zeros((len(Rs), len(Rs), 3, 3), dtype=dtype)
    fourier_factor = (
        (lambda x: np.cos(x)) if k_point is None else (lambda x: np.exp(1j * x))
    )
    for idx_gvec in product(*(range(-i, i + 1) for i in range_G_vector)):
        if idx_gvec == (0, 0, 0):
            continue
        gvec = rlattice.T.dot(idx_gvec)
        k_total = k_point + gvec if k_point is not None else gvec
        k_sq = sum(k_total ** 2)
        if np.sqrt(k_sq) > rec_space_cutoff:
            continue
        k_prefactor = (
            4
            * np.pi
            / volume
            * np.exp(-k_sq / (4 * alpha ** 2))
            / k_sq
            * np.outer(k_total, k_total)
        )
        dipmat += k_prefactor * fourier_factor(np.sum(gvec * Rs, -1))[:, :, None, None]
    dipmat += -np.eye(len(Rs))[:, :, None, None] * np.diag(
        np.repeat(4 * alpha ** 3 / (3 * np.sqrt(np.pi)), 3)
    )
    k_sq = np.sum(k_point ** 2) if k_point is not None else 0
    if np.sqrt(k_sq) > 1e-15:
        dipmat += (
            4
            * np.pi
            / volume
            * np.exp(-k_sq / (4 * alpha ** 2))
            / k_sq
            * np.outer(k_point, k_point)
        )
    else:
        dipmat += np.diag(np.repeat(4 * np.pi / (3 * volume), 3))
    return dipmat


def supercell_circum(latt, radius):
    rlatt = 2 * np.pi * np.linalg.inv(latt.T)
    layer_sep = np.sum(latt * rlatt / np.sqrt(np.sum(rlatt ** 2, 1))[None, :], 0)
    return np.ceil(radius / layer_sep + 0.5).astype(int)


def damping_fermi(R, S_vdw, d):
    return 1 / (1 + np.exp(-d * (R / S_vdw - 1)))


def T_bare(R):
    R_2 = np.sum(R ** 2, -1)
    R_5 = np.where(R_2 > 0, np.sqrt(R_2) ** 5, np.inf)
    return (
        -3 * R[:, :, :, None] * R[:, :, None, :]
        + R_2[:, :, None, None] * np.eye(3)[None, None, :, :]
    ) / R_5[:, :, None, None]


def T_erf_coulomb(R, sigma):
    bare = T_bare(R)
    R_1 = np.sqrt(np.sum(R ** 2, -1))
    R_5 = np.where(R_1 > 0, R_1 ** 5, np.inf)
    RR_R5 = R[:, :, :, None] * R[:, :, None, :] / R_5[:, :, None, None]
    zeta = R_1 / sigma
    theta = 2 * zeta / np.sqrt(np.pi) * np.exp(-(zeta ** 2))
    erf_theta = erf(zeta) - theta
    return (
        erf_theta[:, :, None, None] * bare
        + (2 * (zeta ** 2) * theta)[:, :, None, None] * RR_R5
    )


def T_erfc(R, a):
    R_2 = np.sum(R ** 2, -1)
    R_1 = np.sqrt(R_2)
    R_3 = np.where(R_1 > 0, R_1 ** 3, np.inf)
    R_5 = np.where(R_1 > 0, R_1 ** 5, np.inf)
    B = (
        erfc(a * R_1) + (2 * a * R_1 / np.sqrt(np.pi)) * np.exp(-((a * R_1) ** 2))
    ) / R_3
    C = (
        3 * erfc(a * R_1)
        + (2 * a * R_1 / np.sqrt(np.pi))
        * (3 + 2 * (a * R_1) ** 2)
        * np.exp(-((a * R_1) ** 2))
    ) / R_5
    return -C[:, :, None, None] * R[:, :, :, None] * R[:, :, None, :] + B[
        :, :, None, None
    ] * np.eye(3)


def from_volumes(species, volumes, kind='TS'):
    if kind == 'TS':
        alpha_0, C6, R_vdw = (
            np.array([vdw_params[sp][param] for sp in species])
            for param in 'alpha_0(TS) C6(TS) R_vdw(TS)'.split()
        )
    elif kind == 'BG':
        alpha_0, C6, R_vdw = (
            np.array([vdw_params[sp][param] for sp in species])
            for param in 'alpha_0(BG) C6(BG) R_vdw(TS)'.split()
        )
    elif kind == 'TSsurf':
        alpha_0, C6, R_vdw = (
            np.array(
                [
                    vdw_params[sp][param]
                    or vdw_params[sp][param.replace('TSsurf', 'TS')]
                    for sp in species
                ]
            )
            for param in 'alpha_0(TSsurf) C6(TSsurf) R_vdw(TSsurf)'.split()
        )
    else:
        raise ValueError(f'Unkonwn vdW parameter kind: {kind}')
    volumes = np.array(volumes)
    alpha_0 *= volumes
    C6 *= volumes ** 2
    R_vdw *= volumes ** (1 / 3)
    return alpha_0, C6, R_vdw


def get_kpts(lattice, k_grid, shift=0.5):
    k_grid, lattice = map(np.array, (k_grid, lattice))
    k_idxs = (np.array(list(product(*map(range, k_grid)))) + shift) / k_grid
    k_idxs = np.where(k_idxs > 0.5, k_idxs - 1, k_idxs)
    rlattice = 2 * np.pi * np.linalg.inv(lattice.T)
    k_points = k_idxs.dot(rlattice)
    return k_points


def freq_grid(n, L=0.6):
    x, w = leggauss(n)
    w = 2 * L / (1 - x) ** 2 * w
    x = L * (1 + x) / (1 - x)
    return np.hstack(([0], x[::-1])), np.hstack(([0], w[::-1]))


def _array(obj, *args, **kwargs):
    if obj is not None:
        kwargs.setdefault('dtype', float)
        return np.array(obj, *args, **kwargs)


def _get_vdw_params():
    csv_lines = resource_string(__name__, 'vdw-params.csv').split(b'\n')
    if sys.version_info[0] > 2:
        csv_lines = [l.decode() for l in csv_lines]
    reader = csv.DictReader(csv_lines, quoting=csv.QUOTE_NONNUMERIC)
    vdw_params = {row.pop('symbol'): row for row in reader}
    return vdw_params


vdw_params = _get_vdw_params()
