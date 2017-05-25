# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
import json
import sys
import numpy as np

from .lib import mbd as lib
from .vdw_param import vdw_param as free_atom_db

try:
    from mpi4py import MPI
    ntasks = MPI.COMM_WORLD.Get_size()
    myid = MPI.COMM_WORLD.Get_rank()
except ImportError:
    ntasks = 1
    myid = 0

lib.my_task = myid
lib.n_tasks = ntasks

bohr = lib.bohr


def get_free_atom_data(species):
    return list(map(
        list,
        zip(*[
            (at['alpha_0'], at['C6'], at['R_vdw']) for at in
            [free_atom_db[sp] for sp in species]
        ])
    ))


def printout(s, each=False):
    if each or myid == 0:
        sys.stdout.write('{}\n'.format(s))


def printerr(s, each=False):
    if each or myid == 0:
        sys.stderr.write('{}\n'.format(s))


def get_damping(xc):
    return dict(zip(['ts_d', 'ts_s_r', 'mbd_scs_a', 'mbd_ts_a',
                     'mbd_ts_erf_beta', 'mbd_ts_fermi_beta', 'mbd_rsscs_a',
                     'mbd_rsscs_beta'],
                    lib.get_damping_parameters(xc)))


def scale_hirsh(hirsh, alpha, C6, R_vdw):
    hirsh = np.array(hirsh)
    return [np.array(q)*hirsh**factor
            for q, factor in zip([alpha, C6, R_vdw], [1, 2, 1/3])]


def mbd_rsscs(
        coords, species, volumes, beta, lattice=None, k_grid=None, supercell=None, vacuum=None,
        custom_params=None
):
    lib.param_vacuum_axis = [False]*3 if vacuum is None else vacuum
    mode = 'P' if ntasks > 1 else ''
    params = get_free_atom_data(species)
    if custom_params:
        for i, specie in enumerate(species):
            if specie in custom_params:
                for j, paramname in enumerate(['alpha_0', 'C6', 'R_vdw']):
                    params[j][i] = custom_params[specie][paramname]
    alpha_0, C6, R_vdw = scale_hirsh(volumes, *params)
    if lattice is not None:
        mode += 'C'
    alpha_scs_dyn = lib.run_scs(
        mode, 'fermi,dip,gg', coords,
        lib.alpha_dynamic_ts_all('C', lib.n_grid_omega, alpha_0, c6=C6),
        unit_cell=lattice,
        r_vdw=R_vdw, beta=beta, a=6.
    )
    C6_scs = lib.get_c6_from_alpha(alpha_scs_dyn)
    R_vdw_scs = R_vdw*(alpha_scs_dyn[0]/alpha_0)**(1/3)
    if k_grid is not None:
        mode = mode.replace('C', 'R')
        k_grid = lib.make_k_grid(lib.make_g_grid(*k_grid), lattice)
        ene = lib.get_reciprocal_mbd_energy(
            mode, 'fermi,dip', coords,
            alpha_scs_dyn[0],
            lib.omega_eff(C6_scs, alpha_scs_dyn[0]),
            k_grid,
            unit_cell=lattice,
            r_vdw=R_vdw_scs, beta=beta, a=6.
        )[0]
    elif supercell is not None:
        ene = lib.get_supercell_mbd_energy(
            mode, 'fermi,dip', coords,
            alpha_scs_dyn[0],
            lib.omega_eff(C6_scs, alpha_scs_dyn[0]),
            unit_cell=lattice,
            supercell=supercell,
            r_vdw=R_vdw_scs, beta=beta, a=6.
        )[0]
    else:
        ene = lib.get_single_mbd_energy(
            mode, 'fermi,dip', coords,
            alpha_scs_dyn[0],
            lib.omega_eff(C6_scs, alpha_scs_dyn[0]),
            r_vdw=R_vdw_scs, beta=beta, a=6.
        )[0]
    return ene


def mbd_rsscs_deriv(
        coords, species, volumes, beta, lattice=None, k_grid=None, supercell=None,
        delta=1e-4
):
    forces = np.zeros(coords.shape)
    for i_atom in range(coords.shape[0]):
        for i_xyz in range(3):
            ene = []
            for step in [-1, 1]:
                coords_diff = coords.copy()
                coords_diff[i_atom, i_xyz] += step*delta
                ene.append(mbd_rsscs(
                    coords_diff, species, volumes, beta, lattice, k_grid, supercell
                ))
            forces[i_atom, i_xyz] = -(ene[1]-ene[0])/(2*delta)
    if lattice is None:
        return forces
    stress = np.zeros((3, 3))
    for i_xyz in range(3):
        for j_xyz in range(3):
            ene = []
            for step in [-1, 1]:
                strain = np.eye(3)
                strain[i_xyz, j_xyz] += step*delta
                coords_diff = coords.dot(strain)
                lattice_diff = lattice.dot(strain)
                ene.append(mbd_rsscs(
                    coords_diff, species, volumes, beta, lattice_diff, k_grid, supercell
                ))
            stress[i_xyz, j_xyz] = -(ene[1]-ene[0])/(2*delta)
    return forces, stress


class ArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        try:
            return obj.tolist()
        except AttributeError:
            return super().default(obj)
