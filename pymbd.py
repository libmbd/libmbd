#!/usr/bin/env python3
import os
import json
import sys
import numpy as np
from mbd import mbd

try:
    from mpi4py import MPI
    ntasks = MPI.COMM_WORLD.Get_size()
    myid = MPI.COMM_WORLD.Get_rank()
except ImportError:
    sys.stderr.write('warning: Install mpi4py for MPI support\n')
    ntasks = 1
    myid = 0

mbd.my_task = myid
mbd.n_tasks = ntasks
with open(os.path.join(os.path.dirname(__file__), 'free_atoms.json')) as f:
    free_atom_db = json.load(f)
bohr = mbd.bohr


def get_free_atom_data(species):
    return list(zip(*[(at['alpha_0'], at['C6'], at['R_vdw']) for at in
                      [free_atom_db[sp] for sp in species]]))


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
                    mbd.get_damping_parameters(xc)))


def scale_hirsh(hirsh, alpha, C6, R_vdw):
    hirsh = np.array(hirsh)
    return [np.array(q)*hirsh**factor
            for q, factor in zip([alpha, C6, R_vdw], [1, 2, 1/3])]


def mbd_rsscs(
        coords, species, volumes, beta, lattice=None, k_grid=None, supercell=None, vacuum=None,
        custom_params=None
):
    mbd.param_vacuum_axis = [False]*3 if vacuum is None else vacuum
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
    alpha_scs_dyn = mbd.run_scs(
        mode, 'fermi,dip,gg', coords,
        mbd.alpha_dynamic_ts_all('C', mbd.n_grid_omega, alpha_0, c6=C6),
        unit_cell=lattice,
        r_vdw=R_vdw, beta=beta, a=6.
    )
    C6_scs = mbd.get_c6_from_alpha(alpha_scs_dyn)
    R_vdw_scs = R_vdw*(alpha_scs_dyn[0]/alpha_0)**(1/3)
    if k_grid is not None:
        mode = mode.replace('C', 'R')
        k_grid = mbd.make_k_grid(mbd.make_g_grid(*k_grid), lattice)
        ene = mbd.get_reciprocal_mbd_energy(
            mode, 'fermi,dip', coords,
            alpha_scs_dyn[0],
            mbd.omega_eff(C6_scs, alpha_scs_dyn[0]),
            k_grid,
            unit_cell=lattice,
            r_vdw=R_vdw_scs, beta=beta, a=6.
        )[0]
    elif supercell is not None:
        ene = mbd.get_supercell_mbd_energy(
            mode, 'fermi,dip', coords,
            alpha_scs_dyn[0],
            mbd.omega_eff(C6_scs, alpha_scs_dyn[0]),
            unit_cell=lattice,
            supercell=supercell,
            r_vdw=R_vdw_scs, beta=beta, a=6.
        )[0]
    else:
        ene = mbd.get_single_mbd_energy(
            mode, 'fermi,dip', coords,
            alpha_scs_dyn[0],
            mbd.omega_eff(C6_scs, alpha_scs_dyn[0]),
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


def load_run_script(path):
    import imp
    script = imp.new_module('script')
    try:
        exec(compile(open(path).read(), path, 'exec'), script.__dict__)
    except:
        import traceback
        traceback.print_exc()
        printerr('There was an error while reading run script.')
    return script


if __name__ == '__main__':
    if len(sys.argv) > 1:
        class Context:
            pass
        script = load_run_script(sys.argv[1])
        ctx = Context()
        for key, value in dict(locals()).items():
            if not key.startswith('_'):
                setattr(ctx, key, value)
        script.run(ctx, mbd)
    else:
        mbd.init_grid(15)
        print(mbd_rsscs(
            [[0, 0, 0], [4.0/bohr, 0, 0]],
            ['Ar', 'Ar'],
            [1., 1.],
            0.83
        ))
