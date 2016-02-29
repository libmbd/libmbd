#!/usr/bin/env python3
import json
from mpi4py import MPI
import sys
from pathlib import Path
from mbd import mbd
import numpy as np

ntasks = MPI.COMM_WORLD.Get_size()
myid = MPI.COMM_WORLD.Get_rank()
mbd.my_task = myid
mbd.n_tasks = ntasks

free_atom_db = json.load((Path(__file__).parent/'free_atoms.json').open())
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


def mbd_rsscs(geom, alpha_0, C6, R_vdw, beta, a):
    xyz = [a.xyz/bohr for a in geom]
    mode = 'M' if ntasks > 1 else ''
    alpha_scs_dyn = mbd.run_scs(
        mode, 'fermi,dip,gg', xyz,
        mbd.alpha_dynamic_ts_all('C', mbd.n_grid_omega, alpha_0, c6=C6),
        r_vdw=R_vdw, beta=beta, a=a)
    C6_scs = mbd.get_c6_from_alpha(alpha_scs_dyn)
    R_vdw_scs = R_vdw*(alpha_scs_dyn[0]/alpha_0)**(1/3)
    return mbd.get_single_mbd_energy(mode, 'fermi,dip', xyz,
                                     alpha_scs_dyn[0],
                                     mbd.omega_eff(C6_scs, alpha_scs_dyn[0]),
                                     r_vdw=R_vdw_scs, beta=beta, a=a)[0]


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
    class Context:
        pass
    script = load_run_script(sys.argv[1])
    ctx = Context()
    for key, value in dict(locals()).items():
        if not key.startswith('_'):
            setattr(ctx, key, value)
    script.run(ctx, mbd)
