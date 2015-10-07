#!/usr/bin/env python3
"""pymbd -- Pythom wrapper for MBD.

Usage:
    pymbd.py INPUT [-e EXTENSION]

Options:
    -e, --extension EXTENSION      Format "MODULE:CALLABLE".
"""
from __future__ import print_function
import json
import numpy as np
from mpi4py import MPI
from mbd import mbd
import sys
from contextlib import contextmanager
import os


bohr = 0.529177249
ntasks = MPI.COMM_WORLD.Get_size()
myid = MPI.COMM_WORLD.Get_rank()

free_atom_db = json.load(open(os.path.join(sys.path[0], 'free_atoms.json')))


def printmsg(s):
    if myid == 0:
        sys.stderr.write('{}\n'.format(s))


@contextmanager
def block(name):
    printmsg('Evaluating {}...'.format(name))
    yield


def run_mbd(data, mbd):
    natoms = len(data['geometry'])
    nomega = mbd.omega_grid.shape[0]
    alpha, C6, R_vdw, energy = {}, {}, {}, {}
    damp_params = dict(zip(['ts_d', 'ts_s_r', 'mbd_scs_a', 'mbd_ts_a',
                            'mbd_ts_erf_beta', 'mbd_ts_fermi_beta', 'mbd_rsscs_a',
                            'mbd_rsscs_beta'],
                           mbd.get_damping_parameters('pbe')))
    geom = data['geometry']
    xyz = [np.array(a[1:4])/bohr for a in geom]
    free_atoms = [free_atom_db[a[0]] for a in geom]
    alpha_0_free = [a['alpha_0'] for a in free_atoms]
    C6_free = [a['C6'] for a in free_atoms]
    R_vdw_free = [a['R_vdw'] for a in free_atoms]
    volume_ratio = np.array(data['volume_ratio'])
    unit_cell = np.array(data.get('unit_cell', np.zeros((3, 3))))/bohr
    if unit_cell.any():
        k_pts = mbd.make_k_grid(mbd.make_g_grid(*data['k_grid']), unit_cell)
    else:
        k_pts = [(0, 0, 0)]
    with block('TS'):
        alpha['TS'] = np.zeros((nomega, natoms))
        alpha['TS'][0] = alpha_0_free*volume_ratio
        C6['TS'] = C6_free*volume_ratio**2
        for alphas, omega in zip(alpha['TS'][1:], mbd.omega_grid[1:]):
            alphas[:] = mbd.alpha_dynamic_ts(alpha['TS'][0], C6['TS'], omega)
        R_vdw['TS'] = R_vdw_free*volume_ratio**(1./3)
        energy['TS@TS~fermi@TS'] = \
            mbd.get_ts_energy(xyz,
                              C6['TS'],
                              alpha['TS'][0],
                              'fermi',
                              R_vdw['TS'],
                              d=damp_params['ts_d'],
                              s_r=damp_params['ts_s_r'],
                              unit_cell=unit_cell,
                              my_task=myid, n_tasks=ntasks)
    # with block('MBD@TS'):
    #     energy['MBD@TS~erf@TS,dip'] = \
    #         mbd.get_periodic_mbd_energy(xyz,
    #                                     alpha['TS'][0],
    #                                     mbd.omega_eff(C6['TS'], alpha['TS'][0]),
    #                                     'erf,dip',
    #                                     unit_cell,
    #                                     k_pts,
    #                                     R_vdw['TS'],
    #                                     beta=damp_params['mbd_ts_erf_beta'],
    #                                     a=4.,
    #                                     my_task=myid, n_tasks=ntasks)[0]
    #     energy['MBD@TS~fermi@TS,dip'] = \
    #         mbd.get_periodic_mbd_energy(xyz,
    #                                     alpha['TS'][0],
    #                                     mbd.omega_eff(C6['TS'], alpha['TS'][0]),
    #                                     'fermi@TS,dip',
    #                                     unit_cell,
    #                                     k_pts,
    #                                     R_vdw['TS'],
    #                                     beta=damp_params['mbd_ts_fermi_beta'],
    #                                     a=damp_params['mbd_ts_a'],
    #                                     my_task=myid, n_tasks=ntasks)[0]
    # with block('SCS'):
    #     alpha['SCS'] = mbd.run_scs(xyz,
    #                                alpha['TS'],
    #                                'dip,gg',
    #                                unit_cell=unit_cell,
    #                                my_task=myid, n_tasks=ntasks)
    #     C6['SCS'] = mbd.get_c6_from_alpha(alpha['SCS'])
    #     R_vdw['SCS'] = \
    #         R_vdw_free*(alpha['SCS'][0]/alpha_0_free)**(1./3)
    # with block('TS@SCS'):
    #     energy['TS@SCS~fermi@SCS'] = \
    #         mbd.get_ts_energy(xyz,
    #                           C6['SCS'],
    #                           alpha['SCS'][0],
    #                           'fermi',
    #                           R_vdw['SCS'],
    #                           d=damp_params['ts_d'],
    #                           s_r=damp_params['ts_s_r'],
    #                           unit_cell=unit_cell,
    #                           my_task=myid, n_tasks=ntasks)
    # with block('MBD@SCS'):
    #     energy['MBD@SCS~dip,1mexp@SCS'] = \
    #         mbd.get_periodic_mbd_energy(xyz,
    #                                     alpha['SCS'][0],
    #                                     mbd.omega_eff(C6['SCS'], alpha['SCS'][0]),
    #                                     'dip,1mexp',
    #                                     unit_cell,
    #                                     k_pts,
    #                                     R_vdw['SCS'],
    #                                     beta=1.,
    #                                     a=damp_params['mbd_scs_a'],
    #                                     my_task=myid, n_tasks=ntasks)[0]
    with block('rsSCS'):
        alpha['rsSCS'] = mbd.run_scs(xyz,
                                     alpha['TS'],
                                     'fermi,dip,gg',
                                     R_vdw['TS'],
                                     beta=damp_params['mbd_rsscs_beta'],
                                     a=damp_params['mbd_rsscs_a'],
                                     unit_cell=unit_cell,
                                     my_task=myid, n_tasks=ntasks)
        C6['rsSCS'] = mbd.get_c6_from_alpha(alpha['rsSCS'])
        R_vdw['rsSCS'] = \
            R_vdw_free*(alpha['rsSCS'][0]/alpha_0_free)**(1./3)
    with block('MBD@rsSCS'):
        energy['MBD@rsSCS~fermi@rsSCS,dip'], mode_enes, modes = \
            mbd.get_periodic_mbd_energy(xyz,
                                        alpha['rsSCS'][0],
                                        mbd.omega_eff(C6['rsSCS'], alpha['rsSCS'][0]),
                                        'fermi@rsSCS,dip',
                                        unit_cell,
                                        k_pts,
                                        R_vdw['rsSCS'],
                                        beta=damp_params['mbd_rsscs_beta'],
                                        a=damp_params['mbd_rsscs_a'],
                                        my_task=myid, n_tasks=ntasks)
        energy['MBD(TS)@rsSCS~fermi@rsSCS,dip'] = \
            mbd.get_ts_energy(xyz,
                              C6['rsSCS'],
                              alpha['rsSCS'][0],
                              'fermi2',
                              R_vdw['rsSCS'],
                              s_r=damp_params['mbd_rsscs_beta'],
                              d=damp_params['mbd_rsscs_a'],
                              unit_cell=unit_cell,
                              my_task=myid, n_tasks=ntasks)
        # if not unit_cell.any():
        #     rpa_ene, rpa_orders = \
        #         mbd.get_qho_rpa_energy(xyz,
        #                                alpha['rsSCS'],
        #                                'fermi@rsSCS,dip',
        #                                R_vdw['rsSCS'],
        #                                beta=damp_params['mbd_rsscs_beta'],
        #                                a=damp_params['mbd_rsscs_a'],
        #                                my_task=myid, n_tasks=ntasks)
        #     rpa_orders[0] = rpa_ene
        #     energy['MBD(RPA)@rsSCS~fermi@rsSCS,dip'] = rpa_orders[:10]
        #     energy['MBD(nbody)@rsSCS~fermi@rsSCS,dip'] = \
        #         mbd.mbd_nbody(xyz,
        #                       alpha['rsSCS'][0],
        #                       mbd.omega_eff(C6['rsSCS'], alpha['rsSCS'][0]),
        #                       'fermi@rsSCS,dip',
        #                       R_vdw['rsSCS'],
        #                       beta=damp_params['mbd_rsscs_beta'],
        #                       a=damp_params['mbd_rsscs_a'],
        #                       my_task=myid, n_tasks=ntasks)[:3]
    return energy


def main(data, extension=None):
    if extension:
        module, func = extension.split(':')
        module = __import__(module)
        extension = getattr(module, func)
    else:
        extension = run_mbd
    printmsg('Running on {} nodes...'.format(ntasks))
    results = extension(data, mbd)
    return results


class ArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        try:
            return obj.tolist()
        except AttributeError:
            return super().default(obj)


if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)
    results = main(json.load(open(args['INPUT'])),
                   args['--extension'])
    if myid == 0:
        json.dump(results, sys.stdout, cls=ArrayEncoder, indent=4)
