#!/usr/bin/env python
from __future__ import print_function
import json
import numpy as np
from mpi4py import MPI
from mbd import mbd
import sys
from pathlib import Path


class ArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        try:
            return obj.tolist()
        except AttributeError:
            return super().default(obj)


bohr = 0.529177249
ntasks = MPI.COMM_WORLD.Get_size()
myid = MPI.COMM_WORLD.Get_rank()


def printerr(s):
    if myid == 0:
        sys.stderr.write(s + '\n')


def block(msg):
    def runner(f):
        printerr(msg)
        f()
    return runner


def main(path, extension=None):
    data = json.load(open(path))
    if extension:
        module, func = extension.split(':')
        module = __import__(module)
        extension = getattr(module, func)
    else:
        extension = run_mbd
    for key in data:
        data[key] = np.array(data[key])
    printerr('Running on {} nodes...'.format(ntasks))
    results = extension(data, mbd)
    return results


if __name__ == '__main__':
    config = Path('config.json')
    config = json.load(open(str(config))) if config.exists() else {}
    results = main(sys.argv[1], config.get('extension'))
    if myid == 0:
        json.dump(results, sys.stdout, cls=ArrayEncoder, indent=4)


def run_mbd(data, mbd):
    natoms = len(data['coords'])
    nomega = mbd.omega_grid.shape[0]
    alpha, C6, R_vdw, energy = {}, {}, {}, {}
    damp_params = dict(zip(['ts_d', 'ts_s_r', 'mbd_scs_a', 'mbd_ts_a',
                            'mbd_ts_erf_beta', 'mbd_ts_fermi_beta', 'mbd_rsscs_a',
                            'mbd_rsscs_beta'],
                           mbd.get_damping_parameters('pbe')))

    @block('Evaluating TS...')
    def ts():
        alpha['TS'] = np.zeros((nomega, natoms))
        alpha['TS'][0] = data['alpha_0']*data['volume_ratio']
        C6['TS'] = data['C6']*np.power(data['volume_ratio'], 2)
        for alphas, omega in zip(alpha['TS'][1:], mbd.omega_grid[1:]):
            alphas[:] = mbd.alpha_dynamic_ts(alpha['TS'][0], C6['TS'], omega)
        R_vdw['TS'] = data['R_vdw']*np.power(data['volume_ratio'], 1./3)
        energy['TS@TS~fermi@TS'] = \
            mbd.get_ts_energy(data['coords'],
                              C6['TS'],
                              alpha['TS'][0],
                              'fermi',
                              R_vdw['TS'],
                              d=damp_params['ts_d'],
                              s_r=damp_params['ts_s_r'])

    @block("Evaluating MBD@TS...")
    def mbd_ts():
        energy['MBD@TS~erf@TS,dip'] = \
            mbd.get_mbd_energy(data['coords'],
                               alpha['TS'][0],
                               mbd.omega_eff(C6['TS'], alpha['TS'][0]),
                               'erf,dip',
                               R_vdw['TS'],
                               beta=damp_params['mbd_ts_erf_beta'],
                               a=4.)[0]
        energy['MBD@TS~fermi@TS,dip'] = \
            mbd.get_mbd_energy(data['coords'],
                               alpha['TS'][0],
                               mbd.omega_eff(C6['TS'], alpha['TS'][0]),
                               'fermi@TS,dip',
                               R_vdw['TS'],
                               beta=damp_params['mbd_ts_fermi_beta'],
                               a=damp_params['mbd_ts_a'])[0]

    @block("Evaluating SCS...")
    def scs():
        alpha['SCS'] = mbd.run_scs(data['coords'],
                                   alpha['TS'],
                                   'dip,gg',
                                   my_task=myid, n_tasks=ntasks)
        C6['SCS'] = mbd.get_c6_from_alpha(alpha['SCS'])
        R_vdw['SCS'] = \
            data['R_vdw']*np.power(alpha['SCS'][0]/data['alpha_0'], 1./3)

    @block("Evaluating TS@SCS...")
    def ts_scs():
        energy['TS@SCS~fermi@SCS'] = \
            mbd.get_ts_energy(data['coords'],
                              C6['SCS'],
                              alpha['SCS'][0],
                              'fermi',
                              R_vdw['SCS'],
                              d=damp_params['ts_d'],
                              s_r=damp_params['ts_s_r'])

    @block("Evaluating MBD@SCS...")
    def mbd_scs():
        energy['MBD@SCS~dip,1mexp@SCS'] = \
            mbd.get_mbd_energy(data['coords'],
                               alpha['SCS'][0],
                               mbd.omega_eff(C6['SCS'], alpha['SCS'][0]),
                               'dip,1mexp',
                               R_vdw['SCS'],
                               beta=1.,
                               a=damp_params['mbd_scs_a'])[0]

    @block("Evaluating rsSCS...")
    def rsscs():
        alpha['rsSCS'] = mbd.run_scs(data['coords'],
                                     alpha['TS'],
                                   'fermi,dip,gg',
                                   R_vdw['TS'],
                                   beta=damp_params['mbd_rsscs_beta'],
                                     a=damp_params['mbd_rsscs_a'],
                                     my_task=myid, n_tasks=ntasks)
        C6['rsSCS'] = mbd.get_c6_from_alpha(alpha['rsSCS'])
        R_vdw['rsSCS'] = \
            data['R_vdw']*np.power(alpha['rsSCS'][0]/data['alpha_0'], 1./3)

    @block("Evaluating MBD@rsSCS...")
    def mbd_rsscs():
        energy['MBD@rsSCS~fermi@rsSCS,dip'] = \
            mbd.get_mbd_energy(data['coords'],
                               alpha['rsSCS'][0],
                               mbd.omega_eff(C6['rsSCS'], alpha['rsSCS'][0]),
                               'fermi@rsSCS,dip',
                               R_vdw['rsSCS'],
                               beta=damp_params['mbd_rsscs_beta'],
                               a=damp_params['mbd_rsscs_a'])[0]
        energy['MBD(TS)@rsSCS~fermi@rsSCS,dip'] = \
            mbd.get_ts_energy(data['coords'],
                              C6['rsSCS'],
                              alpha['rsSCS'][0],
                              'fermi2',
                              R_vdw['rsSCS'],
                              s_r=damp_params['mbd_rsscs_beta'],
                              d=damp_params['mbd_rsscs_a'])
        rpa_ene, rpa_orders = \
            mbd.get_qho_rpa_energy(data['coords'],
                                   alpha['rsSCS'],
                                   'fermi@rsSCS,dip',
                                   R_vdw['rsSCS'],
                                   beta=damp_params['mbd_rsscs_beta'],
                                   a=damp_params['mbd_rsscs_a'])
        rpa_orders[0] = rpa_ene
        energy['MBD(RPA)@rsSCS~fermi@rsSCS,dip'] = rpa_orders[:10]
        energy['MBD(nbody)@rsSCS~fermi@rsSCS,dip'] = \
            mbd.mbd_nbody(data['coords'],
                          alpha['rsSCS'][0],
                          mbd.omega_eff(C6['rsSCS'], alpha['rsSCS'][0]),
                          'fermi@rsSCS,dip',
                          R_vdw['rsSCS'],
                          beta=damp_params['mbd_rsscs_beta'],
                          a=damp_params['mbd_rsscs_a'])[:3]

    return energy
