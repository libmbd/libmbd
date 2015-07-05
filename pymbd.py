from __future__ import print_function
import json
import numpy as np
from mpi4py import MPI
from mbd import mbd


def init():
    mbd.n_atoms = natoms
    mbd.coords = np.array(data['coords']).T
    mbd.is_periodic = False
    mbd.lattice_vector = np.zeros((3, 3))
    mbd.electric_field = np.zeros(3)
    mbd.n_tasks = ntasks
    mbd.my_task = myid


def main():

    damp_params = dict(zip(['ts_d', 'ts_s_r', 'mbd_scs_a', 'mbd_ts_a',
                            'mbd_ts_erf_beta', 'mbd_ts_fermi_beta', 'mbd_rsscs_a',
                            'mbd_rsscs_beta'],
                           mbd.get_damping_parameters('pbe')))

    nomega = mbd.omega_grid.shape[0]
    alpha, C6, R_vdw, energy = {}, {}, {}, {}

    alpha['TS'] = np.zeros((nomega, natoms))
    alpha['TS'][0] = data['alpha_0']*data['volume_ratio']
    C6['TS'] = data['C6']*np.power(data['volume_ratio'], 2)
    for alphas, omega in zip(alpha['TS'][1:], mbd.omega_grid[1:]):
        alphas[:] = mbd.alpha_dynamic_ts(alpha['TS'][0], C6['TS'], omega)
    R_vdw['TS'] = data['R_vdw']*np.power(data['volume_ratio'], 1./3)

    energy['TS@TS~fermi@TS'] = mbd.get_ts_energy(
        C6['TS'], alpha['TS'][0], 'fermi', R_vdw['TS'],
        d=damp_params['ts_d'], s_r=damp_params['ts_s_r']
    )
    energy['MBD@TS~erf@TS,dip'] = mbd.get_mbd_energy(
        mbd.omega_eff(C6['TS'], alpha['TS'][0]), alpha['TS'][0],
        'erf,dip', R_vdw['TS'],
        beta=damp_params['mbd_ts_erf_beta'], a=4.
    )[0]
    energy['MBD@TS~fermi@TS,dip'] = mbd.get_mbd_energy(
        mbd.omega_eff(C6['TS'], alpha['TS'][0]), alpha['TS'][0],
        'fermi@TS,dip', R_vdw['TS'],
        beta=damp_params['mbd_ts_fermi_beta'], a=damp_params['mbd_ts_a']
    )[0]

    alpha['SCS'] = np.zeros((nomega, natoms))
    for alphas_ts, alphas_scs in zip(alpha['TS'], alpha['SCS']):
        alpha_3n = mbd.run_scs(alphas_ts, 'dip,gg')
        alphas_scs[:] = mbd.contract_polarizability(alpha_3n)
    C6['SCS'] = mbd.get_c6_from_alpha(alpha['SCS'])
    R_vdw['SCS'] = data['R_vdw']*np.power(alpha['SCS'][0]/data['alpha_0'], 1./3)

    energy['TS@SCS~fermi@SCS'] = mbd.get_ts_energy(
        C6['SCS'], alpha['SCS'][0], 'fermi', R_vdw['SCS'],
        d=damp_params['ts_d'], s_r=damp_params['ts_s_r']
    )
    energy['MBD@SCS~dip,1mexp@SCS'] = mbd.get_mbd_energy(
        mbd.omega_eff(C6['SCS'], alpha['SCS'][0]), alpha['SCS'][0],
        'dip,1mexp', R_vdw['SCS'],
        beta=1., a=damp_params['mbd_scs_a']
    )[0]

    alpha['rsSCS'] = np.zeros((nomega, natoms))
    for alphas_ts, alphas_rsscs in zip(alpha['TS'], alpha['rsSCS']):
        alpha_3n = mbd.run_scs(alphas_ts, 'fermi,dip,gg', R_vdw['TS'],
                               beta=damp_params['mbd_rsscs_beta'],
                               a=damp_params['mbd_rsscs_a'])
        alphas_rsscs[:] = mbd.contract_polarizability(alpha_3n)
    C6['rsSCS'] = mbd.get_c6_from_alpha(alpha['rsSCS'])
    R_vdw['rsSCS'] = data['R_vdw']*np.power(alpha['rsSCS'][0]/data['alpha_0'], 1./3)

    energy['MBD@rsSCS~fermi@rsSCS,dip'] = mbd.get_mbd_energy(
        mbd.omega_eff(C6['rsSCS'], alpha['rsSCS'][0]), alpha['rsSCS'][0],
        'fermi@rsSCS,dip', R_vdw['rsSCS'],
        beta=damp_params['mbd_rsscs_beta'], a=damp_params['mbd_rsscs_a']
    )[0]
    energy['MBD(TS)@rsSCS~fermi@rsSCS,dip'] = mbd.get_ts_energy(
        C6['rsSCS'], alpha['rsSCS'][0], 'fermi2', R_vdw['rsSCS'],
        s_r=damp_params['mbd_rsscs_beta'], d=damp_params['mbd_rsscs_a']
    )
    energy['MBD(pair)@rsSCS~fermi@rsSCS,dip'] = mbd.pairwise_mbd(
        mbd.omega_eff(C6['rsSCS'], alpha['rsSCS'][0]), alpha['rsSCS'][0],
        'fermi@rsSCS,dip', R_vdw['rsSCS'],
        beta=damp_params['mbd_rsscs_beta'], a=damp_params['mbd_rsscs_a']
    )[0]
    energy['MBD(RPA)@rsSCS~fermi@rsSCS,dip'] = mbd.get_qho_rpa_energy(
        alpha['rsSCS'], 'fermi@rsSCS,dip', R_vdw['rsSCS'],
        beta=damp_params['mbd_rsscs_beta'], a=damp_params['mbd_rsscs_a']
    )
    energy['MBD(nbody)@rsSCS~fermi@rsSCS,dip'] = mbd.nbody_mbd(
        mbd.omega_eff(C6['rsSCS'], alpha['rsSCS'][0]), alpha['rsSCS'][0],
        'fermi@rsSCS,dip', R_vdw['rsSCS'],
        beta=damp_params['mbd_rsscs_beta'], a=damp_params['mbd_rsscs_a']
    )

    return energy


data = json.load(open('results.json'))
for key in data:
    data[key] = np.array(data[key])

bohr = 0.529177249
natoms = len(data['coords'])
ntasks = MPI.COMM_WORLD.Get_size()
myid = MPI.COMM_WORLD.Get_rank()

if myid == 0:
    print('Running on %s nodes...' % ntasks)

init()
energy = main()

if myid == 0:
    for key in energy:
        my = energy[key]
        ref = [e['value'] for e in data['energy'] if e['name'] == key][0]
        assert np.linalg.norm(my-ref) < 1e-10
    print('Success: All energies match')
