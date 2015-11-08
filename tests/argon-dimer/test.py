import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent/'../..'))
from pymbd import get_free_atom_data, printerr, myid, ntasks
from mbd import mbd


mbd.my_task = myid
mbd.n_tasks = ntasks
mode = 'M' if ntasks > 1 else ''

bohr = mbd.bohr
mbd.param_dipole_matrix_accuracy = 1e-10
mbd.init_grid(30)

species = ['Ar', 'Ar']
xyz = [(0., 0., 0.), (4., 0., 0.)]/bohr
alpha_0, C6, R_vdw = get_free_atom_data(species)
omega = mbd.omega_eff(C6, alpha_0)

ene = mbd.get_single_mbd_energy(mode, 'fermi,dip', xyz, alpha_0, omega,
                                r_vdw=R_vdw, beta=1., a=6.)[0]
printerr(ene)

ene = mbd.get_qho_rpa_energy(mode, 'fermi,dip', xyz,
                             mbd.alpha_dynamic_ts_all(alpha_0, mbd.n_grid_omega, c6=C6),
                             r_vdw=R_vdw, beta=1., a=6.)[0]
printerr(ene)

ene = mbd.get_ts_energy(mode, 'fermi2', xyz, C6, alpha_0,
                        r_vdw=R_vdw, s_r=1., d=6.)
printerr(ene)
