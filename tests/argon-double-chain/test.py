import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent/'../..'))
from pymbd import get_free_atom_data, printerr, myid, ntasks
from mbd import mbd
import numpy as np

mbd.my_task = myid
mbd.n_tasks = ntasks
mode = 'M' if ntasks > 1 else ''

bohr = mbd.bohr
mbd.param_dipole_matrix_accuracy = 1e-10
mbd.init_grid(30)

species = ['Ar', 'Ar']
xyz = [(0., 0., 0.), (0., 0., 4.)]/bohr
uc = np.array([(4., 0., 0.), (0., 10., 0.), (0., 0., 10.)])/bohr
mbd.param_vacuum_axis = (False, True, True)
alpha_0, C6, R_vdw = get_free_atom_data(species)
omega = mbd.omega_eff(C6, alpha_0)

# mbd.measure_time = True

ns_kpt = [4, 8, 12, 20, 40, 80]
enes_periodic = []
for n_kpt in ns_kpt:
    k_grid = mbd.make_k_grid(mbd.make_g_grid(n_kpt, 1, 1), uc)
    ene = mbd.get_reciprocal_mbd_energy('R' + mode, 'fermi,dip', xyz, alpha_0, omega, k_grid, uc,
                                        r_vdw=R_vdw, beta=1., a=6.)[0]
    enes_periodic.append(ene)
printerr(list(reversed(enes_periodic)))
# printout(mbd.timestamps/mbd.clock_rate())
# printout(mbd.ts_counts)

cutoffs = [25., 50., 100., 200., 400., 800.]
enes_supercell = []
for cutoff in cutoffs:
    mbd.param_mbd_supercell_cutoff = cutoff/bohr
    ene = mbd.get_supercell_mbd_energy(
        'P' + mode, 'fermi,dip', xyz, alpha_0, omega, uc,
        r_vdw=R_vdw, beta=1., a=6.)[0]
    enes_supercell.append(ene)
printerr(list(reversed(enes_supercell)))
# printerr(mbd.timestamps/mbd.clock_rate())
# printerr(mbd.ts_counts)

for _ in range(1000):
    ene = mbd.get_ts_energy('P' + mode, 'fermi2', xyz, C6, alpha_0,
                            r_vdw=R_vdw, s_r=1., d=6., unit_cell=uc)
printerr(ene)
# printout(mbd.timestamps/mbd.clock_rate())
# printout(mbd.ts_counts//1000)
