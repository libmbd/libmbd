import numpy as np
from ase import Atoms
from ase.units import Hartree
from pytest import approx

from pymbd import from_volumes
from pymbd.ase.mbd import MBD
from pymbd.fortran import MBDGeom

ang2 = 1.0 / 0.529177249
ang = 1  # ./0.529177249


def ethylcarbamate():
    return [
        (
            np.array(
                [
                    (4.083, 5.700, 2.856),
                    (0.568, 0.095, 4.217),
                    (0.470, 4.774, 3.551),
                    (4.181, 1.022, 3.522),
                    (5.572, 5.587, 1.892),
                    (-0.920, 0.209, 5.181),
                    (3.663, 3.255, 2.585),
                    (0.988, 2.541, 4.488),
                    (3.834, 4.011, 0.979),
                    (0.816, 1.785, 6.094),
                    (2.223, 1.314, 1.108),
                    (2.428, 4.481, 5.965),
                    (1.177, 0.092, 0.406),
                    (3.474, 5.703, 6.667),
                    (4.911, 5.036, 2.573),
                    (-0.260, 0.759, 4.500),
                    (4.358, 3.787, 1.918),
                    (0.293, 2.009, 5.155),
                    (0.205, 1.729, 1.101),
                    (4.446, 4.067, 5.972),
                    (1.285, 0.947, 0.957),
                    (3.366, 4.848, 6.116),
                    (0.485, 2.901, 1.709),
                    (4.165, 2.895, 5.364),
                    (4.066, 1.426, 0.670),
                    (0.585, 4.369, 6.403),
                ]
            )
            * ang,
            np.array(
                [(5.008, 0.018, -0.070), (1.630, 6.759, 0.064), (-1.987, -0.981, 7.079)]
            )
            * ang,
            (2, 2, 2),
            list('HHHHHHHHHHHHHHCCCCCCNNOOOO'),
            [
                0.703,
                0.703,
                0.726,
                0.726,
                0.731,
                0.731,
                0.727,
                0.727,
                0.754,
                0.754,
                0.750,
                0.750,
                0.755,
                0.755,
                0.809,
                0.809,
                0.827,
                0.827,
                0.834,
                0.834,
                0.840,
                0.840,
                0.886,
                0.886,
                0.892,
                0.892,
            ],
        ),
        (
            np.array(
                [
                    (4.088, 5.753, 2.783),
                    (5.625, 5.562, 1.906),
                    (3.652, 3.273, 2.592),
                    (3.854, 3.998, 0.981),
                    (5.422, 4.834, 3.521),
                    (6.213, 0.125, 0.386),
                    (7.201, 1.360, 1.112),
                    (4.913, 5.058, 2.573),
                    (4.366, 3.792, 1.934),
                    (5.167, 1.729, 1.084),
                    (6.291, 0.963, 0.938),
                    (4.042, 1.399, 0.752),
                    (5.490, 2.915, 1.682),
                ]
            )
            * ang,
            None,
            None,
            list('HHHHHHHCCCNOO'),
            [
                0.581,
                0.607,
                0.642,
                0.646,
                0.607,
                0.596,
                0.597,
                0.762,
                0.799,
                0.845,
                0.824,
                0.974,
                0.896,
            ],
        ),
    ]


print('TEST1')

mol1, mol2 = ethylcarbamate()
ec_dimer = Atoms(mol1[3], positions=mol1[0], cell=mol1[1], pbc=(True, True, True))
ec = Atoms(mol2[3], positions=mol2[0], cell=mol2[1])

calc1 = MBD(scheme='VDW', params='TS', ts_sr=0.83, k_grid=(2, 2, 2))
calc2 = MBD(scheme='VDW', params='TS', ts_sr=0.83, k_grid=None)

ec_dimer.set_calculator(calc1)
ec.set_calculator(calc2)
ec_dimer.calc.set_hirshfeld(mol1[4])
ec.calc.set_hirshfeld(mol2[4])

ene1 = ec_dimer.get_potential_energy()
print(ene1)
ene2 = ec.get_potential_energy()
print(ene2 * 2)
enes = [ene1, ene2]

print('-----------')

enes_old = []
for coords, lattice, k_grid, species, vol_ratios in ethylcarbamate():
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios, kind='TS')
    # enes_old.append(mbd_energy(
    #        coords, alpha_0, C6, R_vdw, 0.83, lattice=lattice, k_grid=k_grid
    #    ))
    if lattice is None:
        pass
    else:
        lattice *= ang2
    enes_old.append(
        MBDGeom(coords * ang2, lattice, k_grid).ts_energy(alpha_0, C6, R_vdw, 0.83)
    )

ene_int = enes_old[0] - 2 * enes_old[1]
# print('energy in eV ', ene_int)
print(ene_int)

ene_int_ase = (enes[0] - 2 * enes[1]) / Hartree
print(ene_int_ase)

assert ene_int == approx(-0.05218213230219945, rel=1e-10)
# TODO ASE accuracy for some reason only down to 1e-6
assert ene_int_ase == approx(-0.05218213230219945, rel=1e-6)
