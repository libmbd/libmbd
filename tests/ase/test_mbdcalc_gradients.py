import numpy as np
from ase import Atoms
from ase.calculators.test import gradient_test
from ase.units import Bohr, Hartree
from pytest import approx

from pymbd.ase.mbd import MBD
from pymbd.fortran import MBDGeom
from pymbd.utils import numerical_gradients

ang2 = 1.0 / 0.529177249
ang = 1  # ./0.529177249


def benzene_dimer():
    return [
        (
            np.array(
                [
                    (-1.047, -1.421, 0.000),
                    (-1.454, -0.855, 1.206),
                    (-1.454, -0.855, -1.206),
                    (-2.266, 0.277, 1.206),
                    (-2.671, 0.845, 0.000),
                    (-2.266, 0.277, -1.206),
                    (-1.133, -1.292, -2.142),
                    (-2.582, 0.716, -2.143),
                    (-3.303, 1.723, 0.000),
                    (-2.582, 0.716, 2.143),
                    (-1.133, -1.292, 2.142),
                    (-0.406, -2.291, 0.000),
                ]
            )
            * ang,
            6 * ['C'] + 6 * ['H'],
            [
                0.825,
                0.821,
                0.821,
                0.815,
                0.814,
                0.815,
                0.624,
                0.611,
                0.610,
                0.611,
                0.624,
                0.643,
            ],
        ),
        (
            np.array(
                [
                    (1.047, 1.421, 0.000),
                    (1.454, 0.855, -1.206),
                    (1.454, 0.855, 1.206),
                    (2.266, -0.277, -1.206),
                    (2.671, -0.845, 0.000),
                    (2.266, -0.277, 1.206),
                    (0.406, 2.291, 0.000),
                    (1.133, 1.292, 2.142),
                    (2.582, -0.716, 2.143),
                    (3.303, -1.723, 0.000),
                    (2.582, -0.716, -2.143),
                    (1.133, 1.292, -2.142),
                ]
            )
            * ang,
            6 * ['C'] + 6 * ['H'],
            [
                0.825,
                0.821,
                0.821,
                0.815,
                0.814,
                0.815,
                0.643,
                0.624,
                0.611,
                0.610,
                0.611,
                0.624,
            ],
        ),
    ]


print('TEST ASE libmbd calculator')

mol1, _ = benzene_dimer()
benzene = Atoms(mol1[1], positions=mol1[0])

calc1 = MBD(scheme='MBD', params='TS', beta=0.83, k_grid=None)

benzene.set_calculator(calc1)
benzene.calc.set_hirshfeld(mol1[2])

ene1 = benzene.get_potential_energy()
forces = benzene.get_forces()
print(ene1)
print(forces)

forces_num = benzene.calc.calculate_numerical_forces(benzene)
print(forces_num)

gradient_test(benzene)

# TODO ASE can only achieve a force accuracy wrt to numerical by 1e-7
for i in range(len(benzene.positions)):
    assert forces[i] == approx(forces_num[i], rel=1e-7, abs=1e-7)

print('TEST libmbd')

coords, species, vol_ratios = benzene_dimer()[0]
ene, gradients = MBDGeom(coords * ang2).mbd_energy_species(
    species, vol_ratios, 0.83, force=True
)

print(ene * Hartree)
print(gradients * Hartree / Bohr)

with MBDGeom(coords * ang2) as geom:
    num_gradients = numerical_gradients(
        geom, 'mbd_energy_species', species, vol_ratios, 0.83
    )

for i in range(len(coords)):
    assert gradients[i] == approx(num_gradients[i], rel=1e-10, abs=1e-10)
