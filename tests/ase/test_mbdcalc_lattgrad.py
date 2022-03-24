import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree
from pytest import approx

from pymbd.ase.mbd import MBD
from pymbd.fortran import MBDGeom
from pymbd.utils import numerical_latt_gradients

ang2 = 1.0 / 0.529177249
ang = 1.0  # /0.529177249


def argon_crystal():
    return (
        np.array([(0.3, 0.1, 0.2), (4.1, -0.2, -0.1)]) * ang,
        np.array([(8.1, 0.1, -0.2), (0.3, 3.9, -0.1), (-0.1, 0.2, 4.2)]) * ang,
        (4, 4, 4),
        ['Ar', 'Ar'],
        [1.0, 1.0],
    )


print('TEST ASE libmbd calculator')

coords, lattice, k_grid, species, vol_ratios = argon_crystal()
atoms = Atoms(species, positions=coords, cell=lattice, pbc=[True, True, True])

calc1 = MBD(scheme='MBD', params='TS', beta=0.83, k_grid=k_grid)

np.set_printoptions(precision=6, suppress=True)

atoms.set_calculator(calc1)
atoms.calc.set_hirshfeld(vol_ratios)
print(atoms.get_pbc())
stress = atoms.get_stress(voigt=False)
print('stress')
print(stress)

stress_num = atoms.calc.calculate_numerical_stress(atoms, d=1e-3, voigt=False)
print('stress num')
print(stress_num)

print(stress / stress_num)

# for i in range(len(atoms.positions)):
#    assert forces[i] == approx(forces_num[i], rel=1e-7, abs=1e-7)

# for i in range(len(atoms.positions)):
#    assert stress[i] == approx(stress_num[i], rel=1e-7, abs=1e-7)


print('TEST libmbd')

coords, lattice, k_grid, species, vol_ratios = argon_crystal()
ene, gradients, latt_gradients = MBDGeom(
    coords * ang2, lattice * ang2, k_grid
).mbd_energy_species(species, vol_ratios, 0.83, force=True)

# with MBDGeom(coords, lattice, k_grid) as geom:
#    num_gradients = numerical_gradients(
#       geom, 'mbd_energy_species', species, vol_ratios, 0.83
#    )
# for i in range(len(coords)):
#    assert gradients[i] == approx(num_gradients[i], rel=1e-10, abs=1e-10)

with MBDGeom(coords * ang2, lattice * ang2, k_grid) as geom:
    num_latt_gradients = numerical_latt_gradients(
        geom, 'mbd_energy_species', species, vol_ratios, 0.83
    )
print('libmbd')
print(latt_gradients * Hartree / Bohr)
print('num libmbd')
print(num_latt_gradients * Hartree / Bohr)
for i in range(3):
    assert latt_gradients[i] == approx(num_latt_gradients[i], rel=1e-10, abs=1e-10)
