import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree
from pytest import approx

from pymbd import ang
from pymbd.ase import MBD, DispersionCorrectionCalculator


class DummyCalculator(Calculator):
    implemented_properties = ['energy', 'forces', 'stress', 'free_energy']

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {
            'energy': 0.0,
            'free_energy': 0.0,
            'forces': np.zeros((len(atoms), 3)),
            'stress': np.zeros(6),
        }


def get_mbd(*, hirbulk, **kwargs):
    return DispersionCorrectionCalculator(
        qm_calculator=DummyCalculator(),
        mm_calculator=MBD(**kwargs),
        hirbulk=hirbulk,
    )


@pytest.mark.parametrize(
    'scheme,ene_ref,kwargs',
    [
        ('MBD', -0.037040868610822564, {'beta': 0.83}),
        ('VDW', -0.05218213230219945, {'ts_sr': 0.83}),
    ],
    ids=['MBD', 'TS'],
)
@pytest.mark.no_scalapack
def test_scheme(ethylcarbamate, scheme, ene_ref, kwargs):
    enes = []
    for coords, lattice, k_grid, species, vol_ratios in ethylcarbamate:
        atoms = Atoms(
            species,
            positions=coords / ang,
            cell=lattice / ang if lattice is not None else None,
            pbc=(True, True, True) if lattice is not None else None,
        )
        calc = get_mbd(
            scheme=scheme, params='TS', k_grid=k_grid, hirbulk=vol_ratios, **kwargs
        )
        atoms.set_calculator(calc)
        enes.append(atoms.get_potential_energy())
    ene_int = (enes[0] - 2 * enes[1]) / Hartree
    # TODO ASE accuracy for some reason down to 1e-6
    assert ene_int == approx(ene_ref, rel=1e-6)


def test_gradients(benzene_dimer):
    mol1, _ = benzene_dimer
    benzene = Atoms(mol1[1], positions=mol1[0] / ang)
    calc1 = get_mbd(scheme='MBD', params='TS', beta=0.83, k_grid=None, hirbulk=mol1[2])
    benzene.set_calculator(calc1)
    forces = benzene.get_forces()
    forces_num = benzene.calc.calculate_numerical_forces(benzene)
    # TODO ASE can only achieve a force accuracy wrt to numerical by 1e-7
    for i in range(len(benzene.positions)):
        assert forces[i] == approx(forces_num[i], rel=1e-7, abs=1e-7)


def test_lattgrad(argon_crystal):
    coords, lattice, k_grid, species, vol_ratios = argon_crystal
    atoms = Atoms(
        species, positions=coords / ang, cell=lattice / ang, pbc=[True, True, True]
    )
    calc = get_mbd(
        scheme='MBD', params='TS', beta=0.83, k_grid=k_grid, hirbulk=vol_ratios
    )
    atoms.set_calculator(calc)
    stress = atoms.get_stress(voigt=False)
    stress_num = atoms.calc.calculate_numerical_stress(atoms, d=1e-3, voigt=False)
    print(stress / stress_num)
    for i in range(len(atoms.positions)):
        assert stress[i] == approx(stress_num[i], rel=1e-7, abs=1e-7)
