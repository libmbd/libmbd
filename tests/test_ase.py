import pytest
from ase import Atoms
from ase.units import Hartree
from pytest import approx

from pymbd import ang
from pymbd.ase.mbd import MBD


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
        calc = MBD(scheme=scheme, params='TS', k_grid=k_grid, **kwargs)
        atoms.set_calculator(calc)
        atoms.calc.set_hirshfeld(vol_ratios)
        enes.append(atoms.get_potential_energy())
    ene_int = (enes[0] - 2 * enes[1]) / Hartree
    # TODO ASE accuracy for some reason down to 1e-6
    assert ene_int == approx(ene_ref, rel=1e-6)


def test_gradients(benzene_dimer):
    mol1, _ = benzene_dimer
    benzene = Atoms(mol1[1], positions=mol1[0] / ang)
    calc1 = MBD(scheme='MBD', params='TS', beta=0.83, k_grid=None)
    benzene.set_calculator(calc1)
    benzene.calc.set_hirshfeld(mol1[2])
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
    calc = MBD(scheme='MBD', params='TS', beta=0.83, k_grid=k_grid)
    atoms.set_calculator(calc)
    atoms.calc.set_hirshfeld(vol_ratios)
    stress = atoms.get_stress(voigt=False)
    stress_num = atoms.calc.calculate_numerical_stress(atoms, d=1e-3, voigt=False)
    print(stress / stress_num)
    for i in range(len(atoms.positions)):
        assert stress[i] == approx(stress_num[i], rel=1e-7, abs=1e-7)
