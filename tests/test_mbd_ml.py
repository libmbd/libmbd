from importlib.util import find_spec

import numpy as np
import pytest
from pytest import approx

pytest.importorskip('ase')

from ase import Atoms  # noqa: E402
from ase.build import molecule  # noqa: E402

from pymbd.mbd_ml import (  # noqa: E402
    compute_stress_from_lattice_gradient,
    mbd_properties_from_structure,
    ratios_from_mbdml,
    resolve_model,
)

BOHR = 0.5291772105638411

# so3lr is not on PyPI and pulls a large git-dependency chain, so the tests that
# need the model are opt-in; everything above the model call runs without it.
needs_so3lr = pytest.mark.skipif(
    find_spec('so3lr') is None,
    reason='so3lr is not installed (see requirements-mbd-ml.txt)',
)


# --- pieces that need no model ------------------------------------------------


@pytest.mark.no_scalapack
def test_stress_from_lattice_gradient_matches_finite_differences():
    # The stress is the derivative of the energy with respect to a strain eps
    # applied as L -> L(1 + eps), R -> R(1 + eps), divided by the cell volume.
    # Check the closed form against a numerical strain derivative of an
    # arbitrary smooth E(L, R), which exercises the formula without needing an
    # actual energy model.
    rng = np.random.default_rng(0)
    lattice = np.eye(3) * 5.0 + rng.normal(scale=0.3, size=(3, 3))
    coords = rng.normal(scale=2.0, size=(4, 3))

    def energy(latt, crd):
        return np.sin(latt).sum() + (crd**2).sum()

    dE_dlattice = np.cos(lattice)
    dE_dcoords = 2 * coords

    ours = compute_stress_from_lattice_gradient(
        lattice, coords, dE_dlattice, dE_dcoords
    )

    volume = abs(np.linalg.det(lattice))
    h = 1e-6
    numerical = np.zeros((3, 3))
    for k in range(3):
        for m in range(3):
            plus, minus = np.zeros((3, 3)), np.zeros((3, 3))
            plus[k, m], minus[k, m] = h, -h
            e_plus = energy(lattice @ (np.eye(3) + plus), coords @ (np.eye(3) + plus))
            e_minus = energy(
                lattice @ (np.eye(3) + minus), coords @ (np.eye(3) + minus)
            )
            numerical[k, m] = (e_plus - e_minus) / (2 * h) / volume

    assert ours == approx(numerical, abs=1e-7)


@pytest.mark.no_scalapack
def test_periodic_structure_without_k_grid_is_rejected():
    # The guard runs before the model is consulted, so this holds whether or not
    # so3lr is installed.
    atoms = Atoms('Ar', positions=[(0, 0, 0)], cell=np.eye(3) * 5, pbc=True)
    with pytest.raises(ValueError, match='k_grid'):
        mbd_properties_from_structure(atoms, beta=0.83)


@pytest.mark.no_scalapack
def test_resolve_model_accepts_a_directory(tmp_path):
    (tmp_path / 'hyperparameters.json').write_text('{}')
    assert resolve_model(str(tmp_path)) == tmp_path.resolve()


@pytest.mark.no_scalapack
def test_mixed_periodicity_is_rejected():
    # A slab is neither a molecule nor a 3D crystal here: the lattice handed to
    # libMBD would be singular, so it must not be silently treated as periodic.
    atoms = Atoms(
        'Ar', positions=[(0, 0, 0)], cell=np.eye(3) * 5, pbc=(True, True, False)
    )
    with pytest.raises(ValueError, match='fully periodic'):
        mbd_properties_from_structure(atoms, beta=0.83, k_grid=(2, 2, 1))


# --- pieces that need the model -----------------------------------------------


@needs_so3lr
@pytest.mark.no_scalapack
def test_resolve_model_rejects_something_that_is_neither(tmp_path):
    with pytest.raises(FileNotFoundError, match='neither'):
        resolve_model(str(tmp_path / 'nope'))


@needs_so3lr
@pytest.mark.no_scalapack
def test_water_ratios():
    # Regression against the shipped model. A change here means the checkpoint,
    # the loader, or the model resolution changed -- all of which have moved
    # under this module before.
    ratios = ratios_from_mbdml(molecule('H2O'))
    assert ratios['a0'] == approx([1.007767, 0.318508, 0.318508], abs=1e-4)
    assert ratios['c6'] == approx([1.082644, 0.157658, 0.157658], abs=1e-4)
    # hydrogens contract in the O-H bonds, oxygen stays near-free
    assert ratios['a0'][0] > ratios['a0'][1]


@needs_so3lr
@pytest.mark.no_scalapack
def test_forces_are_the_gradient_of_the_energy():
    # alpha_0, C6 and R_vdw come from a model that depends on the geometry, so
    # libMBD's fixed-parameter gradient is not the gradient of the energy. With
    # the response term the two agree; without it the error is two orders of
    # magnitude larger, which is what this pins down.
    atoms = molecule('H2O')
    forces = np.asarray(mbd_properties_from_structure(atoms, beta=0.83)['F'])
    bare = np.asarray(
        mbd_properties_from_structure(atoms, beta=0.83, ratio_response=False)['F']
    )

    h = 1e-3  # Angstrom
    numerical = np.zeros(3)
    for k in range(3):
        shifted = []
        for sign in (+1, -1):
            moved = atoms.copy()
            positions = moved.get_positions()
            positions[0, k] += sign * h
            moved.set_positions(positions)
            shifted.append(
                mbd_properties_from_structure(moved, beta=0.83, ratio_response=False)[
                    'E'
                ]
            )
        numerical[k] = -(shifted[0] - shifted[1]) / (2 * h / BOHR)

    assert forces[0] == approx(numerical, abs=1e-5)
    # and the term is worth having: without it the force is off by ~40%
    assert np.abs(bare[0] - numerical).max() > 10 * np.abs(forces[0] - numerical).max()


@needs_so3lr
@pytest.mark.no_scalapack
def test_ratio_response_leaves_the_energy_alone():
    # The response term is a gradient correction only; the energy is the same.
    atoms = molecule('H2O')
    with_response = mbd_properties_from_structure(atoms, beta=0.83)
    without = mbd_properties_from_structure(atoms, beta=0.83, ratio_response=False)
    assert with_response['E'] == approx(without['E'], rel=1e-12)
