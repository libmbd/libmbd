import numpy as np
import pytest
from pytest import approx

from pymbd import from_volumes
from pymbd.fortran import MBDGeom
from pymbd.utils import numerical_gradients, numerical_latt_gradients


@pytest.mark.no_scalapack
@pytest.mark.parametrize('k_point', [None, (0, 0, 0), (0.1, 0.07, -0.03)])
@pytest.mark.parametrize('damping', ['bare', 'fermi,dip,gg'])
def test_dipole_matrix_unwrapped(argon_crystal, k_point, damping):
    coords, lattice, _, species, volumes = argon_crystal
    _, _, radii = from_volumes(species, volumes)
    with MBDGeom(coords, lattice) as geom:
        reference = geom.dipole_matrix(
            damping, beta=0.83, R_vdw=radii, sigma=[1, 1], k_point=k_point
        )
        # different shifts in all three directions of a nonorthogonal cell
        geom.coords = coords + np.array([[19, -23, 17], [-29, 11, -13]]) @ lattice
        shifted = geom.dipole_matrix(
            damping, beta=0.83, R_vdw=radii, sigma=[1, 1], k_point=k_point
        )
    assert shifted == approx(reference, abs=1e-11)


@pytest.mark.parametrize('method', ['ts', 'plain', 'rsscs'])
def test_periodic_energy_and_gradients_unwrapped(argon_crystal, method):
    coords, lattice, k_grid, species, volumes = argon_crystal
    function = 'ts_energy_species' if method == 'ts' else 'mbd_energy_species'
    beta = 0.94 if method == 'ts' else 0.83
    kwargs = {} if method == 'ts' else {'variant': method}
    shifts = np.array([[19, -23, 17], [-29, 11, -13]])
    with MBDGeom(coords, lattice, k_grid) as geom:
        energy, gradient, lattice_gradient = getattr(geom, function)(
            species, volumes, beta, force=True, **kwargs
        )
        geom.coords = coords + shifts @ lattice
        shifted_energy, shifted_gradient, shifted_lattice_gradient = getattr(
            geom, function
        )(species, volumes, beta, force=True, **kwargs)
        assert shifted_energy == approx(energy, abs=1e-12)
        assert shifted_gradient == approx(gradient, abs=1e-11)
        # changing representatives adds a coordinate chain-rule term to dE/dA
        assert shifted_lattice_gradient + shifts.T @ shifted_gradient == approx(
            lattice_gradient, abs=1e-10
        )
        numeric = numerical_gradients(
            geom, function, species, volumes, beta, delta=1e-4, **kwargs
        )
        assert shifted_gradient == approx(numeric, abs=1e-9)
        numeric_lattice = numerical_latt_gradients(
            geom, function, species, volumes, beta, delta=1e-4, **kwargs
        )
        assert shifted_lattice_gradient == approx(numeric_lattice, abs=1e-9)
