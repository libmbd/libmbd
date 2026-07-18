import numpy as np
import pytest
from pytest import approx

from pymbd import (
    atomic_polarizabilities,
    from_volumes,
    mbd_energy_species,
    molecular_polarizability,
    screening,
)


@pytest.mark.no_scalapack
def test_ethylcarbamate(ethylcarbamate):
    enes = [
        mbd_energy_species(
            coords, species, vol_ratios, 0.83, lattice=lattice, k_grid=k_grid
        )
        for coords, lattice, k_grid, species, vol_ratios in ethylcarbamate
    ]
    ene_int = enes[0] - 2 * enes[1]
    assert ene_int == approx(-0.037040868610822564, rel=1e-10)


@pytest.mark.no_scalapack
def test_polarizability_tensors(benzene_dimer):
    coords, species, vol_ratios = benzene_dimer[0]
    na = len(species)
    alpha_a = atomic_polarizabilities(coords, species, vol_ratios, 0.83)
    alpha_mol = molecular_polarizability(coords, species, vol_ratios, 0.83)
    assert alpha_a.shape == (na, 3, 3)
    assert alpha_mol.shape == (3, 3)
    # the molecular tensor is the sum of the atomic tensors
    assert alpha_mol == approx(alpha_a.sum(axis=0))
    # the molecular polarizability tensor is symmetric
    assert alpha_mol == approx(alpha_mol.T)
    # its isotropic part matches the screened static polarizabilities
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    alpha_scs = screening(coords, alpha_0, C6, R_vdw, 0.83)[0]
    assert np.trace(alpha_mol) / 3 == approx(alpha_scs.sum())
