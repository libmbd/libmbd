import numpy as np
import pytest
from pytest import approx

from pymbd import (
    atomic_polarizabilities,
    from_volumes,
    mbd_energy_species,
    molecular_polarizability,
    screening,
    screening_matrix,
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
def test_screening_matrix(benzene_dimer):
    coords, species, vol_ratios = benzene_dimer[0]
    n = len(species)
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    a_nlc = screening_matrix(coords, alpha_0, R_vdw, 0.83)
    assert a_nlc.shape == (3 * n, 3 * n)
    # inverse of a symmetric matrix, so symmetric
    assert a_nlc == approx(a_nlc.T)


@pytest.mark.no_scalapack
def test_polarizabilities(benzene_dimer):
    coords, species, vol_ratios = benzene_dimer[0]
    na = len(species)
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    alpha_scs = screening(coords, alpha_0, C6, R_vdw, 0.83)[0]
    # atomic polarizabilities are per-atom isotropic scalars, equal to the
    # static screened polarizabilities returned by screening()
    alpha_a = atomic_polarizabilities(coords, species, vol_ratios, 0.83)
    assert alpha_a.shape == (na,)
    assert alpha_a == approx(alpha_scs)
    # the molecular polarizability is a symmetric 3x3 tensor whose isotropic
    # part is the sum of the atomic polarizabilities
    alpha_mol = molecular_polarizability(coords, species, vol_ratios, 0.83)
    assert alpha_mol.shape == (3, 3)
    assert alpha_mol == approx(alpha_mol.T)
    assert np.trace(alpha_mol) / 3 == approx(alpha_a.sum())
