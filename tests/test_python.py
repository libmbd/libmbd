import pytest
from pytest import approx

from pymbd import mbd_energy_species


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
