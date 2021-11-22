import numpy as np
import pytest
from pytest import approx

from pymbd import ang, from_volumes
from pymbd.fortran import MBDFortranError, MBDGeom
from pymbd.utils import numerical_gradients, numerical_latt_gradients


def test_argon_dimer_plain():
    ene = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)]).mbd_energy(
        [11, 11], [63.525, 63.525], [3.55, 3.55], 0.83, variant='plain'
    )
    assert ene == approx(-0.00024329110270970844, rel=1e-10)


@pytest.mark.no_scalapack
def test_argon_dimer_dipole_matrix():
    dip = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)]).dipole_matrix('bare')
    assert (dip != 0).sum() == 6


def test_argon_dimer_rsscs():
    ene = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)]).mbd_energy_species(
        ['Ar', 'Ar'], [1, 1], 0.83
    )
    assert ene == approx(-0.0002462647623815428, rel=1e-10)


def test_argon_dimer_rsscs_rpa():
    geom = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)], do_rpa=True, get_rpa_orders=True)
    ene, orders = geom.mbd_energy_species(['Ar', 'Ar'], [1, 1], 0.83)
    assert ene == approx(-0.0002462647623815428, rel=1e-10)
    assert orders[1] == approx(-0.0002461558113413099, rel=1e-10)
    assert orders[2] == approx(0)
    assert orders[3] == approx(-1.0885208380438466e-07, rel=1e-10)


def test_argon_dimer_ts():
    ene = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)]).ts_energy(
        [11, 11], [63.525, 63.525], [3.55, 3.55], 0.94
    )
    assert ene == approx(-0.000318123017869182, rel=1e-10)


def test_benzene_dimer(benzene_dimer):
    mon1, mon2 = benzene_dimer
    dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
    enes = [
        MBDGeom(coords).mbd_energy_species(species, vol_ratios, 0.83)
        for coords, species, vol_ratios in (mon1, mon2, dim)
    ]
    ene_int = enes[2] - enes[1] - enes[0]
    assert ene_int == approx(-0.006312323931302544, rel=1e-10)


def test_benzene_gradients(benzene_dimer):
    coords, species, vol_ratios = benzene_dimer[0]
    ene, gradients = MBDGeom(coords).mbd_energy_species(
        species, vol_ratios, 0.83, force=True
    )
    with MBDGeom(coords) as geom:
        num_gradients = numerical_gradients(
            geom, 'mbd_energy_species', species, vol_ratios, 0.83
        )
    for i in range(len(coords)):
        assert gradients[i] == approx(num_gradients[i], rel=1e-10, abs=1e-10)


@pytest.mark.no_scalapack
def test_benzene_dimer_python(benzene_dimer):
    mon1, mon2 = benzene_dimer
    dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
    enes = [
        MBDGeom(coords).mbd_energy_species(species, vol_ratios, 0.83)
        for coords, species, vol_ratios in (mon1, mon2, dim)
    ]
    ene_int = enes[2] - enes[1] - enes[0]
    assert ene_int == approx(-0.006312323931302544, rel=1e-10)


def test_benzene_gradients_plain(benzene_dimer):
    coords, species, vol_ratios = benzene_dimer[0]
    ene, gradients = MBDGeom(coords).mbd_energy_species(
        species, vol_ratios, 0.83, variant='plain', force=True
    )
    with MBDGeom(coords) as geom:
        num_gradients = numerical_gradients(
            geom, 'mbd_energy_species', species, vol_ratios, 0.83, variant='plain'
        )
    for i in range(len(coords)):
        assert gradients[i] == approx(num_gradients[i], rel=1e-10, abs=1e-10)


def test_benzene_dimer_scs(benzene_dimer):
    mon1, mon2 = benzene_dimer
    dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
    enes = [
        MBDGeom(coords).mbd_energy_species(
            species, vol_ratios, 1, a=2.56, variant='scs'
        )
        for coords, species, vol_ratios in (mon1, mon2, dim)
    ]
    ene_int = enes[2] - enes[1] - enes[0]
    assert ene_int == approx(-0.007462380657774048, rel=1e-10)


def test_benzene_dimer_ts(benzene_dimer):
    mon1, mon2 = benzene_dimer
    dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
    enes = [
        MBDGeom(coords).ts_energy_species(species, vol_ratios, 0.94)
        for coords, species, vol_ratios in (mon1, mon2, dim)
    ]
    ene_int = enes[2] - enes[1] - enes[0]
    assert ene_int == approx(-0.008490052683234028, rel=1e-10)


def test_benzene(benzene_dimer):
    coords, species, vol_ratios = benzene_dimer[0]
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    ene = MBDGeom(coords).mbd_energy(alpha_0, C6, R_vdw, 0.83, variant='plain')
    assert ene == approx(-0.007002398506090302, rel=1e-10)


def test_benzene_rpa(benzene_dimer):
    coords, species, vol_ratios = benzene_dimer[0]
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    ene = MBDGeom(coords, do_rpa=True).mbd_energy(
        alpha_0, C6, R_vdw, 0.83, variant='plain'
    )
    assert ene == approx(-0.007002398506090302, rel=1e-9)


def test_benzene_rpa_scaled(benzene_dimer):
    coords, species, vol_ratios = benzene_dimer[0]
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    ene = MBDGeom(coords, do_rpa=True, rpa_rescale_eigs=True).mbd_energy(
        alpha_0, C6, R_vdw, 0.83, variant='plain'
    )
    assert ene != approx(-0.007002398506090302, rel=1e-9)
    assert ene == approx(-0.007002398506090302, rel=1e-7)


def test_ethylcarbamate(ethylcarbamate):
    enes = [
        MBDGeom(coords, lattice, k_grid).mbd_energy_species(species, vol_ratios, 0.83)
        for coords, lattice, k_grid, species, vol_ratios in ethylcarbamate
    ]
    ene_int = enes[0] - 2 * enes[1]
    assert ene_int == approx(-0.037040868610822564, rel=1e-10)


def test_argon_crystal(argon_crystal):
    coords, lattice, k_grid, species, vol_ratios = argon_crystal
    ene = MBDGeom(coords, lattice, k_grid).mbd_energy_species(species, vol_ratios, 0.83)
    assert ene == approx(-0.0021037562496878173, rel=1e-10)


def test_argon_crystal_rpa(argon_crystal):
    coords, lattice, k_grid, species, vol_ratios = argon_crystal
    ene = MBDGeom(coords, lattice, k_grid, do_rpa=True).mbd_energy_species(
        species, vol_ratios, 0.83, variant='plain'
    )
    assert ene == approx(-0.0021036969146744147, rel=1e-10)


@pytest.mark.no_scalapack
def test_argon_crystal_modes(argon_crystal):
    coords, lattice, k_grid, species, vol_ratios = argon_crystal
    geom = MBDGeom(
        coords, lattice, custom_k_pts=[(0, 0, 0), (0.1, 0, 0)], get_spectrum=True
    )
    _, eigs, C = geom.mbd_energy_species(species, vol_ratios, 0.83)
    assert abs(C[:, :, 0].imag).sum() == approx(0)
    assert abs(C[:, :, 1].imag).sum() != approx(0)


def test_argon_crystal_gradients(argon_crystal):
    coords, lattice, k_grid, species, vol_ratios = argon_crystal
    ene, gradients, latt_gradients = MBDGeom(
        coords, lattice, k_grid
    ).mbd_energy_species(species, vol_ratios, 0.83, force=True)
    with MBDGeom(coords, lattice, k_grid) as geom:
        num_gradients = numerical_gradients(
            geom, 'mbd_energy_species', species, vol_ratios, 0.83
        )
    for i in range(len(coords)):
        assert gradients[i] == approx(num_gradients[i], rel=1e-10, abs=1e-10)
    with MBDGeom(coords, lattice, k_grid) as geom:
        num_latt_gradients = numerical_latt_gradients(
            geom, 'mbd_energy_species', species, vol_ratios, 0.83
        )
    for i in range(3):
        assert latt_gradients[i] == approx(num_latt_gradients[i], rel=1e-10, abs=1e-10)


def test_lithium(bulk_lithium):
    coords, lattice, k_grid, species, vol_ratios = bulk_lithium
    with pytest.raises(MBDFortranError):
        MBDGeom(coords, lattice, k_grid).mbd_energy_species(species, vol_ratios, 0.83)


def test_ethylcarbamate_scs(ethylcarbamate):
    enes = [
        MBDGeom(coords, lattice, k_grid).mbd_energy_species(
            species, vol_ratios, 1, a=2.56, variant='scs'
        )
        for coords, lattice, k_grid, species, vol_ratios in ethylcarbamate
    ]
    ene_int = enes[0] - 2 * enes[1]
    assert ene_int == approx(-0.03633331132194684, rel=1e-10)


def test_ethylcarbamate_ts(ethylcarbamate):
    enes = [
        MBDGeom(coords, lattice).ts_energy_species(species, vol_ratios, 0.83)
        for coords, lattice, _, species, vol_ratios in ethylcarbamate
    ]
    ene_int = enes[0] - 2 * enes[1]
    assert ene_int == approx(-0.05218213230219945, rel=1e-10)


@pytest.mark.no_scalapack
def test_mbd_coulomb(peptide_meoh):
    a = 14.4
    beta = 2.0
    enes = []
    for coords, species, vol_ratios in peptide_meoh:
        geom = MBDGeom(coords, get_spectrum=True)
        _, eigs, C = geom.mbd_energy_species(species, vol_ratios, beta=0.83)
        omega_t = np.sqrt(eigs)
        alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
        omega = 4 * C6 / (3 * alpha_0 ** 2)
        charges = np.ones_like(alpha_0)
        masses = 1 / (alpha_0 * omega ** 2)
        ecoul = geom.coulomb_energy(
            charges, masses, omega_t, 'fermi', R_vdw, beta, a, C
        )
        edip = geom.dipole_energy(
            alpha_0, omega, omega_t, 'fermi,dip', R_vdw, beta, a, C
        )
        C = np.identity(len(omega_t))
        omega_non = np.repeat(omega, 3)
        ecoul_non = geom.coulomb_energy(
            charges, masses, omega_non, 'fermi', R_vdw, beta, a, C
        )
        edip_non = geom.dipole_energy(
            alpha_0, omega, omega_t, 'fermi,dip', R_vdw, beta, a, C
        )
        enes.append(ecoul - edip - (ecoul_non - edip_non))
    ene_int = enes[2] - enes[0] - enes[1]
    assert ene_int == approx(0.0002460638172163822 / 627.503, rel=1e-10)
