# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
import numpy as np
import pytest
from pytest import approx

from pymbd import ang, from_volumes, mbd_energy_species
from pymbd.fortran import MBDFortranException, MBDGeom, with_scalapack
from pymbd.utils import numerical_gradients, numerical_latt_gradients

no_scalapack = pytest.mark.skipif(with_scalapack, reason="doesn't support ScaLAPACK")

benzene_dimer = [
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

ethylcarbamate = [
    (
        np.array(
            [
                (4.083, 5.700, 2.856),
                (0.568, 0.095, 4.217),
                (0.470, 4.774, 3.551),
                (4.181, 1.022, 3.522),
                (5.572, 5.587, 1.892),
                (-0.920, 0.209, 5.181),
                (3.663, 3.255, 2.585),
                (0.988, 2.541, 4.488),
                (3.834, 4.011, 0.979),
                (0.816, 1.785, 6.094),
                (2.223, 1.314, 1.108),
                (2.428, 4.481, 5.965),
                (1.177, 0.092, 0.406),
                (3.474, 5.703, 6.667),
                (4.911, 5.036, 2.573),
                (-0.260, 0.759, 4.500),
                (4.358, 3.787, 1.918),
                (0.293, 2.009, 5.155),
                (0.205, 1.729, 1.101),
                (4.446, 4.067, 5.972),
                (1.285, 0.947, 0.957),
                (3.366, 4.848, 6.116),
                (0.485, 2.901, 1.709),
                (4.165, 2.895, 5.364),
                (4.066, 1.426, 0.670),
                (0.585, 4.369, 6.403),
            ]
        )
        * ang,
        np.array(
            [(5.008, 0.018, -0.070), (1.630, 6.759, 0.064), (-1.987, -0.981, 7.079)]
        )
        * ang,
        (2, 2, 2),
        list('HHHHHHHHHHHHHHCCCCCCNNOOOO'),
        [
            0.703,
            0.703,
            0.726,
            0.726,
            0.731,
            0.731,
            0.727,
            0.727,
            0.754,
            0.754,
            0.750,
            0.750,
            0.755,
            0.755,
            0.809,
            0.809,
            0.827,
            0.827,
            0.834,
            0.834,
            0.840,
            0.840,
            0.886,
            0.886,
            0.892,
            0.892,
        ],
    ),
    (
        np.array(
            [
                (4.088, 5.753, 2.783),
                (5.625, 5.562, 1.906),
                (3.652, 3.273, 2.592),
                (3.854, 3.998, 0.981),
                (5.422, 4.834, 3.521),
                (6.213, 0.125, 0.386),
                (7.201, 1.360, 1.112),
                (4.913, 5.058, 2.573),
                (4.366, 3.792, 1.934),
                (5.167, 1.729, 1.084),
                (6.291, 0.963, 0.938),
                (4.042, 1.399, 0.752),
                (5.490, 2.915, 1.682),
            ]
        )
        * ang,
        None,
        None,
        list('HHHHHHHCCCNOO'),
        [
            0.581,
            0.607,
            0.642,
            0.646,
            0.607,
            0.596,
            0.597,
            0.762,
            0.799,
            0.845,
            0.824,
            0.974,
            0.896,
        ],
    ),
]

argon_crystal = (
    np.array([(0.3, 0.1, 0.2), (4.1, -0.2, -0.1)]) * ang,
    np.array([(8.1, 0.1, -0.2), (0.3, 3.9, -0.1), (-0.1, 0.2, 4.2)]) * ang,
    (4, 4, 4),
    ['Ar', 'Ar'],
    [1.0, 1.0],
)

peptide_meoh = [
    (
        np.array(
            [
                (2.137, 0.252, 0.453),
                (2.857, 0.879, 0.544),
                (2.656, -1.053, 0.687),
                (1.823, -1.742, 0.582),
                (3.422, -1.322, -0.039),
                (3.064, -1.154, 1.693),
            ]
        )
        * ang,
        list('OHCHHH'),
        [0.9114, 0.5960, 0.7523, 0.5886, 0.5850, 0.5850],
    ),
    (
        np.array(
            [
                (-0.849, -0.339, 2.491),
                (0.184, -0.011, 2.416),
                (-0.882, -1.342, 2.912),
                (-1.390, 0.316, 3.168),
                (-1.564, -0.353, 1.159),
                (-2.749, -0.651, 1.056),
                (-0.801, -0.027, 0.088),
                (0.161, 0.240, 0.218),
                (-1.385, -0.002, -1.234),
                (-1.891, -0.942, -1.440),
                (-2.119, 0.796, -1.330),
                (-0.594, 0.149, -1.963),
            ]
        )
        * ang,
        list('CHHHCONHCHHH'),
        [
            0.7657,
            0.6027,
            0.6062,
            0.6077,
            0.8343,
            0.9815,
            0.8325,
            0.5931,
            0.7592,
            0.6286,
            0.6133,
            0.5698,
        ],
    ),
    (
        np.array(
            [
                (-0.849, -0.339, 2.491),
                (0.184, -0.011, 2.416),
                (-0.882, -1.342, 2.912),
                (-1.390, 0.316, 3.168),
                (-1.564, -0.353, 1.159),
                (-2.749, -0.651, 1.056),
                (-0.801, -0.027, 0.088),
                (0.161, 0.240, 0.218),
                (-1.385, -0.002, -1.234),
                (-1.891, -0.942, -1.440),
                (-2.119, 0.796, -1.330),
                (-0.594, 0.149, -1.963),
                (2.137, 0.252, 0.453),
                (2.857, 0.879, 0.544),
                (2.656, -1.053, 0.687),
                (1.823, -1.742, 0.582),
                (3.422, -1.322, -0.039),
                (3.064, -1.154, 1.693),
            ]
        )
        * ang,
        list('CHHHCONHCHHHOHCHHH'),
        [
            0.7767,
            0.6594,
            0.6193,
            0.6167,
            0.8414,
            0.9898,
            0.8462,
            0.7213,
            0.7668,
            0.6367,
            0.6211,
            0.5915,
            0.8615,
            0.5511,
            0.7415,
            0.6022,
            0.5701,
            0.5759,
        ],
    ),
]

bulk_lithium = [
    (
        np.array([(0.0, 0.0, 0.0)]) * ang,
        np.array(
            [
                (-1.7385, 1.7385, 1.7385),
                (1.7385, -1.7385, 1.7385),
                (1.7385, 1.7385, -1.7385),
            ]
        )
        * ang,
        (4, 4, 4),
        ['Li'],
        [1.0],
    )
]


def test_main():
    from pymbd.__main__ import ene, ene_expected

    assert ene == approx(ene_expected, rel=1e-10)


def test_argon_dimer_plain():
    ene = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)]).mbd_energy(
        [11, 11], [63.525, 63.525], [3.55, 3.55], 0.83, variant='plain'
    )
    assert ene == approx(-0.00024329110270970844, rel=1e-10)


@no_scalapack
def test_argon_dimer_dipole_matrix():
    dip = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)]).dipole_matrix('bare')
    assert (dip != 0).sum() == 6


def test_argon_dimer_rsscs():
    ene = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)]).mbd_energy_species(
        ['Ar', 'Ar'], [1, 1], 0.83
    )
    assert ene == approx(-0.0002462647623815428, rel=1e-10)


@no_scalapack
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


def test_benzene_dimer():
    mon1, mon2 = benzene_dimer
    dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
    enes = [
        MBDGeom(coords).mbd_energy_species(species, vol_ratios, 0.83)
        for coords, species, vol_ratios in (mon1, mon2, dim)
    ]
    ene_int = enes[2] - enes[1] - enes[0]
    assert ene_int == approx(-0.006312323931302544, rel=1e-10)


def test_benzene_gradients():
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


@no_scalapack
def test_benzene_dimer_python():
    mon1, mon2 = benzene_dimer
    dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
    enes = [
        MBDGeom(coords).mbd_energy_species(species, vol_ratios, 0.83)
        for coords, species, vol_ratios in (mon1, mon2, dim)
    ]
    ene_int = enes[2] - enes[1] - enes[0]
    assert ene_int == approx(-0.006312323931302544, rel=1e-10)


def test_benzene_gradients_plain():
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


def test_benzene_dimer_scs():
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


def test_benzene_dimer_ts():
    mon1, mon2 = benzene_dimer
    dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
    enes = [
        MBDGeom(coords).ts_energy_species(species, vol_ratios, 0.94)
        for coords, species, vol_ratios in (mon1, mon2, dim)
    ]
    ene_int = enes[2] - enes[1] - enes[0]
    assert ene_int == approx(-0.008490052683234028, rel=1e-10)


def test_benzene():
    coords, species, vol_ratios = benzene_dimer[0]
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    ene = MBDGeom(coords).mbd_energy(alpha_0, C6, R_vdw, 0.83, variant='plain')
    assert ene == approx(-0.007002398506090302, rel=1e-10)


@no_scalapack
def test_benzene_rpa():
    coords, species, vol_ratios = benzene_dimer[0]
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    ene = MBDGeom(coords, do_rpa=True).mbd_energy(
        alpha_0, C6, R_vdw, 0.83, variant='plain'
    )
    assert ene == approx(-0.007002398506090302, rel=1e-9)


@no_scalapack
def test_benzene_rpa_scaled():
    coords, species, vol_ratios = benzene_dimer[0]
    alpha_0, C6, R_vdw = from_volumes(species, vol_ratios)
    ene = MBDGeom(coords, do_rpa=True, rpa_rescale_eigs=True).mbd_energy(
        alpha_0, C6, R_vdw, 0.83, variant='plain'
    )
    assert ene != approx(-0.007002398506090302, rel=1e-9)
    assert ene == approx(-0.007002398506090302, rel=1e-7)


def test_ethylcarbamate():
    enes = [
        MBDGeom(coords, lattice, k_grid).mbd_energy_species(species, vol_ratios, 0.83)
        for coords, lattice, k_grid, species, vol_ratios in ethylcarbamate
    ]
    ene_int = enes[0] - 2 * enes[1]
    assert ene_int == approx(-0.037040868610822564, rel=1e-10)


def test_argon_crystal():
    coords, lattice, k_grid, species, vol_ratios = argon_crystal
    ene = MBDGeom(coords, lattice, k_grid).mbd_energy_species(species, vol_ratios, 0.83)
    assert ene == approx(-0.0021037562496878173, rel=1e-10)


@no_scalapack
def test_argon_crystal_rpa():
    coords, lattice, k_grid, species, vol_ratios = argon_crystal
    ene = MBDGeom(coords, lattice, k_grid, do_rpa=True).mbd_energy_species(
        species, vol_ratios, 0.83, variant='plain'
    )
    assert ene == approx(-0.0021036969146744147, rel=1e-10)


@no_scalapack
def test_argon_crystal_modes():
    coords, lattice, k_grid, species, vol_ratios = argon_crystal
    geom = MBDGeom(
        coords, lattice, custom_k_pts=[(0, 0, 0), (0.1, 0, 0)], get_spectrum=True
    )
    _, eigs, C = geom.mbd_energy_species(species, vol_ratios, 0.83)
    assert abs(C[:, :, 0].imag).sum() == approx(0)
    assert abs(C[:, :, 1].imag).sum() != approx(0)


def test_argon_crystal_gradients():
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


def test_lithium():
    with pytest.raises(MBDFortranException):
        [
            MBDGeom(coords, lattice, k_grid).mbd_energy_species(
                species, vol_ratios, 0.83
            )
            for coords, lattice, k_grid, species, vol_ratios in bulk_lithium
        ]


@no_scalapack
def test_ethylcarbamate_python():
    enes = [
        mbd_energy_species(
            coords, species, vol_ratios, 0.83, lattice=lattice, k_grid=k_grid
        )
        for coords, lattice, k_grid, species, vol_ratios in ethylcarbamate
    ]
    ene_int = enes[0] - 2 * enes[1]
    assert ene_int == approx(-0.037040868610822564, rel=1e-10)


def test_ethylcarbamate_scs():
    enes = [
        MBDGeom(coords, lattice, k_grid).mbd_energy_species(
            species, vol_ratios, 1, a=2.56, variant='scs'
        )
        for coords, lattice, k_grid, species, vol_ratios in ethylcarbamate
    ]
    ene_int = enes[0] - 2 * enes[1]
    assert ene_int == approx(-0.03633331132194684, rel=1e-10)


def test_ethylcarbamate_ts():
    enes = [
        MBDGeom(coords, lattice).ts_energy_species(species, vol_ratios, 0.83)
        for coords, lattice, _, species, vol_ratios in ethylcarbamate
    ]
    ene_int = enes[0] - 2 * enes[1]
    assert ene_int == approx(-0.052171811689150846, rel=1e-10)


@no_scalapack
def test_mbd_coulomb():
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
