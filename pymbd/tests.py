import numpy as np
from numpy import sqrt
import unittest

import pymbd
from pymbd import bohr


benzene_dimer = [(
    np.array([
        (-1.047, -1.421, 0.000), (-1.454, -0.855, 1.206), (-1.454, -0.855, -1.206),
        (-2.266, 0.277, 1.206), (-2.671, 0.845, 0.000), (-2.266, 0.277, -1.206),
        (-1.133, -1.292, -2.142), (-2.582, 0.716, -2.143), (-3.303, 1.723, 0.000),
        (-2.582, 0.716, 2.143), (-1.133, -1.292, 2.142), (-0.406, -2.291, 0.000)
    ])/bohr,
    6*['C'] + 6*['H'],
    [0.825, 0.821, 0.821, 0.815, 0.814, 0.815, 0.624, 0.611, 0.610, 0.611, 0.624, 0.643]
), (
    np.array([
        (1.047, 1.421, 0.000), (1.454, 0.855, -1.206), (1.454, 0.855, 1.206),
        (2.266, -0.277, -1.206), (2.671, -0.845, 0.000), (2.266, -0.277, 1.206),
        (0.406, 2.291, 0.000), (1.133, 1.292, 2.142), (2.582, -0.716, 2.143),
        (3.303, -1.723, 0.000), (2.582, -0.716, -2.143), (1.133, 1.292, -2.142)
    ])/bohr,
    6*['C'] + 6*['H'],
    [0.825, 0.821, 0.821, 0.815, 0.814, 0.815, 0.643, 0.624, 0.611, 0.610, 0.611, 0.624]
)]
ethylcarbamate = [(
    np.array([
        (4.083, 5.700, 2.856), (0.568, 0.095, 4.217), (0.470, 4.774, 3.551),
        (4.181, 1.022, 3.522), (5.572, 5.587, 1.892), (-0.920, 0.209, 5.181),
        (3.663, 3.255, 2.585), (0.988, 2.541, 4.488), (3.834, 4.011, 0.979),
        (0.816, 1.785, 6.094), (2.223, 1.314, 1.108), (2.428, 4.481, 5.965),
        (1.177, 0.092, 0.406), (3.474, 5.703, 6.667), (4.911, 5.036, 2.573),
        (-0.260, 0.759, 4.500), (4.358, 3.787, 1.918), (0.293, 2.009, 5.155),
        (0.205, 1.729, 1.101), (4.446, 4.067, 5.972), (1.285, 0.947, 0.957),
        (3.366, 4.848, 6.116), (0.485, 2.901, 1.709), (4.165, 2.895, 5.364),
        (4.066, 1.426, 0.670), (0.585, 4.369, 6.403),
    ])/bohr,
    np.array([(5.008, 0.018, -0.070), (1.630, 6.759, 0.064), (-1.987, -0.981, 7.079)])/bohr,
    (2, 2, 2),
    list('HHHHHHHHHHHHHHCCCCCCNNOOOO'),
    [0.703, 0.703, 0.726, 0.726, 0.731, 0.731, 0.727, 0.727, 0.754, 0.754, 0.750,
     0.750, 0.755, 0.755, 0.809, 0.809, 0.827, 0.827, 0.834, 0.834, 0.840, 0.840,
     0.886, 0.886, 0.892, 0.892]
), (
    np.array([
        (4.088, 5.753, 2.783), (5.625, 5.562, 1.906), (3.652, 3.273, 2.592),
        (3.854, 3.998, 0.981), (5.422, 4.834, 3.521), (6.213, 0.125, 0.386),
        (7.201, 1.360, 1.112), (4.913, 5.058, 2.573), (4.366, 3.792, 1.934),
        (5.167, 1.729, 1.084), (6.291, 0.963, 0.938), (4.042, 1.399, 0.752),
        (5.490, 2.915, 1.682)
    ])/bohr,
    None,
    None,
    list('HHHHHHHCCCNOO'),
    [0.581, 0.607, 0.642, 0.646, 0.607, 0.596, 0.597, 0.762, 0.799, 0.845,
     0.824, 0.974, 0.896]
)]


class TestFreqGrid(unittest.TestCase):
    pass


for n in range(15, 40, 5):
    def test(self, n=n):
        grid1 = -pymbd.lib.gauss_legendre(n)[0]
        grid2 = np.polynomial.legendre.leggauss(n)[0]
        self.assertAlmostEqual(0, np.linalg.norm(grid1-grid2), 8)
    if n > 25:
        test = unittest.expectedFailure(test)
    setattr(TestFreqGrid, 'test_frequency_grid_{}'.format(n), test)


class TestMBD(unittest.TestCase):
    def setUp(self):
        pymbd.lib.init_grid(15)

    def test_argon_dimer(self):
        ene = pymbd.mbd_rsscs(
            [[0, 0, 0], [4.0/pymbd.bohr, 0, 0]],
            ['Ar', 'Ar'],
            [1., 1.],
            0.83
        )
        self.assertAlmostEqual(ene, -0.0002462647623815428, 10)

    def test_argon_dimer_plain(self):
        ene = pymbd.lib.get_single_mbd_energy(
            '', 'fermi,dip', [(0, 0, 0), (0, 0, 4/bohr)], [11, 11], [0.7, 0.7],
            r_vdw=[3.55, 3.55], beta=0.83, a=6
        )[0]
        self.assertAlmostEqual(ene, -0.00024329110270970844, 10)

    def test_benzene_dimer(self):
        mon1, mon2 = benzene_dimer
        dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
        enes = [
            pymbd.mbd_rsscs(coords, species, vols, 0.83)
            for coords, species, vols in (mon1, mon2, dim)
        ]
        ene_int = enes[2]-enes[1]-enes[0]
        self.assertAlmostEqual(ene_int, -0.006312323931302544, 10)

    def test_benzene_dimer_scs(self):
        mon1, mon2 = benzene_dimer
        dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
        enes = [
            pymbd.mbd_scs(coords, species, vols, 2.56)
            for coords, species, vols in (mon1, mon2, dim)
        ]
        ene_int = enes[2]-enes[1]-enes[0]
        self.assertAlmostEqual(ene_int, -0.004517208661233951, 10)

    def test_ethylcarbamate(self):
        enes = [
            pymbd.mbd_rsscs(coords, species, vols, 0.83, lattice=lattice, k_grid=k_grid)
            for coords, lattice, k_grid, species, vols in ethylcarbamate
        ]
        ene_int = enes[0]-2*enes[1]
        self.assertAlmostEqual(ene_int, -0.037040868610822564, 10)

    def test_ethylcarbamate_scs(self):
        enes = [
            pymbd.mbd_scs(coords, species, vols, 2.56, lattice=lattice, k_grid=k_grid)
            for coords, lattice, k_grid, species, vols in ethylcarbamate
        ]
        ene_int = enes[0]-2*enes[1]
        self.assertAlmostEqual(ene_int, -0.03185124226203406, 10)

    def tearDown(self):
        pymbd.lib.destroy_grid()


class TestCoulomb(unittest.TestCase):
    def test_basic(self):
        a = 1/sqrt(2)
        C = np.array([
            [ 0,  0, -a,  0, -a,  0],  # noqa
            [ 0,  a,  0,  a,  0,  0],  # noqa
            [-a,  0,  0,  0,  0, -a],  # noqa
            [ 0,  0,  a,  0, -a,  0],  # noqa
            [ 0, -a,  0,  a,  0,  0],  # noqa
            [-a,  0,  0,  0,  0,  a],  # noqa
        ], dtype=float)
        coords = [[0, 0, 0], [0, 0, 1.1]]
        n = len(coords)
        e1, eatt, erep = pymbd.lib_coul.fullcoulomb(
            C, coords, n*[1], n*[1], 3*n*[2], n*[2]
        )
        self.assertAlmostEqual(e1, -0.058346493647401076, 10)
        self.assertAlmostEqual(eatt, 1.767623829681647, 10)
        self.assertAlmostEqual(erep, 0.800186426943337, 10)
