import numpy as np
from numpy import sqrt
import unittest

import pymbd


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
