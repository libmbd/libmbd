import numpy as np
import pymbd
import unittest


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
