from unittest import TestCase
import mbd


class TestMBD(TestCase):
    def test_argon_dimer(self):
        mbd.lib.init_grid(15)
        ene = mbd.mbd_rsscs(
            [[0, 0, 0], [4.0/mbd.bohr, 0, 0]],
            ['Ar', 'Ar'],
            [1., 1.],
            0.83
        )
        self.assertAlmostEqual(ene, -0.0002462647623815428, 10)
