# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
import numpy as np
from pytest import approx

import pymbd
from pymbd import ang

benzene_dimer = [(
    np.array([
        (-1.047, -1.421, 0.000), (-1.454, -0.855, 1.206), (-1.454, -0.855, -1.206),
        (-2.266, 0.277, 1.206), (-2.671, 0.845, 0.000), (-2.266, 0.277, -1.206),
        (-1.133, -1.292, -2.142), (-2.582, 0.716, -2.143), (-3.303, 1.723, 0.000),
        (-2.582, 0.716, 2.143), (-1.133, -1.292, 2.142), (-0.406, -2.291, 0.000)
    ])*ang,
    6*['C'] + 6*['H'],
    [0.825, 0.821, 0.821, 0.815, 0.814, 0.815, 0.624, 0.611, 0.610, 0.611, 0.624, 0.643]
), (
    np.array([
        (1.047, 1.421, 0.000), (1.454, 0.855, -1.206), (1.454, 0.855, 1.206),
        (2.266, -0.277, -1.206), (2.671, -0.845, 0.000), (2.266, -0.277, 1.206),
        (0.406, 2.291, 0.000), (1.133, 1.292, 2.142), (2.582, -0.716, 2.143),
        (3.303, -1.723, 0.000), (2.582, -0.716, -2.143), (1.133, 1.292, -2.142)
    ])*ang,
    6*['C'] + 6*['H'],
    [0.825, 0.821, 0.821, 0.815, 0.814, 0.815, 0.643, 0.624, 0.611, 0.610, 0.611, 0.624]
)]


def test_argon_dimer_plain():
    ene = pymbd.mbd_energy(
        [(0, 0, 0), (0, 0, 4*ang)], [11, 11], [0.7, 0.7], [3.55, 3.55], 0.83,
        func='calc_mbd_energy'
    )
    assert ene == approx(-0.00024329110270970844, rel=1e-10)


def test_argon_dimer_rsscs():
    ene = pymbd.mbd_energy_species(
        [(0, 0, 0), (0, 0, 4*ang)], ['Ar', 'Ar'], [1, 1], 0.83
    )
    assert ene == approx(-0.0002462647623815428, rel=1e-10)


def test_benzene_dimer():
    mon1, mon2 = benzene_dimer
    dim = (np.vstack((mon1[0], mon2[0])), mon1[1] + mon2[1], mon1[2] + mon2[2])
    enes = [
        pymbd.mbd_energy_species(coords, species, vols, 0.83)
        for coords, species, vols in (mon1, mon2, dim)
    ]
    ene_int = enes[2]-enes[1]-enes[0]
    assert ene_int == approx(-0.006312323931302544, rel=1e-10)
