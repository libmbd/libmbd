# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from pytest import approx

import pymbd
from pymbd import ang


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
