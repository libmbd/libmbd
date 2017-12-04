# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from pytest import approx

import pymbd


def test_argon_dimer_plain():
    ene = pymbd.calculate([(0, 0, 0), (0, 0, 4)], [11, 11], [0.5, 0.5], [3.7, 3.7], 0.83, 6)
    assert ene == approx(-2.9343329137621055e-06, rel=1e-10)
