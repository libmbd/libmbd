# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
import subprocess
import sys

import pytest

from pymbd.benchmark import parse


@pytest.mark.parametrize(
    'args', [[], ['--finite'], ['--method=ts']], ids=lambda x: ' '.join(x)
)
@pytest.mark.no_scalapack
def test_benchmark(args):
    stdout = subprocess.run(
        [sys.executable, '-u', '-m', 'pymbd.benchmark', *args],
        check=True,
        stdout=subprocess.PIPE,
    ).stdout.decode()
    assert parse(stdout)['energy']
