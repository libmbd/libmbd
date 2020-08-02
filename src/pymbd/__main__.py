# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
import sys

from pymbd import ang
from pymbd.fortran import MBDGeom, with_mpi

ene_expected = -0.0002462647623815428
ene = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)]).mbd_energy_species(
    ['Ar', 'Ar'], [1, 1], 0.83
)
if with_mpi:
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.Get_rank()
else:
    rank = 0
if rank == 0:
    print(f'Expected energy:   {ene_expected}')
    print(f'Calculated energy: {ene}')
if ene - ene_expected > 1e-10:
    sys.exit(1)
