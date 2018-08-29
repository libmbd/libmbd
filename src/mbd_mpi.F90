! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_mpi

#ifdef WITH_MPIFH
include 'mpif.h'
#else
use mpi
#endif

end module
