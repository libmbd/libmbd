! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
program main

use mbd, only: run_tests

implicit none

integer :: err

external :: MPI_INIT, MPI_FINALIZE

call MPI_INIT(err)
call run_tests()
call MPI_FINALIZE(err)

end program
