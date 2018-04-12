! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_defaults

use mbd_common, only: dp

implicit none

real(dp), parameter :: TS_ENERGY_ACCURACY = 1d-6
real(dp), parameter :: TS_FORCES_ACCURACY = 1d-7
integer, parameter :: N_FREQUENCY_GRID = 15
real(dp), parameter :: K_GRID_SHIFT = 0.5d0
real(dp), parameter :: TS_DAMPING_D = 20.d0
real(dp), parameter :: MBD_DAMPING_A = 6.d0

end module
