! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_constants

implicit none

integer, parameter :: dp = kind(0.d0)
real(dp), parameter :: pi = acos(-1.d0)

real(dp), parameter :: TS_ENERGY_ACCURACY = 1d-6
real(dp), parameter :: TS_FORCES_ACCURACY = 1d-7
integer, parameter :: N_FREQUENCY_GRID = 15
real(dp), parameter :: K_GRID_SHIFT = 0.5d0
real(dp), parameter :: TS_DAMPING_D = 20d0
real(dp), parameter :: MBD_DAMPING_A = 6d0

integer, parameter :: MBD_EXC_NEG_EIGVALS = 1
integer, parameter :: MBD_EXC_NEG_POL = 2
integer, parameter :: MBD_EXC_LINALG = 3
integer, parameter :: MBD_EXC_UNIMPL = 4
integer, parameter :: MBD_EXC_DAMPING = 5

real(dp), parameter :: ZERO_REAL = 0d0
complex(dp), parameter :: ZERO_COMPLEX = (0d0, 0d0)

end module
