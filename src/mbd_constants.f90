! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_constants
!! Constants used throughout.

implicit none

integer, parameter :: dp = kind(0.d0)
real(dp), parameter :: pi = acos(-1.d0)
real(dp), parameter :: ang = 1.8897259886d0
    !! Value of angstrom in atomic units

integer, parameter :: MBD_EXC_NEG_EIGVALS = 1
    !! Negative eigenvalue exception
integer, parameter :: MBD_EXC_NEG_POL = 2
    !! Negative polarizability exception
integer, parameter :: MBD_EXC_LINALG = 3
    !! Exception in LAPACK or ScaLAPACK
integer, parameter :: MBD_EXC_UNIMPL = 4
    !! Functionality is not implemented
integer, parameter :: MBD_EXC_DAMPING = 5
    !! Damping-function exception
integer, parameter :: MBD_EXC_INPUT = 6
    !! Invalid input

real(dp), parameter :: ZERO_REAL = 0d0
complex(dp), parameter :: ZERO_COMPLEX = (0d0, 0d0)
complex(dp), parameter :: IMI = (0d0, 1d0)

end module
