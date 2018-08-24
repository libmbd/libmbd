! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
program mbd_api_tests

use mbd_api, only: mbd_input, mbd_calculation, mbd_get_damping_parameters, &
    mbd_get_free_vdw_params

#ifdef WITH_MPI
use mpi
#endif

implicit none

integer, parameter :: dp = kind(0d0)
real(dp), parameter :: ang = 1.8897259886d0

type(mbd_input) :: inp
type(mbd_calculation) :: calc
real(dp) :: energy
real(dp), allocatable :: gradients(:, :)
logical :: failed

#ifdef WITH_MPI
integer :: err

call MPI_INIT(err)
#endif

failed = .false.

inp%free_values = mbd_get_free_vdw_params(['Ar', 'Ar'], 'ts')
call mbd_get_damping_parameters('pbe', inp%mbd_beta, inp%ts_d)
call calc%init(inp)
call calc%update_coords(reshape([0d0, 0d0, 0d0, 0d0, 0d0, 4d0*ang], [3, 2]))
call calc%update_vdw_params_custom([11d0, 11d0], [63.525d0, 63.525d0], [3.55d0, 3.55d0])
call calc%get_energy(energy)
call check('Ar2 energy', energy, -2.4329456747018696d-4, 1d-10)
allocate (gradients(3, 2))
call calc%get_gradients(gradients)
call check('Ar2 sum(abs(gradients))', sum(abs(gradients)), 2.3279742219399908d-4, 1d-10)
call calc%update_vdw_params_from_ratios([1d0, 1d0])
call calc%get_energy(energy)
call check('Ar2 energy 2', energy, -0.0002462647623815428d0, 1d-10)
call calc%destroy()

#ifdef WITH_MPI
call MPI_FINALIZE(err)
#endif

if (failed) stop 1

contains

subroutine check(label, val, ref, rel)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: val, ref, rel

    if (abs((val-ref)/ref) < rel) then
        write (6, *) label, ' OK'
    else
        failed = .true.
        write (6, *) label, val, ref, ' FAILED!'
    end if
end subroutine

end program

