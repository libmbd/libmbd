! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
program mbd_api_tests

use mbd_constants
use mbd, only: mbd_input_t, mbd_calc_t
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

logical :: failed
character(len=50) :: test_name

#ifdef WITH_MPI
integer :: err

call MPI_INIT(err)
#endif

call get_command_argument(1, test_name)
failed = .false.
call test(test_name)

#ifdef WITH_MPI
call MPI_FINALIZE(err)
#endif

if (failed) stop 1

contains

subroutine test(test_name)
    character(len=*), intent(in) :: test_name

    type(mbd_input_t) :: inp
    type(mbd_calc_t) :: calc
    real(dp) :: energy
    real(dp), allocatable :: gradients(:, :)
    integer :: code
    character(200) :: origin, msg

    inp%atom_types = ['Ar', 'Ar']
    inp%coords = reshape([0d0, 0d0, 0d0, 0d0, 0d0, 4d0*ang], [3, 2])
    inp%log_level = MBD_LOG_LVL_DEBUG
    select case (test_name)
    case ('exception')
        inp%xc = 'xxx'
        call calc%init(inp)
        call calc%get_exception(code, origin, msg)
        if (code /= MBD_EXC_DAMPING) failed = .true.
    case ('energy')
        inp%xc = 'pbe'
        call calc%init(inp)
        call calc%update_vdw_params_custom([11d0, 11d0], [63.525d0, 63.525d0], [3.55d0, 3.55d0])
        call calc%evaluate_vdw_method(energy)
        call check(energy, -2.4329456747018696d-4, 1d-10)
    case ('gradients')
        inp%xc = 'pbe'
        call calc%init(inp)
        call calc%update_vdw_params_custom([11d0, 11d0], [63.525d0, 63.525d0], [3.55d0, 3.55d0])
        call calc%evaluate_vdw_method(energy)
        allocate (gradients(3, 2))
        call calc%get_gradients(gradients)
        call check(sum(abs(gradients)), 2.3279742219399908d-4, 1d-10)
    case ('energy_from_ratios')
        inp%xc = 'pbe'
        call calc%init(inp)
        call calc%update_vdw_params_from_ratios([1d0, 1d0])
        call calc%evaluate_vdw_method(energy)
        call check(energy, -0.0002462647623815428d0, 1d-10)
    case ('ts_gradients')
        inp%xc = 'pbe'
        inp%method = 'ts'
        call calc%init(inp)
        call calc%update_vdw_params_from_ratios([1d0, 1d0])
        call calc%evaluate_vdw_method(energy)
        allocate (gradients(3, 2))
        call calc%get_gradients(gradients)
        call check(sum(abs(gradients)), 3.8405254013557403d-4, 1d-10)
    end select
    call calc%destroy()
end subroutine

subroutine check(val, ref, rel)
    real(dp), intent(in) :: val, ref, rel

    if (.not. abs((val-ref)/ref) < rel) then
        failed = .true.
    end if
end subroutine

end program
