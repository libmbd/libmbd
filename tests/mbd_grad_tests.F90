! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
program mbd_grad_tests

use mbd_constants
use mbd_geom, only: geom_t
use mbd_grad_test_cases
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

character(len=50) :: test_name

#ifdef WITH_MPI
integer :: err

call MPI_INIT(err)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
#else
rank = 0
#endif

call get_command_argument(1, test_name)
n_failed = 0
call exec_test(test_name)

#ifdef WITH_MPI
call MPI_FINALIZE(err)
#endif

if (n_failed /= 0) stop 1

contains

subroutine exec_test(test_name)
    character(len=*), intent(in) :: test_name

    allocate (geom)
    geom%log%level = MBD_LOG_LVL_INFO
    geom%parallel_mode = 'atoms'
    select case (test_name)
    case ('T_bare_deriv'); call test_T_bare_deriv()
    case ('T_GG_deriv_expl'); call test_T_GG_deriv_expl()
    case ('T_GG_deriv_impl'); call test_T_GG_deriv_impl()
    case ('T_erfc_deriv_expl'); call test_T_erfc_deriv_expl()
    case ('T_fermi_deriv_impl'); call test_T_fermi_deriv_impl()
    case ('mbd_deriv_expl'); call test_mbd_deriv_expl()
    case ('mbd_ewald_deriv_expl'); call test_mbd_ewald_deriv_expl()
    case ('mbd_ewald_deriv_stress'); call test_mbd_ewald_deriv_stress()
    case ('scs_deriv_expl'); call test_scs_deriv_expl()
    case ('scs_ewald_deriv_expl'); call test_scs_ewald_deriv_expl()
    case ('scs_ewald_deriv_stress'); call test_scs_ewald_deriv_stress()
    case ('scs_deriv_impl_alpha'); call test_scs_deriv_impl_alpha()
    case ('scs_ewald_deriv_impl_alpha'); call test_scs_ewald_deriv_impl_alpha()
    case ('scs_deriv_impl_vdw'); call test_scs_deriv_impl_vdw()
    case ('scs_ewald_deriv_impl_vdw'); call test_scs_ewald_deriv_impl_vdw()
    case ('mbd_deriv_impl_alpha'); call test_mbd_deriv_impl_alpha()
    case ('mbd_ewald_deriv_impl_alpha'); call test_mbd_ewald_deriv_impl_alpha()
    case ('mbd_deriv_impl_omega'); call test_mbd_deriv_impl_omega()
    case ('mbd_ewald_deriv_impl_omega'); call test_mbd_ewald_deriv_impl_omega()
    case ('mbd_deriv_impl_vdw'); call test_mbd_deriv_impl_vdw()
    case ('mbd_ewald_deriv_impl_vdw'); call test_mbd_ewald_deriv_impl_vdw()
    case ('mbd_ewald_deriv_impl_q'); call test_mbd_ewald_deriv_impl_q()
    case ('mbd_rsscs_deriv_expl'); call test_mbd_rsscs_deriv_expl()
    case ('mbd_rsscs_ewald_deriv_expl'); call test_mbd_rsscs_ewald_deriv_expl()
    case ('mbd_rsscs_ewald_deriv_stress'); call test_mbd_rsscs_ewald_deriv_stress()
    case ('mbd_rsscs_deriv_impl_alpha'); call test_mbd_rsscs_deriv_impl_alpha()
    case ('mbd_rsscs_deriv_impl_C6'); call test_mbd_rsscs_deriv_impl_C6()
    case ('mbd_rsscs_deriv_impl_vdw'); call test_mbd_rsscs_deriv_impl_vdw()
    case ('mbd_rsscs_ewald_deriv_impl_alpha'); call test_mbd_rsscs_ewald_deriv_impl_alpha()
    case ('mbd_rsscs_ewald_deriv_impl_C6'); call test_mbd_rsscs_ewald_deriv_impl_C6()
    case ('mbd_rsscs_ewald_deriv_impl_vdw'); call test_mbd_rsscs_ewald_deriv_impl_vdw()
    end select
    if (geom%exc%code /= 0) print *, 'Exception!'
    deallocate (geom)
end subroutine

end program
