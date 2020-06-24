! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
program mbd_unit_tests

use mbd_constants
use mbd_geom, only: geom_t
use mbd_unit_test_cases
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

integer :: n_all

#ifdef WITH_MPI
integer :: err

call MPI_INIT(err)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
#else
    rank = 0
#endif

n_failed = 0
n_all = 0
call exec_test('T_bare derivative')
call exec_test('T_GG derivative explicit')
call exec_test('T_GG derivative implicit')
call exec_test('T_erfc derivative explicit')
call exec_test('T_fermi derivative implicit')
call exec_test('MBD derivative explicit')
call exec_test('MBD Ewald derivative explicit')
call exec_test('MBD Ewald derivative stress')
call exec_test('SCS derivative explicit')
call exec_test('SCS derivative implicit alpha')
call exec_test('SCS derivative implicit Rvdw')
call exec_test('SCS Ewald derivative explicit')
call exec_test('SCS Ewald derivative stress')
call exec_test('SCS Ewald derivative implicit alpha')
call exec_test('SCS Ewald derivative implicit Rvdw')
call exec_test('MBD derivative implicit alpha')
call exec_test('MBD derivative implicit omega')
call exec_test('MBD derivative implicit Rvdw')
call exec_test('MBD Ewald derivative implicit alpha')
call exec_test('MBD Ewald derivative implicit omega')
call exec_test('MBD Ewald derivative implicit Rvdw')
call exec_test('MBD Ewald derivative implicit q')
call exec_test('MBD@rsscs derivative explicit')
call exec_test('MBD@rsscs Ewald derivative explicit')
call exec_test('MBD@rsscs Ewald derivative stress')
call exec_test('MBD@rsscs derivative implicit alpha')
call exec_test('MBD@rsscs derivative implicit C6')
call exec_test('MBD@rsscs derivative implicit Rvdw')
call exec_test('MBD@rsscs Ewald derivative implicit alpha')
call exec_test('MBD@rsscs Ewald derivative implicit C6')
call exec_test('MBD@rsscs Ewald derivative implicit Rvdw')
if (rank == 0) write (6, *) &
    trim(tostr(n_failed)) // '/' // trim(tostr(n_all)) // ' tests failed'
if (n_failed /= 0) stop 1

#ifdef WITH_MPI
call MPI_FINALIZE(err)
#endif

contains

subroutine exec_test(test_name)
    character(len=*), intent(in) :: test_name

    if (rank == 0) write (6, '(A,A,A)', advance='no') &
        'Executing test "', test_name, '"... '
    allocate (geom)
    select case (test_name)
    case ('T_bare derivative'); call test_T_bare_deriv()
    case ('T_GG derivative explicit'); call test_T_GG_deriv_expl()
    case ('T_GG derivative implicit'); call test_T_GG_deriv_impl()
    case ('T_erfc derivative explicit'); call test_T_erfc_deriv_expl()
    case ('T_fermi derivative implicit'); call test_T_fermi_deriv_impl()
    case ('MBD derivative explicit'); call test_mbd_deriv_expl()
    case ('MBD Ewald derivative explicit'); call test_mbd_ewald_deriv_expl()
    case ('MBD Ewald derivative stress'); call test_mbd_ewald_deriv_stress()
    case ('SCS derivative explicit'); call test_scs_deriv_expl()
    case ('SCS Ewald derivative explicit'); call test_scs_ewald_deriv_expl()
    case ('SCS Ewald derivative stress'); call test_scs_ewald_deriv_stress()
    case ('SCS derivative implicit alpha'); call test_scs_deriv_impl_alpha()
    case ('SCS Ewald derivative implicit alpha'); call test_scs_ewald_deriv_impl_alpha()
    case ('SCS derivative implicit Rvdw'); call test_scs_deriv_impl_vdw()
    case ('SCS Ewald derivative implicit Rvdw'); call test_scs_ewald_deriv_impl_vdw()
    case ('MBD derivative implicit alpha'); call test_mbd_deriv_impl_alpha()
    case ('MBD Ewald derivative implicit alpha'); call test_mbd_ewald_deriv_impl_alpha()
    case ('MBD derivative implicit omega'); call test_mbd_deriv_impl_omega()
    case ('MBD Ewald derivative implicit omega'); call test_mbd_ewald_deriv_impl_omega()
    case ('MBD derivative implicit Rvdw'); call test_mbd_deriv_impl_vdw()
    case ('MBD Ewald derivative implicit Rvdw'); call test_mbd_ewald_deriv_impl_vdw()
    case ('MBD Ewald derivative implicit q'); call test_mbd_ewald_deriv_impl_q()
    case ('MBD@rsscs derivative explicit'); call test_mbd_rsscs_deriv_expl()
    case ('MBD@rsscs Ewald derivative explicit'); call test_mbd_rsscs_ewald_deriv_expl()
    case ('MBD@rsscs Ewald derivative stress'); call test_mbd_rsscs_ewald_deriv_stress()
    case ('MBD@rsscs derivative implicit alpha'); call test_mbd_rsscs_deriv_impl_alpha()
    case ('MBD@rsscs derivative implicit C6'); call test_mbd_rsscs_deriv_impl_C6()
    case ('MBD@rsscs derivative implicit Rvdw'); call test_mbd_rsscs_deriv_impl_vdw()
    case ('MBD@rsscs Ewald derivative implicit alpha'); call test_mbd_rsscs_ewald_deriv_impl_alpha()
    case ('MBD@rsscs Ewald derivative implicit C6'); call test_mbd_rsscs_ewald_deriv_impl_C6()
    case ('MBD@rsscs Ewald derivative implicit Rvdw'); call test_mbd_rsscs_ewald_deriv_impl_vdw()
    end select
    if (geom%exc%code /= 0) print *, 'Exception!'
    deallocate (geom)
    n_all = n_all + 1
end subroutine

end program
