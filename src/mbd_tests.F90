! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#define MODULE_UNIT_TESTS
#include "mbd.F90"

program mbd_tests

use mbd
use mbd_common, only: diff7, findval

#ifdef WITH_MPI
use mpi
#endif

implicit none

integer :: n_failed, n_all, rank
type(mbd_calc), target :: calc

#ifdef WITH_MPI
integer :: err

call MPI_INIT(err)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
#else
    rank = 0
#endif

call calc%init_grid()
n_failed = 0
n_all = 0
call exec_test('T_bare derivative')
call exec_test('T_GG derivative explicit')
call exec_test('T_GG derivative implicit')
call exec_test('T_fermi derivative implicit')
call exec_test('MBD derivative explicit')
call exec_test('SCS derivative explicit')
call exec_test('SCS derivative implicit alpha')
call exec_test('SCS derivative implicit Rvdw')
call exec_test('MBD derivative implicit alpha')
call exec_test('MBD derivative implicit C6')
call exec_test('MBD derivative implicit Rvdw')
call exec_test('MBD@rsscs derivative explicit')
call exec_test('MBD@rsscs derivative implicit alpha')
call exec_test('MBD@rsscs derivative implicit C6')
call exec_test('MBD@rsscs derivative implicit Rvdw')
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
    select case (test_name)
    case ('T_bare derivative'); call test_T_bare_deriv()
    case ('T_GG derivative explicit'); call test_T_GG_deriv_expl()
    case ('T_GG derivative implicit'); call test_T_GG_deriv_impl()
    case ('T_fermi derivative implicit'); call test_T_fermi_deriv_impl()
    case ('MBD derivative explicit'); call test_mbd_deriv_expl()
    case ('SCS derivative explicit'); call test_scs_deriv_expl()
    case ('SCS derivative implicit alpha'); call test_scs_deriv_impl_alpha()
    case ('SCS derivative implicit Rvdw'); call test_scs_deriv_impl_vdw()
    case ('MBD derivative implicit alpha'); call test_mbd_deriv_impl_alpha()
    case ('MBD derivative implicit C6'); call test_mbd_deriv_impl_C6()
    case ('MBD derivative implicit Rvdw'); call test_mbd_deriv_impl_vdw()
    case ('MBD@rsscs derivative explicit'); call test_mbd_rsscs_deriv_expl()
    case ('MBD@rsscs derivative implicit alpha'); call test_mbd_rsscs_deriv_impl_alpha()
    case ('MBD@rsscs derivative implicit C6'); call test_mbd_rsscs_deriv_impl_C6()
    case ('MBD@rsscs derivative implicit Rvdw'); call test_mbd_rsscs_deriv_impl_vdw()
    end select
    n_all = n_all + 1
end subroutine

logical function failed(diff, thre)
    real(dp), intent(in) :: diff, thre

    failed = abs(diff) > thre
    if (rank == 0) write (6, '(A,G10.3,A,G10.3,A)', advance='no') &
        'diff:', diff, ', threshold:', thre, ': '
    if (failed) n_failed = n_failed + 1
    if (rank == 0) then
        if (failed) then
            write (6, *) 'FAILED!'
        else
            write (6, *) 'OK'
        end if
    end if
end function

subroutine test_T_bare_deriv()
    real(dp) :: r(3), r_diff(3)
    type(mat33) :: T
    real(dp) :: diff(3, 3)
    real(dp) :: T_diff_anl(3, 3, 3)
    real(dp) :: T_diff_num(3, 3, -3:3)
    integer :: a, b, c, i_step
    real(dp) :: delta

    delta = 1d-2
    r = [1.12d0, -2.12d0, 0.12d0]
    T = T_bare_v2(r, deriv=.true.)
    T_diff_anl = T%dr(:, :, :)
    diff = 0d0
    do c = 1, 3
        do i_step = -3, 3
            if (i_step == 0) cycle
            r_diff = r
            r_diff(c) = r_diff(c)+i_step*delta
            T = T_bare_v2(r_diff, deriv=.false.)
            T_diff_num(:, :, i_step) = T%val
        end do
        forall (a = 1:3, b = 1:3)
            T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
        end forall
        diff = max(diff, abs(T_diff_num(:, :, 0)-T_diff_anl(:, :, c))/T_diff_num(:, :, 0))
    end do
    if (failed(maxval(abs(diff)), 1d-10)) then
    end if
end subroutine test_T_bare_deriv

subroutine test_T_GG_deriv_expl()
    real(dp) :: r(3), r_diff(3)
    type(mat33) :: T
    real(dp) :: diff(3, 3)
    real(dp) :: T_diff_anl(3, 3, 3)
    real(dp) :: T_diff_num(3, 3, -3:3)
    integer :: a, b, c, i_step
    real(dp) :: delta
    real(dp) :: sigma

    delta = 1d-2
    r = [1.02d0, -2.22d0, 0.15d0]
    sigma = 1.2d0
    T = T_erf_coulomb(r, sigma, deriv=.true.)
    T_diff_anl = T%dr
    diff = 0d0
    do c = 1, 3
        do i_step = -3, 3
            if (i_step == 0) cycle
            r_diff = r
            r_diff(c) = r_diff(c)+i_step*delta
            T = T_erf_coulomb(r_diff, sigma, deriv=.false.)
            T_diff_num(:, :, i_step) = T%val
        end do
        forall (a = 1:3, b = 1:3)
            T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
        end forall
        diff = max(diff, abs(T_diff_num(:, :, 0)-T_diff_anl(:, :, c))/T_diff_num(:, :, 0))
    end do
    if (failed(maxval(abs(diff)), 1d-10)) then
    end if
end subroutine test_T_GG_deriv_expl

subroutine test_T_GG_deriv_impl()
    real(dp) :: r(3)
    type(mat33) :: T
    real(dp) :: diff(3, 3)
    real(dp) :: T_diff_anl(3, 3)
    real(dp) :: T_diff_num(3, 3, -3:3)
    integer :: a, b, i_step
    real(dp) :: delta
    real(dp) :: sigma, dsigma_dr, sigma_diff

    delta = 1d-3
    r = [1.02d0, -2.22d0, 0.15d0]
    sigma = 1.2d0
    dsigma_dr = -0.3d0
    T = T_erf_coulomb(r, sigma, deriv=.true.)
    T_diff_anl = T%dsigma(:, :)*dsigma_dr
    do i_step = -3, 3
        if (i_step == 0) cycle
        sigma_diff = sigma+i_step*delta*dsigma_dr
        T = T_erf_coulomb(r, sigma_diff, deriv=.false.)
        T_diff_num(:, :, i_step) = T%val
    end do
    forall (a = 1:3, b = 1:3)
        T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
    end forall
    diff = (T_diff_num(:, :, 0)-T_diff_anl)/T_diff_num(:, :, 0)
    if (failed(maxval(abs(diff)), 1d-10)) then
        call print_matrix('delta dTGG', diff)
    end if
end subroutine test_T_GG_deriv_impl

subroutine test_T_fermi_deriv_impl()
    real(dp) :: r(3)
    type(mat33) :: T
    real(dp) :: diff(3, 3)
    real(dp) :: T_diff_anl(3, 3)
    real(dp) :: T_diff_num(3, 3, -3:3)
    integer :: a, b, i_step
    real(dp) :: delta
    real(dp) :: rvdw, drvdw_dr, rvdw_diff

    delta = 1d-3
    r = [1.02d0, -2.22d0, 0.15d0]
    rvdw = 2.5d0
    drvdw_dr = -0.3d0
    T = T_damped(damping_fermi(r, rvdw, 6d0, .true.), T_bare_v2(r, .true.))
    T_diff_anl = T%dvdw(:, :)*drvdw_dr
    do i_step = -3, 3
        if (i_step == 0) cycle
        rvdw_diff =rvdw+i_step*delta*drvdw_dr
        T = T_damped(damping_fermi(r, rvdw_diff, 6d0, .false.), T_bare_v2(r, .false.))
        T_diff_num(:, :, i_step) = T%val
    end do
    forall (a = 1:3, b = 1:3)
        T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
    end forall
    diff = (T_diff_num(:, :, 0)-T_diff_anl)/T_diff_num(:, :, 0)
    if (failed(maxval(abs(diff)), 1d-10)) then
        call print_matrix('delta dTfermi', diff)
    end if
end subroutine test_T_fermi_deriv_impl

subroutine test_mbd_deriv_expl()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: gradients(:, :)
    real(dp), allocatable :: diff(:, :)
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
    type(mbd_result) :: res(-3:3)
    type(mbd_gradients) :: dene
    real(dp), allocatable :: gradients_anl(:, :)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 0.01d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms, 3))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    allocate (dene%dcoords(n_atoms, 3))
    res(0) = get_single_mbd_energy(sys, alpha_0, C6, damp, dene)
    gradients_anl = dene%dcoords
    deallocate (dene%dcoords)
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                sys%coords = coords
                sys%coords(i_xyz, i_atom) = sys%coords(i_xyz, i_atom)+i_step*delta
                res(i_step) = get_single_mbd_energy(sys, alpha_0, C6, damp, dene)
            end do
            gradients(i_atom, i_xyz) = diff7(res%energy, delta)
        end do
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', diff)
    end if
end subroutine test_mbd_deriv_expl

subroutine test_scs_deriv_expl()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: gradients(:, :, :), gradients_anl(:, :, :)
    real(dp), allocatable :: diff(:, :, :)
    real(dp), allocatable :: alpha_0(:)
    integer :: i_atom, n_atoms, i_xyz, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    real(dp), allocatable :: alpha_scs(:, :)
    type(mbd_gradients), allocatable :: dalpha_scs(:)

    delta = 0.05d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    my_nratoms = size(sys%blacs%i_atom)
    my_ncatoms = size(sys%blacs%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms, 3))
    allocate (gradients_anl(my_nratoms, my_ncatoms, 3))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    do my_i_atom = 1, my_nratoms
        allocate (dalpha_scs(my_i_atom)%dcoords(my_ncatoms, 3))
    end do
    alpha_scs(:, 0) = run_scs(sys, alpha_0, damp, dalpha_scs)
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :, :) = dalpha_scs(my_i_atom)%dcoords
        deallocate (dalpha_scs(my_i_atom)%dcoords)
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(sys%blacs%j_atom, j_atom)
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                sys%coords = coords
                sys%coords(i_xyz, j_atom) = sys%coords(i_xyz, j_atom) + &
                    i_step*delta
                alpha_scs(:, i_step) = run_scs(sys, alpha_0, damp, dalpha_scs)
            end do
            if (my_j_atom > 0) then
                do my_i_atom = 1, my_nratoms
                    i_atom = sys%blacs%i_atom(my_i_atom)
                    gradients(my_i_atom, my_j_atom, i_xyz) = &
                        diff7(alpha_scs(i_atom, :), delta)
                end do
            end if
        end do
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-5)) then
        call print_matrix('diff x', diff(:, :, 1))
        call print_matrix('diff y', diff(:, :, 2))
        call print_matrix('diff z', diff(:, :, 3))
    end if
end subroutine test_scs_deriv_expl

subroutine test_scs_deriv_impl_alpha
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:, :), &
        gradients_anl(:, :), diff(:, :), alpha_0(:), alpha_0_diff(:), &
        alpha_scs(:, :)
    integer :: i_atom, n_atoms, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    type(mbd_gradients), allocatable :: dalpha_scs(:)

    delta = 0.1d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    my_nratoms = size(sys%blacs%i_atom)
    my_ncatoms = size(sys%blacs%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms))
    allocate (gradients_anl(my_nratoms, my_ncatoms))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    do my_i_atom = 1, my_nratoms
        allocate (dalpha_scs(my_i_atom)%dalpha(my_ncatoms))
    end do
    alpha_scs(:, 0) = run_scs(sys, alpha_0, damp, dalpha_scs)
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :) = dalpha_scs(my_i_atom)%dalpha
        deallocate (dalpha_scs(my_i_atom)%dalpha)
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(sys%blacs%j_atom, j_atom)
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(j_atom) = alpha_0_diff(j_atom) + i_step*delta
            alpha_scs(:, i_step) = run_scs(sys, alpha_0_diff, damp, dalpha_scs)
        end do
        if (my_j_atom > 0) then
            do my_i_atom = 1, my_nratoms
                i_atom = sys%blacs%i_atom(my_i_atom)
                gradients(my_i_atom, my_j_atom) = diff7(alpha_scs(i_atom, :), delta)
            end do
    end if
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        call print_matrix('diff', diff)
    end if
end subroutine test_scs_deriv_impl_alpha

subroutine test_scs_deriv_impl_vdw
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:, :), &
        gradients_anl(:, :), diff(:, :), alpha_0(:), alpha_scs(:, :), rvdw(:)
    integer :: i_atom, n_atoms, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    type(mbd_gradients), allocatable :: dalpha_scs(:)

    delta = 0.1d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    my_nratoms = size(sys%blacs%i_atom)
    my_ncatoms = size(sys%blacs%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms))
    allocate (gradients_anl(my_nratoms, my_ncatoms))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    rvdw = [3.55d0, 3.5d0, 3.56d0]
    damp%r_vdw = rvdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    do my_i_atom = 1, my_nratoms
        allocate (dalpha_scs(my_i_atom)%dr_vdw(my_ncatoms))
    end do
    alpha_scs(:, 0) = run_scs(sys, alpha_0, damp, dalpha_scs)
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :) = dalpha_scs(my_i_atom)%dr_vdw
        deallocate (dalpha_scs(my_i_atom)%dr_vdw)
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(sys%blacs%j_atom, j_atom)
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = rvdw
            damp%r_vdw(j_atom) = damp%r_vdw(j_atom) + i_step*delta
            alpha_scs(:, i_step) = run_scs(sys, alpha_0, damp, dalpha_scs)
        end do
        if (my_j_atom > 0) then
            do my_i_atom = 1, my_nratoms
                i_atom = sys%blacs%i_atom(my_i_atom)
                gradients(my_i_atom, my_j_atom) = diff7(alpha_scs(i_atom, :), delta)
            end do
        end if
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        call print_matrix('diff', diff(:, :))
    end if
end subroutine test_scs_deriv_impl_vdw

subroutine test_mbd_deriv_impl_alpha()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), alpha_0_diff(:), C6(:)
    type(mbd_result) :: res(-3:3)
    type(mbd_gradients) :: dene
    integer :: i_atom, n_atoms, i_step

    delta = 0.1d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    allocate (dene%dalpha(n_atoms))
    res(0) = get_single_mbd_energy(sys, alpha_0, C6, damp, dene)
    call move_alloc(dene%dalpha, gradients_anl)
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(i_atom) = alpha_0_diff(i_atom) + i_step*delta
            res(i_step) = get_single_mbd_energy(sys, alpha_0_diff, C6, damp, dene)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        call print_matrix('diff', reshape(diff, [n_atoms, 1]))
    end if
end subroutine test_mbd_deriv_impl_alpha

subroutine test_mbd_deriv_impl_C6()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), C6_diff(:), C6(:)
    type(mbd_result) :: res(-3:3)
    type(mbd_gradients) :: dene
    integer :: i_atom, n_atoms, i_step

    delta = 0.03d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    allocate (dene%dC6(n_atoms))
    res(0) = get_single_mbd_energy(sys, alpha_0, C6, damp, dene)
    call move_alloc(dene%dC6, gradients_anl)
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            C6_diff = C6
            C6_diff(i_atom) = C6_diff(i_atom) + i_step*delta
            res(i_step) = get_single_mbd_energy(sys, alpha_0, C6_diff, damp, dene)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 2d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine test_mbd_deriv_impl_C6

subroutine test_mbd_deriv_impl_vdw()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), C6(:), r_vdw(:)
    type(mbd_result) :: res(-3:3)
    type(mbd_gradients) :: dene
    integer :: i_atom, n_atoms, i_step

    delta = 1d-3
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    damp%version = 'fermi,dip'
    r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%r_vdw = r_vdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    allocate (dene%dr_vdw(n_atoms))
    res(0) = get_single_mbd_energy(sys, alpha_0, C6, damp, dene)
    call move_alloc(dene%dr_vdw, gradients_anl)
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = r_vdw
            damp%r_vdw(i_atom) = damp%r_vdw(i_atom) + i_step*delta
            res(i_step) = get_single_mbd_energy(sys, alpha_0, C6, damp, dene)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine test_mbd_deriv_impl_vdw

subroutine test_mbd_rsscs_deriv_expl()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: gradients(:, :), gradients_anl(:, :)
    real(dp), allocatable :: diff(:, :)
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
    type(mbd_result) :: res(-3:3)
    type(mbd_gradients) :: dene
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 0.01d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms, 3))
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    coords(1, 3) = 1d0
    sys%coords = coords
    call sys%init(calc)
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    allocate (dene%dcoords(n_atoms, 3))
    res(0) = mbd_scs_energy(sys, 'rsscs', alpha_0, C6, damp, dene)
    call move_alloc(dene%dcoords, gradients_anl)
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                sys%coords = coords
                sys%coords(i_xyz, i_atom) = sys%coords(i_xyz, i_atom) + &
                    i_step*delta
                res(i_step) = mbd_scs_energy(sys, 'rsscs', alpha_0, C6, damp, dene)
            end do
            gradients(i_atom, i_xyz) = diff7(res%energy, delta)
        end do
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', diff)
    end if
end subroutine test_mbd_rsscs_deriv_expl

subroutine test_mbd_rsscs_deriv_impl_alpha()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), alpha_0_diff(:), C6(:)
    type(mbd_result) :: res(-3:3)
    type(mbd_gradients) :: dene
    integer :: i_atom, n_atoms, i_step

    delta = 3d-2
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    allocate (dene%dalpha(n_atoms))
    res(0) = mbd_scs_energy(sys, 'rsscs', alpha_0, C6, damp, dene)
    call move_alloc(dene%dalpha, gradients_anl)
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(i_atom) = alpha_0_diff(i_atom) + i_step*delta
            res(i_step) = mbd_scs_energy(sys, 'rsscs', alpha_0_diff, C6, damp, dene)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-7)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine test_mbd_rsscs_deriv_impl_alpha

subroutine test_mbd_rsscs_deriv_impl_C6()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), C6_diff(:), C6(:)
    type(mbd_result) :: res(-3:3)
    type(mbd_gradients) :: dene
    integer :: i_atom, n_atoms, i_step

    delta = 0.01d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    allocate (dene%dC6(n_atoms))
    res(0) = mbd_scs_energy(sys, 'rsscs', alpha_0, C6, damp, dene)
    call move_alloc(dene%dC6, gradients_anl)
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            C6_diff = C6
            C6_diff(i_atom) = C6_diff(i_atom) + i_step*delta
            res(i_step) = mbd_scs_energy(sys, 'rsscs', alpha_0, C6_diff, damp, dene)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 5d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine test_mbd_rsscs_deriv_impl_C6

subroutine test_mbd_rsscs_deriv_impl_vdw()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), C6(:), r_vdw(:)
    type(mbd_result) :: res(-3:3)
    type(mbd_gradients) :: dene
    integer :: i_atom, n_atoms, i_step

    delta = 1d-2
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    sys%coords = coords
    call sys%init(calc)
    damp%version = 'fermi,dip'
    r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%r_vdw = r_vdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    allocate (dene%dr_vdw(n_atoms))
    res(0) = mbd_scs_energy(sys, 'rsscs', alpha_0, C6, damp, dene)
    call move_alloc(dene%dr_vdw, gradients_anl)
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = r_vdw
            damp%r_vdw(i_atom) = damp%r_vdw(i_atom) + i_step*delta
            res(i_step) = mbd_scs_energy(sys, 'rsscs', alpha_0, C6, damp, dene)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call sys%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine test_mbd_rsscs_deriv_impl_vdw

end program
