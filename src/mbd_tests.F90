! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#define MODULE_UNIT_TESTS
#include "mbd.F90"

program mbd_tests

use mbd
use mbd_common, only: diff3, diff5

implicit none

integer :: err

external :: MPI_INIT, MPI_FINALIZE

integer :: n_failed, n_all
type(mbd_calc), target :: calc

call MPI_INIT(err)

call init_grid(calc)
n_failed = 0
n_all = 0
call exec_test('T_bare derivative')
call exec_test('T_GG derivative explicit')
call exec_test('T_GG derivative implicit')
call exec_test('T_fermi derivative implicit')
call exec_test('MBD derivative explicit')
call exec_test('SCS derivative explicit')
call exec_test('MBD derivative implicit alpha')
call exec_test('MBD derivative implicit C6')
call exec_test('MBD derivative implicit Rvdw')
write (6, *) &
    trim(tostr(n_failed)) // '/' // trim(tostr(n_all)) // ' tests failed'
if (n_failed /= 0) stop 1

call MPI_FINALIZE(err)

contains

subroutine exec_test(test_name)
    character(len=*), intent(in) :: test_name

    integer :: n_failed_in

    write (6, '(A,A,A)', advance='no') 'Executing test "', test_name, '"... '
    n_failed_in = n_failed
    select case (test_name)
    case ('T_bare derivative'); call test_T_bare_deriv()
    case ('T_GG derivative explicit'); call test_T_GG_deriv_expl()
    case ('T_GG derivative implicit'); call test_T_GG_deriv_impl()
    case ('T_fermi derivative implicit'); call test_T_fermi_deriv_impl()
    case ('MBD derivative explicit'); call test_mbd_deriv_expl()
    case ('SCS derivative explicit'); call test_scs_deriv_expl()
    case ('MBD derivative implicit alpha'); call test_mbd_deriv_impl_alpha()
    case ('MBD derivative implicit C6'); call test_mbd_deriv_impl_C6()
    case ('MBD derivative implicit Rvdw'); call test_mbd_deriv_impl_vdw()
    end select
    n_all = n_all + 1
    if (n_failed == n_failed_in) write (6, *) 'OK'
end subroutine

subroutine failed()
    n_failed = n_failed + 1
    write (6, *) 'FAILED!'
end subroutine

subroutine test_T_bare_deriv()
    real(dp) :: r(3), r_diff(3)
    type(mat33) :: T
    real(dp) :: diff(3, 3)
    real(dp) :: T_diff_anl(3, 3, 3)
    real(dp) :: T_diff_num(3, 3, -2:2)
    integer :: a, b, c, i_step
    real(dp) :: delta

    delta = 1d-3
    r = [1.12d0, -2.12d0, 0.12d0]
    T = T_bare_v2(r, deriv=.true.)
    T_diff_anl = T%dr(:, :, :)
    do c = 1, 3
        do i_step = -2, 2
            if (i_step == 0) cycle
            r_diff = r
            r_diff(c) = r_diff(c)+i_step*delta
            T = T_bare_v2(r_diff, deriv=.false.)
            T_diff_num(:, :, i_step) = T%val
        end do
        forall (a = 1:3, b = 1:3)
            T_diff_num(a, b, 0) = diff5(T_diff_num(a, b, :), delta)
        end forall
        diff = (T_diff_num(:, :, 0)-T_diff_anl(:, :, c))/T_diff_num(:, :, 0)
        if (any(abs(diff) > 1d-10)) then
            call failed()
            call print_matrix('delta dT(:, :, ' // trim(tostr(c)) // ')', diff)
        end if
    end do
end subroutine test_T_bare_deriv

subroutine test_T_GG_deriv_expl()
    real(dp) :: r(3), r_diff(3)
    type(mat33) :: T
    real(dp) :: diff(3, 3)
    real(dp) :: T_diff_anl(3, 3, 3)
    real(dp) :: T_diff_num(3, 3, -2:2)
    integer :: a, b, c, i_step
    real(dp) :: delta
    real(dp) :: sigma

    delta = 1d-3
    r = [1.02d0, -2.22d0, 0.15d0]
    sigma = 1.2d0
    T = T_erf_coulomb(r, sigma, deriv=.true.)
    T_diff_anl = T%dr
    do c = 1, 3
        do i_step = -2, 2
            if (i_step == 0) cycle
            r_diff = r
            r_diff(c) = r_diff(c)+i_step*delta
            T = T_erf_coulomb(r_diff, sigma, deriv=.false.)
            T_diff_num(:, :, i_step) = T%val
        end do
        forall (a = 1:3, b = 1:3)
            T_diff_num(a, b, 0) = diff5(T_diff_num(a, b, :), delta)
        end forall
        diff = (T_diff_num(:, :, 0)-T_diff_anl(:, :, c))/T_diff_num(:, :, 0)
        if (any(abs(diff) > 1d-10)) then
            call failed()
            call print_matrix('delta dTGG_{ab,' // trim(tostr(c)) // '}', diff)
        end if
    end do
end subroutine test_T_GG_deriv_expl

subroutine test_T_GG_deriv_impl()
    real(dp) :: r(3)
    type(mat33) :: T
    real(dp) :: diff(3, 3)
    real(dp) :: T_diff_anl(3, 3)
    real(dp) :: T_diff_num(3, 3, -2:2)
    integer :: a, b, i_step
    real(dp) :: delta
    real(dp) :: sigma, dsigma_dr, sigma_diff

    delta = 1d-3
    r = [1.02d0, -2.22d0, 0.15d0]
    sigma = 1.2d0
    dsigma_dr = -0.3d0
    T = T_erf_coulomb(r, sigma, deriv=.true.)
    T_diff_anl = T%dsigma(:, :)*dsigma_dr
    do i_step = -2, 2
        if (i_step == 0) cycle
        sigma_diff = sigma+i_step*delta*dsigma_dr
        T = T_erf_coulomb(r, sigma_diff, deriv=.false.)
        T_diff_num(:, :, i_step) = T%val
    end do
    forall (a = 1:3, b = 1:3)
        T_diff_num(a, b, 0) = diff5(T_diff_num(a, b, :), delta)
    end forall
    diff = (T_diff_num(:, :, 0)-T_diff_anl)/T_diff_num(:, :, 0)
    if (any(abs(diff) > 1d-10)) then
        call failed()
        call print_matrix('delta dTGG', diff)
    end if
end subroutine test_T_GG_deriv_impl

subroutine test_T_fermi_deriv_impl()
    real(dp) :: r(3)
    type(mat33) :: T
    real(dp) :: diff(3, 3)
    real(dp) :: T_diff_anl(3, 3)
    real(dp) :: T_diff_num(3, 3, -2:2)
    integer :: a, b, i_step
    real(dp) :: delta
    real(dp) :: rvdw, drvdw_dr, rvdw_diff

    delta = 1d-3
    r = [1.02d0, -2.22d0, 0.15d0]
    rvdw = 2.5d0
    drvdw_dr = -0.3d0
    T = damping_fermi(r, rvdw, 6.d0, .true.).prod.T_bare_v2(r, .true.)
    T_diff_anl = T%dvdw(:, :)*drvdw_dr
    do i_step = -2, 2
        if (i_step == 0) cycle
        rvdw_diff =rvdw+i_step*delta*drvdw_dr
        T = damping_fermi(r, rvdw_diff, 6.d0, .false.).prod.T_bare_v2(r, .false.)
        T_diff_num(:, :, i_step) = T%val
    end do
    forall (a = 1:3, b = 1:3)
        T_diff_num(a, b, 0) = diff5(T_diff_num(a, b, :), delta)
    end forall
    diff = (T_diff_num(:, :, 0)-T_diff_anl)/T_diff_num(:, :, 0)
    if (any(abs(diff) > 1d-10)) then
        call failed()
        call print_matrix('delta dTfermi', diff)
    end if
end subroutine test_T_fermi_deriv_impl

subroutine test_mbd_deriv_expl()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: forces(:, :)
    real(dp), allocatable :: diff(:, :)
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
    real(dp) :: ene(-2:2)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 1d-2
    n_atoms = 3
    allocate (coords(n_atoms, 3), source=0.d0)
    allocate (forces(n_atoms, 3))
    coords(1, 3) = 1.d0*ang
    coords(2, 1) = 4.d0*ang
    coords(3, 2) = 4.d0*ang
    sys%calc => calc
    sys%coords = coords
    sys%do_force = .true.
    damp%version = 'fermi,dip'
    damp%r_vdw%val = [3.55d0, 3.55d0, 3.55d0]
    damp%beta = 0.83
    alpha_0 = [11.d0, 11.d0, 11.d0]
    C6 = [65d0, 65d0, 65d0]
    ene(0) = get_single_mbd_energy(sys, vecn(alpha_0), vecn(C6), damp)
    sys%do_force = .false.
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -2, 2
                if (i_step == 0) cycle
                sys%coords = coords
                sys%coords(i_atom, i_xyz) = sys%coords(i_atom, i_xyz)+i_step*delta
                ene(i_step) = get_single_mbd_energy(sys, vecn(alpha_0), vecn(C6), damp)
            end do
            forces(i_atom, i_xyz) = diff5(ene, delta)
        end do
    end do
    diff = (forces-sys%work%forces)/forces
    if (any(abs(diff) > 1d-8)) then
        call failed()
        call print_matrix('delta forces', diff)
    end if
end subroutine test_mbd_deriv_expl

subroutine test_scs_deriv_expl()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: forces(:, :, :)
    real(dp), allocatable :: diff(:, :, :)
    real(dp), allocatable :: alpha_0(:)
    integer :: i_atom, n_atoms, i_xyz, i_step, j_atom
    type(vecn) :: alpha_scs(-2:2)

    delta = 1.5d-2
    n_atoms = 3
    allocate (coords(n_atoms, 3), source=0.d0)
    allocate (forces(n_atoms, n_atoms, 3))
    coords(1, 3) = 1.d0*ang
    coords(2, 1) = 4.d0*ang
    coords(3, 2) = 4.d0*ang
    sys%calc => calc
    sys%coords = coords
    sys%do_force = .true.
    damp%version = 'fermi,dip,gg'
    damp%r_vdw%val = [3.55d0, 3.55d0, 3.55d0]
    damp%beta = 0.83
    alpha_0 = [11.d0, 11.d0, 11.d0]
    alpha_scs(0) = run_scs(sys, vecn(alpha_0), damp)
    sys%do_force = .false.
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -2, 2
                if (i_step == 0) cycle
                sys%coords = coords
                sys%coords(i_atom, i_xyz) = sys%coords(i_atom, i_xyz)+i_step*delta
                alpha_scs(i_step) = run_scs(sys, vecn(alpha_0), damp)
            end do
            do j_atom = 1, n_atoms
                forces(j_atom, i_atom, i_xyz) = &
                    diff5([(alpha_scs(i_step)%val(j_atom), i_step = -2, 2)], delta)
            end do
        end do
    end do
    diff = (forces-alpha_scs(0)%dr)/forces
    if (any(abs(diff) > 1d-7)) then
        call failed()
        call print_matrix('diff x', diff(:, :, 1))
        call print_matrix('diff y', diff(:, :, 2))
        call print_matrix('diff z', diff(:, :, 3))
    end if
end subroutine test_scs_deriv_expl

subroutine test_mbd_deriv_impl_alpha()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: forces(:, :)
    real(dp), allocatable :: diff(:, :)
    type(vecn) :: alpha_0
    real(dp), allocatable :: alpha_0_diff(:)
    real(dp), allocatable :: C6(:)
    real(dp) :: ene(-2:2)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 1d-2
    n_atoms = 3
    allocate (coords(n_atoms, 3), source=0.d0)
    allocate (forces(n_atoms, 3))
    coords(2, 1) = 4.d0*ang
    coords(3, 2) = 4.d0*ang
    sys%calc => calc
    sys%coords = coords
    sys%do_force = .true.
    damp%version = 'fermi,dip'
    damp%r_vdw%val = [3.55d0, 3.55d0, 3.55d0]
    damp%beta = 0.83
    alpha_0%val= [11.d0, 11.d0, 11.d0]
    allocate (alpha_0%dr(n_atoms, n_atoms, 3), source=200d0)
    C6 = [65d0, 65d0, 65d0]
    ene(0) = get_single_mbd_energy(sys, alpha_0, vecn(C6), damp)
    sys%do_force = .false.
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -2, 2
                if (i_step == 0) cycle
                sys%coords = coords
                sys%coords(i_atom, i_xyz) = sys%coords(i_atom, i_xyz)+i_step*delta
                alpha_0_diff = alpha_0%val + alpha_0%dr(:, i_atom, i_xyz)*i_step*delta
                ene(i_step) = get_single_mbd_energy(sys, vecn(alpha_0_diff), vecn(C6), damp)
            end do
            forces(i_atom, i_xyz) = diff5(ene, delta)
        end do
    end do
    diff = (forces-sys%work%forces)/forces
    if (any(abs(diff) > 1d-7)) then
        call failed()
        call print_matrix('delta forces', diff)
    end if
end subroutine test_mbd_deriv_impl_alpha

subroutine test_mbd_deriv_impl_C6()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: forces(:, :)
    real(dp), allocatable :: diff(:, :)
    type(vecn) :: C6
    real(dp), allocatable :: C6_diff(:)
    real(dp), allocatable :: alpha_0(:)
    real(dp) :: ene(-2:2)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 1d-2
    n_atoms = 3
    allocate (coords(n_atoms, 3), source=0.d0)
    allocate (forces(n_atoms, 3))
    coords(2, 1) = 4.d0*ang
    coords(3, 2) = 4.d0*ang
    sys%calc => calc
    sys%coords = coords
    sys%do_force = .true.
    damp%version = 'fermi,dip'
    damp%r_vdw%val = [3.55d0, 3.55d0, 3.55d0]
    damp%beta = 0.83
    alpha_0 = [11.d0, 11.d0, 11.d0]
    C6%val = [65d0, 65d0, 65d0]
    allocate (C6%dr(n_atoms, n_atoms, 3), source=200d0)
    ene(0) = get_single_mbd_energy(sys, vecn(alpha_0), C6, damp)
    sys%do_force = .false.
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -2, 2
                if (i_step == 0) cycle
                sys%coords = coords
                sys%coords(i_atom, i_xyz) = sys%coords(i_atom, i_xyz)+i_step*delta
                C6_diff = C6%val + C6%dr(:, i_atom, i_xyz)*i_step*delta
                ene(i_step) = get_single_mbd_energy(sys, vecn(alpha_0), vecn(C6_diff), damp)
            end do
            forces(i_atom, i_xyz) = diff5(ene, delta)
        end do
    end do
    diff = (forces-sys%work%forces)/forces
    if (any(abs(diff) > 1d-9)) then
        call failed()
        call print_matrix('delta forces', diff)
    end if
end subroutine test_mbd_deriv_impl_C6

subroutine test_mbd_deriv_impl_vdw()
    real(dp) :: delta
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: forces(:, :)
    real(dp), allocatable :: diff(:, :)
    real(dp), allocatable :: rvdw(:)
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
    real(dp) :: ene(-2:2)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 1d-3
    n_atoms = 3
    allocate (coords(n_atoms, 3), source=0.d0)
    allocate (forces(n_atoms, 3))
    coords(2, 1) = 4.d0*ang
    coords(3, 2) = 4.d0*ang
    sys%calc => calc
    sys%coords = coords
    sys%do_force = .true.
    damp%version = 'fermi,dip'
    rvdw = [3.55d0, 3.55d0, 3.55d0]
    damp%r_vdw%val = rvdw
    allocate (damp%r_vdw%dr(n_atoms, n_atoms, 3), source=5d0)
    damp%beta = 0.83
    alpha_0 = [11.d0, 11.d0, 11.d0]
    C6 = [65d0, 65d0, 65d0]
    ene(0) = get_single_mbd_energy(sys, vecn(alpha_0), vecn(C6), damp)
    sys%do_force = .false.
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -2, 2
                if (i_step == 0) cycle
                sys%coords = coords
                sys%coords(i_atom, i_xyz) = sys%coords(i_atom, i_xyz)+i_step*delta
                damp%r_vdw%val = rvdw + damp%r_vdw%dr(:, i_atom, i_xyz)*i_step*delta
                ene(i_step) = get_single_mbd_energy(sys, vecn(alpha_0), vecn(C6), damp)
            end do
            forces(i_atom, i_xyz) = diff5(ene, delta)
        end do
    end do
    diff = (forces-sys%work%forces)/forces
    if (any(abs(diff) > 1d-9)) then
        call failed()
        call print_matrix('delta forces', diff)
    end if
end subroutine test_mbd_deriv_impl_vdw

end program
