! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_grad_test_cases

use mbd_constants
use mbd_damping, only: damping_t, damping_fermi
use mbd_dipole, only: dipole_matrix, T_bare, T_erf_coulomb, damping_grad, T_erfc
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_t, grad_matrix_re_t, grad_request_t, grad_scalar_t
use mbd_hamiltonian, only: get_mbd_hamiltonian_energy
use mbd_methods, only: get_mbd_scs_energy
use mbd_scs, only: run_scs
use mbd_utils, only: diff7, findval, tostr, result_t

implicit none

type(geom_t), allocatable :: geom
integer :: n_failed, rank

contains

logical function failed(diff, thre)
    real(dp), intent(in) :: diff, thre

    failed = .not. abs(diff) < thre
    if (rank == 0) write (6, '(A,G10.3,A,G10.3,A)') 'diff:', diff, ', threshold:', thre, ': '
    if (failed) n_failed = n_failed + 1
end function

subroutine print_matrix(label, A, prec)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: A(:, :)
    integer, optional, intent(in) :: prec

    integer :: m, n, i, j, prec_
    character(len=10) :: fm

    if (present(prec)) then
        prec_ = prec
    else
        prec_ = 3
    end if
    m = size(A, 1)
    n = size(A, 2)
    write (fm, '("(g",i2,".",i1,")")') prec_+8, prec_
    write (6, '(A,":")') label
    do i = 1, m
        do j = 1, n
            write (6, fm, advance="no") A(i, j)
        end do
        write (6, *)
    end do
end subroutine

subroutine test_T_bare_deriv()
    real(dp) :: r(3), r_diff(3), T(3, 3), diff(3, 3), T_diff_num(3, 3, -3:3), delta
    type(grad_matrix_re_t) :: dT
    integer :: a, b, c, i_step

    delta = 1d-2
    r = [1.12d0, -2.12d0, 0.12d0]
    T = T_bare(r, dT, .true.)
    diff = 0d0
    do c = 1, 3
        do i_step = -3, 3
            if (i_step == 0) cycle
            r_diff = r
            r_diff(c) = r_diff(c)+i_step*delta
            T_diff_num(:, :, i_step) = T_bare(r_diff)
        end do
        do concurrent (a = 1:3, b = 1:3)
            T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
        end do
        diff = max(diff, abs(T_diff_num(:, :, 0)-dT%dr(:, :, c))/T_diff_num(:, :, 0))
    end do
    if (failed(maxval(abs(diff)), 1d-10)) then
    end if
end subroutine

subroutine test_T_GG_deriv_expl()
    real(dp) :: r(3), r_diff(3), T(3, 3), diff(3, 3), T_diff_num(3, 3, -3:3), delta, sigma
    type(grad_matrix_re_t) :: dT
    integer :: a, b, c, i_step

    delta = 1d-2
    r = [1.02d0, -2.22d0, 0.15d0]
    sigma = 1.2d0
    T = T_erf_coulomb(r, sigma, dT, grad_request_t(dcoords=.true.))
    diff = 0d0
    do c = 1, 3
        do i_step = -3, 3
            if (i_step == 0) cycle
            r_diff = r
            r_diff(c) = r_diff(c)+i_step*delta
            T_diff_num(:, :, i_step) = T_erf_coulomb(r_diff, sigma)
        end do
        do concurrent (a = 1:3, b = 1:3)
            T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
        end do
        diff = max(diff, abs(T_diff_num(:, :, 0)-dT%dr(:, :, c))/T_diff_num(:, :, 0))
    end do
    if (failed(maxval(abs(diff)), 1d-10)) then
    end if
end subroutine

subroutine test_T_GG_deriv_impl()
    real(dp) :: r(3), T(3, 3), diff(3, 3), T_diff_num(3, 3, -3:3), delta, sigma, sigma_diff
    type(grad_matrix_re_t) :: dT
    integer :: a, b, i_step

    delta = 1d-3
    r = [1.02d0, -2.22d0, 0.15d0]
    sigma = 1.2d0
    T = T_erf_coulomb(r, sigma, dT, grad_request_t(dsigma=.true.))
    do i_step = -3, 3
        if (i_step == 0) cycle
        sigma_diff = sigma+i_step*delta
        T_diff_num(:, :, i_step) = T_erf_coulomb(r, sigma_diff)
    end do
    do concurrent (a = 1:3, b = 1:3)
        T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
    end do
    diff = (T_diff_num(:, :, 0)-dT%dsigma)/T_diff_num(:, :, 0)
    if (failed(maxval(abs(diff)), 1d-10)) then
        call print_matrix('delta dTGG', diff)
    end if
end subroutine

subroutine test_T_fermi_deriv_impl()
    real(dp) :: r(3), T(3, 3), T0(3, 3), &
        diff(3, 3), T_diff_num(3, 3, -3:3), delta, rvdw, rvdw_diff, f
    type(grad_matrix_re_t) :: dT, dT0
    type(grad_scalar_t) :: df
    integer :: a, b, i_step

    delta = 1d-3
    r = [1.02d0, -2.22d0, 0.15d0]
    rvdw = 2.5d0
    f = damping_fermi(r, rvdw, 6d0, df, grad_request_t(dr_vdw=.true.))
    T0 = T_bare(r)
    T = damping_grad(f, df, T0, dT0, dT, grad_request_t(dr_vdw=.true.))
    do i_step = -3, 3
        if (i_step == 0) cycle
        rvdw_diff =rvdw+i_step*delta
        T_diff_num(:, :, i_step) = damping_fermi(r, rvdw_diff, 6d0)*T_bare(r)
    end do
    do concurrent (a = 1:3, b = 1:3)
        T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
    end do
    diff = (T_diff_num(:, :, 0)-dT%dvdw)/T_diff_num(:, :, 0)
    if (failed(maxval(abs(diff)), 1d-10)) then
        call print_matrix('delta dTfermi', diff)
    end if
end subroutine

subroutine test_T_erfc_deriv_expl()
    real(dp) :: r(3), r_diff(3), T(3, 3), diff(3, 3), T_diff_num(3, 3, -3:3), delta, gamm
    type(grad_matrix_re_t) :: dT
    integer :: a, b, c, i_step

    delta = 1d-2
    r = [1.02d0, -2.22d0, 0.15d0]
    gamm = 1.2d0
    T = T_erfc(r, gamm, dT, grad_request_t(dcoords=.true.))
    diff = 0d0
    do c = 1, 3
        do i_step = -3, 3
            if (i_step == 0) cycle
            r_diff = r
            r_diff(c) = r_diff(c)+i_step*delta
            T_diff_num(:, :, i_step) = T_erfc(r_diff, gamm)
        end do
        do concurrent (a = 1:3, b = 1:3)
            T_diff_num(a, b, 0) = diff7(T_diff_num(a, b, :), delta)
        end do
        diff = max(diff, abs(T_diff_num(:, :, 0)-dT%dr(:, :, c))/T_diff_num(:, :, 0))
    end do
    if (failed(maxval(abs(diff)), 1d-9)) then
    end if
end subroutine

subroutine test_mbd_deriv_expl()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: gradients(:, :)
    real(dp), allocatable :: diff(:, :)
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: omega(:)
    type(result_t) :: res(-3:3)
    real(dp), allocatable :: gradients_anl(:, :)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 0.01d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms, 3))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    omega = [.7d0, .65d0, .75d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, grad_request_t(dcoords=.true.))
    gradients_anl = res(0)%dE%dcoords
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%coords = coords
                geom%coords(i_xyz, i_atom) = geom%coords(i_xyz, i_atom)+i_step*delta
                res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, grad_request_t())
            end do
            gradients(i_atom, i_xyz) = diff7(res%energy, delta)
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', diff)
    end if
end subroutine

subroutine test_mbd_ewald_deriv_expl()
    real(dp) :: delta, k_point(3)
    type(damping_t) :: damp
    real(dp), allocatable :: &
        coords(:, :), gradients(:, :), diff(:, :), alpha_0(:), omega(:), &
        gradients_anl(:, :)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 0.01d0
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms, 3))
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    k_point = [0.4d0, 0d0, 0d0]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    omega = [.7d0, .65d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(dcoords=.true.), k_point)
    gradients_anl = res(0)%dE%dcoords
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%coords = coords
                geom%coords(i_xyz, i_atom) = geom%coords(i_xyz, i_atom)+i_step*delta
                res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
                    grad_request_t(), k_point)
            end do
            gradients(i_atom, i_xyz) = diff7(res%energy, delta)
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', diff)
        call print_matrix('anl', gradients_anl)
        call print_matrix('num', gradients)
    end if
end subroutine

subroutine test_mbd_ewald_deriv_impl_q()
    real(dp) :: delta, k_point(3), k_point_diff(3)
    type(damping_t) :: damp
    real(dp), allocatable :: &
        gradients(:), diff(:), alpha_0(:), omega(:), gradients_anl(:)
    type(result_t) :: res(-3:3)
    integer :: n_atoms, i_xyz, i_step

    delta = 0.01d0
    n_atoms = 2
    allocate (geom%coords(3, n_atoms), source=0d0)
    allocate (gradients(3))
    geom%coords(3, 1) = 1d0
    geom%coords(1, 2) = 1d0
    geom%coords(2, 2) = 4d0
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    k_point = [0.4d0, 0d0, 0d0]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    omega = [.7d0, .65d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(dq=.true.), k_point)
    gradients_anl = res(0)%dE%dq
    do i_xyz = 1, 3
        do i_step = -3, 3
            if (i_step == 0) cycle
            k_point_diff = k_point
            k_point_diff(i_xyz) = k_point_diff(i_xyz)+i_step*delta
            res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
                grad_request_t(), k_point_diff)
        end do
        gradients(i_xyz) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        print *, 'delta gradients', diff
        print *, 'anl', gradients_anl
        print *, 'num', gradients
    end if
end subroutine

subroutine test_mbd_ewald_deriv_stress()
    real(dp) :: delta, k_point(3)
    type(damping_t) :: damp
    real(dp), allocatable :: &
        lattice(:, :), gradients(:, :), diff(:, :), alpha_0(:), omega(:), &
        gradients_anl(:, :)
    type(result_t) :: res(-3:3)
    integer :: i_latt, n_atoms, i_xyz, i_step

    delta = 0.01d0
    n_atoms = 2
    allocate (geom%coords(3, n_atoms), source=0d0)
    allocate (gradients(3, 3))
    geom%coords(3, 1) = 1d0
    geom%coords(1, 2) = 1d0
    geom%coords(2, 2) = 4d0
    lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    geom%lattice = lattice
    k_point = [0.4d0, 0d0, 0d0]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    omega = [.7d0, .65d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(dlattice=.true.), k_point)
    gradients_anl = res(0)%dE%dlattice
    do i_latt = 1, 3
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%lattice = lattice
                geom%lattice(i_xyz, i_latt) = geom%lattice(i_xyz, i_latt)+i_step*delta
                res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
                    grad_request_t(), k_point)
            end do
            gradients(i_latt, i_xyz) = diff7(res%energy, delta)
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-5)) then
        call print_matrix('delta gradients', diff)
        call print_matrix('anl', gradients_anl)
        call print_matrix('num', gradients)
    end if
end subroutine

subroutine test_scs_deriv_expl()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: gradients(:, :, :), gradients_anl(:, :, :)
    real(dp), allocatable :: diff(:, :, :)
    real(dp), allocatable :: alpha_0(:)
    integer :: i_atom, n_atoms, i_xyz, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    real(dp), allocatable :: alpha_scs(:, :)
    type(grad_t), allocatable :: dalpha_scs(:)

    delta = 0.05d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    my_nratoms = size(geom%idx%i_atom)
    my_ncatoms = size(geom%idx%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms, 3))
    allocate (gradients_anl(my_nratoms, my_ncatoms, 3))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    alpha_scs(:, 0) = &
        run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t(dcoords=.true.))
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :, :) = dalpha_scs(my_i_atom)%dcoords
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(geom%idx%j_atom, j_atom)
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%coords = coords
                geom%coords(i_xyz, j_atom) = geom%coords(i_xyz, j_atom) + &
                    i_step*delta
                alpha_scs(:, i_step) = &
                    run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t())
            end do
            if (my_j_atom > 0) then
                do my_i_atom = 1, my_nratoms
                    i_atom = geom%idx%i_atom(my_i_atom)
                    gradients(my_i_atom, my_j_atom, i_xyz) = &
                        diff7(alpha_scs(i_atom, :), delta)
                end do
            end if
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-5)) then
        call print_matrix('diff x', diff(:, :, 1))
        call print_matrix('diff y', diff(:, :, 2))
        call print_matrix('diff z', diff(:, :, 3))
    end if
end subroutine

subroutine test_scs_ewald_deriv_expl()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: gradients(:, :, :), gradients_anl(:, :, :)
    real(dp), allocatable :: diff(:, :, :)
    real(dp), allocatable :: alpha_0(:)
    integer :: i_atom, n_atoms, i_xyz, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    real(dp), allocatable :: alpha_scs(:, :)
    type(grad_t), allocatable :: dalpha_scs(:)

    delta = 0.05d0
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    call geom%init()
    my_nratoms = size(geom%idx%i_atom)
    my_ncatoms = size(geom%idx%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms, 3))
    allocate (gradients_anl(my_nratoms, my_ncatoms, 3))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    alpha_scs(:, 0) = &
        run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t(dcoords=.true.))
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :, :) = dalpha_scs(my_i_atom)%dcoords
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(geom%idx%j_atom, j_atom)
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%coords = coords
                geom%coords(i_xyz, j_atom) = geom%coords(i_xyz, j_atom) + &
                    i_step*delta
                alpha_scs(:, i_step) = &
                    run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t())
            end do
            if (my_j_atom > 0) then
                do my_i_atom = 1, my_nratoms
                    i_atom = geom%idx%i_atom(my_i_atom)
                    gradients(my_i_atom, my_j_atom, i_xyz) = &
                        diff7(alpha_scs(i_atom, :), delta)
                end do
            end if
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-5)) then
        call print_matrix('diff x', diff(:, :, 1))
        call print_matrix('diff y', diff(:, :, 2))
        call print_matrix('diff z', diff(:, :, 3))
    end if
end subroutine

subroutine test_scs_ewald_deriv_stress()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: gradients(:, :, :), gradients_anl(:, :, :), &
        diff(:, :, :), alpha_0(:), lattice(:, :)
    integer :: i_atom, n_atoms, i_xyz, i_step, my_i_atom, my_nratoms, &
        my_ncatoms, i_latt
    real(dp), allocatable :: alpha_scs(:, :)
    type(grad_t), allocatable :: dalpha_scs(:)

    delta = 0.05d0
    n_atoms = 2
    allocate (geom%coords(3, n_atoms), source=0d0)
    geom%coords(3, 1) = 1d0
    geom%coords(1, 2) = 1d0
    geom%coords(2, 2) = 4d0
    lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    geom%lattice = lattice
    call geom%init()
    my_nratoms = size(geom%idx%i_atom)
    my_ncatoms = size(geom%idx%j_atom)
    allocate (gradients(my_nratoms, 3, 3))
    allocate (gradients_anl(my_nratoms, 3, 3))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    alpha_scs(:, 0) = &
        run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t(dlattice=.true.))
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :, :) = dalpha_scs(my_i_atom)%dlattice
    end do
    do i_latt = 1, 3
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%lattice = lattice
                geom%lattice(i_xyz, i_latt) = geom%lattice(i_xyz, i_latt)+i_step*delta
                alpha_scs(:, i_step) = &
                    run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t())
            end do
            do my_i_atom = 1, my_nratoms
                i_atom = geom%idx%i_atom(my_i_atom)
                gradients(my_i_atom, i_latt, i_xyz) = diff7(alpha_scs(i_atom, :), delta)
            end do
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-4)) then
        call print_matrix('diff x', diff(:, :, 1))
        call print_matrix('diff y', diff(:, :, 2))
        call print_matrix('diff z', diff(:, :, 3))
    end if
end subroutine

subroutine test_scs_deriv_impl_alpha
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:, :), &
        gradients_anl(:, :), diff(:, :), alpha_0(:), alpha_0_diff(:), &
        alpha_scs(:, :)
    integer :: i_atom, n_atoms, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    type(grad_t), allocatable :: dalpha_scs(:)

    delta = 0.1d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    my_nratoms = size(geom%idx%i_atom)
    my_ncatoms = size(geom%idx%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms))
    allocate (gradients_anl(my_nratoms, my_ncatoms))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    alpha_scs(:, 0) = &
        run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t(dalpha=.true.))
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :) = dalpha_scs(my_i_atom)%dalpha
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(geom%idx%j_atom, j_atom)
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(j_atom) = alpha_0_diff(j_atom) + i_step*delta
            alpha_scs(:, i_step) = &
                run_scs(geom, alpha_0_diff, damp, dalpha_scs, grad_request_t())
        end do
        if (my_j_atom > 0) then
            do my_i_atom = 1, my_nratoms
                i_atom = geom%idx%i_atom(my_i_atom)
                gradients(my_i_atom, my_j_atom) = diff7(alpha_scs(i_atom, :), delta)
            end do
    end if
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        call print_matrix('diff', diff)
    end if
end subroutine

subroutine test_scs_ewald_deriv_impl_alpha
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:, :), &
        gradients_anl(:, :), diff(:, :), alpha_0(:), alpha_0_diff(:), &
        alpha_scs(:, :)
    integer :: i_atom, n_atoms, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    type(grad_t), allocatable :: dalpha_scs(:)

    delta = 0.1d0
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    call geom%init()
    my_nratoms = size(geom%idx%i_atom)
    my_ncatoms = size(geom%idx%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms))
    allocate (gradients_anl(my_nratoms, my_ncatoms))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    alpha_scs(:, 0) = &
        run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t(dalpha=.true.))
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :) = dalpha_scs(my_i_atom)%dalpha
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(geom%idx%j_atom, j_atom)
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(j_atom) = alpha_0_diff(j_atom) + i_step*delta
            alpha_scs(:, i_step) = &
                run_scs(geom, alpha_0_diff, damp, dalpha_scs, grad_request_t())
        end do
        if (my_j_atom > 0) then
            do my_i_atom = 1, my_nratoms
                i_atom = geom%idx%i_atom(my_i_atom)
                gradients(my_i_atom, my_j_atom) = diff7(alpha_scs(i_atom, :), delta)
            end do
    end if
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        call print_matrix('diff', diff)
    end if
end subroutine

subroutine test_scs_deriv_impl_vdw
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:, :), &
        gradients_anl(:, :), diff(:, :), alpha_0(:), alpha_scs(:, :), rvdw(:)
    integer :: i_atom, n_atoms, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    type(grad_t), allocatable :: dalpha_scs(:)

    delta = 0.1d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    my_nratoms = size(geom%idx%i_atom)
    my_ncatoms = size(geom%idx%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms))
    allocate (gradients_anl(my_nratoms, my_ncatoms))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    rvdw = [3.55d0, 3.5d0, 3.56d0]
    damp%r_vdw = rvdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    alpha_scs(:, 0) = &
        run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t(dr_vdw=.true.))
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :) = dalpha_scs(my_i_atom)%dr_vdw
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(geom%idx%j_atom, j_atom)
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = rvdw
            damp%r_vdw(j_atom) = damp%r_vdw(j_atom) + i_step*delta
            alpha_scs(:, i_step) = &
                run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t())
        end do
        if (my_j_atom > 0) then
            do my_i_atom = 1, my_nratoms
                i_atom = geom%idx%i_atom(my_i_atom)
                gradients(my_i_atom, my_j_atom) = diff7(alpha_scs(i_atom, :), delta)
            end do
        end if
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        call print_matrix('diff', diff(:, :))
    end if
end subroutine

subroutine test_scs_ewald_deriv_impl_vdw
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:, :), &
        gradients_anl(:, :), diff(:, :), alpha_0(:), alpha_scs(:, :), rvdw(:)
    integer :: i_atom, n_atoms, i_step, j_atom, my_i_atom, my_nratoms, &
        my_ncatoms, my_j_atom
    type(grad_t), allocatable :: dalpha_scs(:)

    delta = 0.1d0
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    call geom%init()
    my_nratoms = size(geom%idx%i_atom)
    my_ncatoms = size(geom%idx%j_atom)
    allocate (gradients(my_nratoms, my_ncatoms))
    allocate (gradients_anl(my_nratoms, my_ncatoms))
    allocate (alpha_scs(n_atoms, -3:3), dalpha_scs(my_nratoms))
    damp%version = 'fermi,dip,gg'
    rvdw = [3.55d0, 3.5d0]
    damp%r_vdw = rvdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    alpha_scs(:, 0) = &
        run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t(dr_vdw=.true.))
    do my_i_atom = 1, my_nratoms
        gradients_anl(my_i_atom, :) = dalpha_scs(my_i_atom)%dr_vdw
    end do
    do j_atom = 1, n_atoms
        my_j_atom = findval(geom%idx%j_atom, j_atom)
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = rvdw
            damp%r_vdw(j_atom) = damp%r_vdw(j_atom) + i_step*delta
            alpha_scs(:, i_step) = &
                run_scs(geom, alpha_0, damp, dalpha_scs, grad_request_t())
        end do
        if (my_j_atom > 0) then
            do my_i_atom = 1, my_nratoms
                i_atom = geom%idx%i_atom(my_i_atom)
                gradients(my_i_atom, my_j_atom) = diff7(alpha_scs(i_atom, :), delta)
            end do
        end if
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-5)) then
        call print_matrix('diff', diff(:, :))
    end if
end subroutine

subroutine test_mbd_deriv_impl_alpha()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), alpha_0_diff(:), omega(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 0.1d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    omega = [.7d0, .65d0, .75d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(dalpha=.true.))
    gradients_anl = res(0)%dE%dalpha
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(i_atom) = alpha_0_diff(i_atom) + i_step*delta
            res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0_diff, omega, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        call print_matrix('diff', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_ewald_deriv_impl_alpha()
    real(dp) :: delta, k_point(3)
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), alpha_0_diff(:), omega(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 0.1d0
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    k_point = [0.4d0, 0d0, 0d0]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    omega = [.7d0, .65d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(dalpha=.true.), k_point)
    gradients_anl = res(0)%dE%dalpha
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(i_atom) = alpha_0_diff(i_atom) + i_step*delta
            res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0_diff, omega, damp, &
                grad_request_t(), k_point)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-6)) then
        call print_matrix('diff', reshape(diff, [n_atoms, 1]))
        call print_matrix('anl', reshape(gradients_anl, [n_atoms, 1]))
        call print_matrix('num', reshape(gradients, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_deriv_impl_omega()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), omega_diff(:), omega(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 0.03d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    omega = [.7d0, .65d0, .75d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(domega=.true.))
    gradients_anl = res(0)%dE%domega
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            omega_diff = omega
            omega_diff(i_atom) = omega_diff(i_atom) + i_step*delta
            res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0, omega_diff, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 2d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_ewald_deriv_impl_omega()
    real(dp) :: delta, k_point(3)
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), omega_diff(:), omega(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 0.03d0
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    k_point = [0.4d0, 0d0, 0d0]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    omega = [.7d0, .65d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(domega=.true.), k_point)
    gradients_anl = res(0)%dE%domega
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            omega_diff = omega
            omega_diff(i_atom) = omega_diff(i_atom) + i_step*delta
            res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0, omega_diff, damp, &
                grad_request_t(), k_point)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 2d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_deriv_impl_vdw()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), omega(:), r_vdw(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 1d-3
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    damp%version = 'fermi,dip'
    r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%r_vdw = r_vdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    omega = [.7d0, .65d0, .75d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(dr_vdw=.true.))
    gradients_anl = res(0)%dE%dr_vdw
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = r_vdw
            damp%r_vdw(i_atom) = damp%r_vdw(i_atom) + i_step*delta
            res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_ewald_deriv_impl_vdw()
    real(dp) :: delta, k_point(3)
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), omega(:), r_vdw(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 1d-3
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    k_point = [0.4d0, 0d0, 0d0]
    call geom%init()
    damp%version = 'fermi,dip'
    r_vdw = [3.55d0, 3.5d0]
    damp%r_vdw = r_vdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    omega = [.7d0, .65d0]
    res(0) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(dr_vdw=.true.), k_point)
    gradients_anl = res(0)%dE%dr_vdw
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = r_vdw
            damp%r_vdw(i_atom) = damp%r_vdw(i_atom) + i_step*delta
            res(i_step) = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
                grad_request_t(), k_point)
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
        call print_matrix('gradients anl', reshape(gradients_anl, [n_atoms, 1]))
        call print_matrix('gradients num', reshape(gradients, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_rsscs_deriv_expl()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: gradients(:, :), gradients_anl(:, :)
    real(dp), allocatable :: diff(:, :)
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 0.01d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms, 3))
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    coords(1, 3) = 1d0
    geom%coords = coords
    call geom%init()
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dcoords=.true.))
    gradients_anl = res(0)%dE%dcoords
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%coords = coords
                geom%coords(i_xyz, i_atom) = geom%coords(i_xyz, i_atom) + &
                    i_step*delta
                res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
                    grad_request_t())
            end do
            gradients(i_atom, i_xyz) = diff7(res%energy, delta)
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', diff)
    end if
end subroutine

subroutine test_mbd_rsscs_ewald_deriv_expl()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: gradients(:, :), gradients_anl(:, :)
    real(dp), allocatable :: diff(:, :)
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_xyz, i_step

    delta = 0.01d0
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms, 3))
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    geom%k_grid = [2, 2, 2]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    C6 = [65d0, 60d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dcoords=.true.))
    gradients_anl = res(0)%dE%dcoords
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%coords = coords
                geom%coords(i_xyz, i_atom) = geom%coords(i_xyz, i_atom) + &
                    i_step*delta
                res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
                    grad_request_t())
            end do
            gradients(i_atom, i_xyz) = diff7(res%energy, delta)
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', diff)
    end if
end subroutine

subroutine test_mbd_rsscs_ewald_deriv_stress()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: gradients(:, :), gradients_anl(:, :), &
        diff(:, :), alpha_0(:), C6(:), lattice(:, :)
    type(result_t) :: res(-3:3)
    integer :: n_atoms, i_xyz, i_step, i_latt

    delta = 0.01d0
    n_atoms = 2
    allocate (geom%coords(3, n_atoms), source=0d0)
    allocate (gradients(3, 3))
    geom%coords(3, 1) = 1d0
    geom%coords(1, 2) = 1d0
    geom%coords(2, 2) = 4d0
    lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    geom%lattice = lattice
    geom%k_grid = [2, 2, 2]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    C6 = [65d0, 60d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dlattice=.true.))
    gradients_anl = res(0)%dE%dlattice
    do i_latt = 1, 3
        do i_xyz = 1, 3
            do i_step = -3, 3
                if (i_step == 0) cycle
                geom%lattice = lattice
                geom%lattice(i_xyz, i_latt) = geom%lattice(i_xyz, i_latt)+i_step*delta
                res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
                    grad_request_t())
            end do
            gradients(i_latt, i_xyz) = diff7(res%energy, delta)
        end do
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 3d-6)) then
        call print_matrix('delta gradients', diff)
        call print_matrix('gradients anl', gradients_anl)
        call print_matrix('gradients num', gradients)
    end if
end subroutine

subroutine test_mbd_rsscs_deriv_impl_alpha()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), alpha_0_diff(:), C6(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 3d-2
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dalpha=.true.))
    gradients_anl = res(0)%dE%dalpha
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(i_atom) = alpha_0_diff(i_atom) + i_step*delta
            res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0_diff, C6, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-7)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_rsscs_deriv_impl_C6()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), C6_diff(:), C6(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 0.01d0
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(1, 3) = 1d0
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dC6=.true.))
    gradients_anl = res(0)%dE%dC6
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            C6_diff = C6
            C6_diff(i_atom) = C6_diff(i_atom) + i_step*delta
            res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6_diff, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 5d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_rsscs_deriv_impl_vdw()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), C6(:), r_vdw(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 1d-2
    n_atoms = 3
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(2, 1) = 4d0
    coords(3, 2) = 4d0
    geom%coords = coords
    call geom%init()
    damp%version = 'fermi,dip'
    r_vdw = [3.55d0, 3.5d0, 3.56d0]
    damp%r_vdw = r_vdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0, 12d0]
    C6 = [65d0, 60d0, 70d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dr_vdw=.true.))
    gradients_anl = res(0)%dE%dr_vdw
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = r_vdw
            damp%r_vdw(i_atom) = damp%r_vdw(i_atom) + i_step*delta
            res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_rsscs_ewald_deriv_impl_alpha()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), alpha_0_diff(:), C6(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 3d-2
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    geom%k_grid = [2, 2, 2]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    C6 = [65d0, 60d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dalpha=.true.))
    gradients_anl = res(0)%dE%dalpha
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            alpha_0_diff = alpha_0
            alpha_0_diff(i_atom) = alpha_0_diff(i_atom) + i_step*delta
            res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0_diff, C6, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-7)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_rsscs_ewald_deriv_impl_C6()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), C6_diff(:), C6(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 0.01d0
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    geom%k_grid = [2, 2, 2]
    call geom%init()
    damp%version = 'fermi,dip'
    damp%r_vdw = [3.55d0, 3.5d0]
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    C6 = [65d0, 60d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dC6=.true.))
    gradients_anl = res(0)%dE%dC6
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            C6_diff = C6
            C6_diff(i_atom) = C6_diff(i_atom) + i_step*delta
            res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6_diff, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 5d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

subroutine test_mbd_rsscs_ewald_deriv_impl_vdw()
    real(dp) :: delta
    type(damping_t) :: damp
    real(dp), allocatable :: coords(:, :), gradients(:), &
        gradients_anl(:), diff(:), alpha_0(:), C6(:), r_vdw(:)
    type(result_t) :: res(-3:3)
    integer :: i_atom, n_atoms, i_step

    delta = 1d-2
    n_atoms = 2
    allocate (coords(3, n_atoms), source=0d0)
    allocate (gradients(n_atoms))
    coords(3, 1) = 1d0
    coords(1, 2) = 1d0
    coords(2, 2) = 4d0
    geom%coords = coords
    geom%lattice = reshape([6d0, 1d0, 0d0, -1d0, 9d0, 1d0, 0d0, 1d0, 7d0], [3, 3])
    geom%k_grid = [2, 2, 2]
    call geom%init()
    damp%version = 'fermi,dip'
    r_vdw = [3.55d0, 3.5d0]
    damp%r_vdw = r_vdw
    damp%beta = 0.83d0
    alpha_0 = [11d0, 10d0]
    C6 = [65d0, 60d0]
    res(0) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dr_vdw=.true.))
    gradients_anl = res(0)%dE%dr_vdw
    do i_atom = 1, n_atoms
        do i_step = -3, 3
            if (i_step == 0) cycle
            damp%r_vdw = r_vdw
            damp%r_vdw(i_atom) = damp%r_vdw(i_atom) + i_step*delta
            res(i_step) = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
                grad_request_t())
        end do
        gradients(i_atom) = diff7(res%energy, delta)
    end do
    call geom%destroy()
    diff = (gradients-gradients_anl)/gradients_anl
    if (failed(maxval(abs(diff)), 1d-8)) then
        call print_matrix('delta gradients', reshape(diff, [n_atoms, 1]))
    end if
end subroutine

end module
