! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_c_api

use iso_c_binding, only: c_ptr, c_int, c_double, c_f_pointer, c_loc, c_bool, &
    c_null_ptr, c_null_char, c_char
use mbd, only: mbd_system, mbd_calc, mbd_damping, get_mbd_energy, init_grid, &
    mbd_scs_energy, mbd_scs_energy, dipole_matrix, get_ts_energy, mbd_result, &
    mbd_gradients
use mbd_common, only: dp
use mbd_types, only: mat3n3n
use mbd_repulsion, only: full_coulomb, get_dipole_energy

implicit none

#ifdef WITH_SCALAPACK
logical(c_bool), bind(c) :: with_scalapack = .true.
#else
logical(c_bool), bind(c) :: with_scalapack = .false.
#endif

type, bind(c) :: mbd_calc_c
    integer(c_int) :: n_freq = 0
    type(c_ptr) :: omega_grid = c_null_ptr
    type(c_ptr) :: omega_grid_w = c_null_ptr
    type(c_ptr) :: mbd_calc_f = c_null_ptr
end type

type, bind(c) :: mbd_system_c
    type(c_ptr) :: mbd_system_f = c_null_ptr
end type

contains

type(c_ptr) function mbd_init_calc(n_freq) bind(c)
    integer(c_int), intent(in), value :: n_freq

    type(mbd_calc), pointer :: calc
    type(mbd_calc_c), pointer :: calc_c

    allocate (calc)
    calc%param%n_frequency_grid = n_freq
    call init_grid(calc)
    allocate (calc_c)
    calc_c%n_freq = ubound(calc%omega_grid, 1)
    calc_c%omega_grid = c_loc(calc%omega_grid)
    calc_c%omega_grid_w = c_loc(calc%omega_grid_w)
    calc_c%mbd_calc_f = c_loc(calc)
    mbd_init_calc = c_loc(calc_c)
end function mbd_init_calc

subroutine mbd_set_parallel(calc_cp) bind(c)
    type(c_ptr), value :: calc_cp

    type(mbd_calc), pointer :: calc

    calc => get_mbd_calc(calc_cp)
end subroutine mbd_set_parallel

subroutine mbd_destroy_calc(calc_cp) bind(c)
    type(c_ptr), value :: calc_cp

    type(mbd_calc), pointer :: calc
    type(mbd_calc_c), pointer :: calc_c

    call c_f_pointer(calc_cp, calc_c)
    call c_f_pointer(calc_c%mbd_calc_f, calc)
    deallocate (calc)
    deallocate (calc_c)
end subroutine mbd_destroy_calc

type(c_ptr) function mbd_init_system(calc_cp, n_atoms, coords, lattice, k_grid) bind(c)
    type(c_ptr), value :: calc_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: coords(3, n_atoms)
    real(c_double), intent(in), optional :: lattice(3, 3)
    integer(c_int), intent(in), optional :: k_grid(3)

    type(mbd_calc), pointer :: calc
    type(mbd_system), pointer :: sys
    type(mbd_system_c), pointer :: sys_c

    calc => get_mbd_calc(calc_cp)
    allocate (sys)
    sys%calc => calc
    sys%coords = coords
    call sys%blacs_grid%init()
    if (present(lattice)) then
        sys%periodic = .true.
        sys%lattice = lattice
    end if
    if (present(k_grid)) sys%k_grid = k_grid
    allocate (sys_c)
    sys_c%mbd_system_f = c_loc(sys)
    mbd_init_system = c_loc(sys_c)
end function mbd_init_system

subroutine mbd_destroy_system(sys_cp) bind(c)
    type(c_ptr), value :: sys_cp

    type(mbd_system_c), pointer :: sys_c
    type(mbd_system), pointer :: sys

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, sys)
    call sys%blacs_grid%destroy()
    deallocate (sys)
    deallocate (sys_c)
end subroutine mbd_destroy_system

type(c_ptr) function mbd_init_damping(n_atoms, version_c, r_vdw, sigma, beta, a) bind(c)
    integer(c_int), intent(in), value :: n_atoms
    character(kind=c_char), intent(in) :: version_c(*)
    real(c_double), intent(in), optional :: r_vdw(n_atoms)
    real(c_double), intent(in), optional :: sigma(n_atoms)
    real(c_double), intent(in), value :: beta
    real(c_double), intent(in), value :: a

    type(mbd_damping), pointer :: damping

    allocate (damping)
    damping%version = f_string(version_c)
    if (present(r_vdw)) damping%r_vdw = r_vdw
    if (present(sigma)) damping%sigma = sigma
    damping%beta = beta
    damping%a = a
    damping%ts_sr = beta
    damping%ts_d = a
    mbd_init_damping = c_loc(damping)
end function mbd_init_damping

subroutine mbd_destroy_damping(damping_p) bind(c)
    type(c_ptr), value :: damping_p

    type(mbd_damping), pointer :: damping

    call c_f_pointer(damping_p, damping)
    if (allocated(damping%r_vdw)) deallocate (damping%r_vdw)
    if (allocated(damping%sigma)) deallocate (damping%sigma)
    deallocate (damping)
end subroutine mbd_destroy_damping

function get_mbd_system(sys_cp)
    type(c_ptr), intent(in), value :: sys_cp
    type(mbd_system), pointer :: get_mbd_system

    type(mbd_system_c), pointer :: sys_c

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, get_mbd_system)
end function

function get_mbd_calc(calc_cp)
    type(c_ptr), intent(in), value :: calc_cp
    type(mbd_calc), pointer :: get_mbd_calc

    type(mbd_calc_c), pointer :: calc_c

    call c_f_pointer(calc_cp, calc_c)
    call c_f_pointer(calc_c%mbd_calc_f, get_mbd_calc)
end function

real(c_double) function calc_ts_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)

    type(mbd_system_c), pointer :: sys_c
    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, sys)
    call c_f_pointer(damping_p, damping)
    calc_ts_energy = get_ts_energy(sys, alpha_0, C6, damping)
end function calc_ts_energy

real(c_double) function calc_mbd_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)

    type(mbd_system_c), pointer :: sys_c
    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping
    type(mbd_result) :: res
    type(mbd_gradients) :: dene

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, sys)
    call c_f_pointer(damping_p, damping)
    if (present(gradients)) allocate (dene%dcoords(n_atoms, 3))
    res = get_mbd_energy(sys, alpha_0, C6, damping, dene)
    calc_mbd_energy = res%energy
    if (present(gradients)) gradients = transpose(dene%dcoords)
end function calc_mbd_energy

real(c_double) function calc_rpa_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)

    type(mbd_system), pointer :: sys
    type(mbd_system) :: sys2
    type(mbd_damping), pointer :: damping
    type(mbd_result) :: res
    type(mbd_gradients) :: dene

    sys => get_mbd_system(sys_cp)
    call c_f_pointer(damping_p, damping)
    sys2 = sys
    sys2%do_rpa = .true.
    res = get_mbd_energy(sys2, alpha_0, C6, damping, dene)
    calc_rpa_energy = res%energy
end function calc_rpa_energy

real(c_double) function calc_mbd_rsscs_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients, eigvals, eigvecs) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: eigvals(3*n_atoms)
    real(c_double), intent(out), optional :: eigvecs(3*n_atoms, 3*n_atoms)

    type(mbd_system_c), pointer :: sys_c
    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping
    type(mbd_result) :: res
    type(mbd_gradients) :: dene

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, sys)
    call c_f_pointer(damping_p, damping)
    if (present(gradients)) allocate (dene%dcoords(n_atoms, 3))
    sys%get_eigs = present(eigvals)
    sys%get_modes = present(eigvecs)
    res = mbd_scs_energy(sys, 'rsscs', alpha_0, C6, damping, dene)
    if (sys%has_exc()) return
    calc_mbd_rsscs_energy = res%energy
    if (present(gradients)) gradients = transpose(dene%dcoords)
    if (present(eigvals)) eigvals = res%mode_eigs
    if (present(eigvecs)) eigvecs = res%modes
end function calc_mbd_rsscs_energy

real(c_double) function calc_mbd_scs_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)

    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping
    type(mbd_result) :: res
    type(mbd_gradients) :: dene

    sys => get_mbd_system(sys_cp)
    call c_f_pointer(damping_p, damping)
    res = mbd_scs_energy(sys, 'scs', alpha_0, C6, damping, dene)
    calc_mbd_scs_energy = res%energy
end function calc_mbd_scs_energy

subroutine calc_dipole_matrix(sys_cp, damping_p, k_point, dipmat_p) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(in), optional :: k_point(3)
    type(c_ptr), intent(in), value :: dipmat_p

    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damp
    type(mat3n3n) :: dipmat
    real(dp), pointer :: dipmat_re(:, :)
    complex(dp), pointer :: dipmat_cplx(:, :)
    integer :: n_atoms

    sys => get_mbd_system(sys_cp)
    n_atoms = size(sys%coords, 2)
    call c_f_pointer(damping_p, damp)
    dipmat = dipole_matrix(sys, damp, .false., k_point)
    if (present(k_point)) then
        call c_f_pointer(dipmat_p, dipmat_cplx, [3*n_atoms, 3*n_atoms])
        dipmat_cplx = transpose(dipmat%cplx)
    else
        call c_f_pointer(dipmat_p, dipmat_re, [3*n_atoms, 3*n_atoms])
        dipmat_re = transpose(dipmat%re)
    end if
end subroutine calc_dipole_matrix

subroutine calc_full_coulomb(n, coords, C, w, w0, a0, rvdw0, alpha, beta, &
        version, dampswitch, ecoul, en, ee, nn) bind(c)
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: &
        coords(n, 3), C(3*n, 3*n), w(3*n), w0(n), a0(n), rvdw0(n)
    real(c_double), value, intent(in) :: alpha, beta, dampswitch
    character(c_char), intent(in) :: version(20)
    real(c_double), intent(out) :: ecoul, en, ee, nn

    call full_coulomb( &
        n, coords, C, w, w0, a0, rvdw0, &
        alpha, beta, f_string(version), dampswitch, ecoul, en, ee, nn &
    )
end subroutine calc_full_coulomb

real(c_double) function calc_get_dipole_energy(n, version, R, a0, w, w_t, r0, &
        beta, alpha, C) result(ene) bind(c)
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: &
        R(n, 3), C(3*n, 3*n), w_t(3*n), w(n), a0(n), r0(n)
    real(c_double), value, intent(in) :: alpha, beta
    character(c_char), intent(in) :: version(20)

    ene = get_dipole_energy( &
        f_string(version), R, a0, w, w_t,r0, beta, alpha, C &
    )
end function calc_get_dipole_energy

function f_string(str_c) result(str_f)
    character(kind=c_char), intent(in) :: str_c(*)
    character(len=:), allocatable :: str_f

    integer :: i
    i = 0
    do
        if (str_c(i+1) == c_null_char) exit
        i = i + 1
    end do
    allocate (character(len=i) :: str_f)
    str_f = transfer(str_c(1:i), str_f)
end function f_string

end module mbd_c_api
