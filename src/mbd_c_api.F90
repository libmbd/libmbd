! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_c_api

use iso_c_binding
use mbd_constants
use mbd_system_type, only: mbd_system, mbd_calc
use mbd_core, only: mbd_energy, mbd_scs_energy, mbd_scs_energy, mbd_result
use mbd_dipole, only: dipole_matrix
use mbd_damping_type, only: mbd_damping
use mbd_gradients_type, only: mbd_gradients, mbd_grad => mbd_grad_switch
use mbd_ts, only: ts_energy
use mbd_matrix_type, only: mbd_matrix_real, mbd_matrix_complex
use mbd_coulomb, only: dipole_energy, coulomb_energy

implicit none

#ifdef WITH_MPI
logical(c_bool), bind(c) :: cmbd_with_mpi = .true.
#else
logical(c_bool), bind(c) :: cmbd_with_mpi = .false.
#endif
#ifdef WITH_SCALAPACK
logical(c_bool), bind(c) :: cmbd_with_scalapack = .true.
#else
logical(c_bool), bind(c) :: cmbd_with_scalapack = .false.
#endif

type, bind(c) :: cmbd_calc
    integer(c_int) :: n_freq = 0
    type(c_ptr) :: omega_grid = c_null_ptr
    type(c_ptr) :: omega_grid_w = c_null_ptr
    type(c_ptr) :: mbd_calc_f = c_null_ptr
end type

type, bind(c) :: cmbd_system
    type(c_ptr) :: mbd_system_f = c_null_ptr
end type

contains

type(c_ptr) function cmbd_init_calc(n_freq) bind(c)
    integer(c_int), intent(in), value :: n_freq

    type(mbd_calc), pointer :: calc
    type(cmbd_calc), pointer :: calc_c

    allocate (calc)
    calc%param%n_frequency_grid = n_freq
    call calc%init_grid()
    allocate (calc_c)
    calc_c%n_freq = ubound(calc%omega_grid, 1)
    calc_c%omega_grid = c_loc(calc%omega_grid)
    calc_c%omega_grid_w = c_loc(calc%omega_grid_w)
    calc_c%mbd_calc_f = c_loc(calc)
    cmbd_init_calc = c_loc(calc_c)
end function cmbd_init_calc

subroutine cmbd_destroy_calc(calc_cp) bind(c)
    type(c_ptr), value :: calc_cp

    type(mbd_calc), pointer :: calc
    type(cmbd_calc), pointer :: calc_c

    call c_f_pointer(calc_cp, calc_c)
    call c_f_pointer(calc_c%mbd_calc_f, calc)
    deallocate (calc)
    deallocate (calc_c)
end subroutine cmbd_destroy_calc

subroutine cmbd_get_exception(calc_cp, code, origin, msg) bind(c)
    type(c_ptr), value :: calc_cp
    integer(c_int), intent(out) :: code
    character(kind=c_char), intent(out) :: origin(50), msg(150)

    type(mbd_calc), pointer :: calc

    calc => get_mbd_calc(calc_cp)
    code = calc%exc%code
    call f_c_string(calc%exc%origin, origin)
    call f_c_string(calc%exc%msg, msg)
    calc%exc%code = 0
    calc%exc%origin = ''
    calc%exc%msg = ''
end subroutine

type(c_ptr) function cmbd_init_system( &
        calc_cp, n_atoms, coords, lattice, k_grid) bind(c)
    type(c_ptr), value :: calc_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: coords(3, n_atoms)
    real(c_double), intent(in), optional :: lattice(3, 3)
    integer(c_int), intent(in), optional :: k_grid(3)

    type(mbd_calc), pointer :: calc
    type(mbd_system), pointer :: sys
    type(cmbd_system), pointer :: sys_c

    calc => get_mbd_calc(calc_cp)
    allocate (sys)
    sys%coords = coords
    call sys%init(calc)
    if (present(lattice)) sys%lattice = lattice
    if (present(k_grid)) sys%k_grid = k_grid
    allocate (sys_c)
    sys_c%mbd_system_f = c_loc(sys)
    cmbd_init_system = c_loc(sys_c)
end function cmbd_init_system

subroutine cmbd_destroy_system(sys_cp) bind(c)
    type(c_ptr), value :: sys_cp

    type(cmbd_system), pointer :: sys_c
    type(mbd_system), pointer :: sys

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, sys)
    call sys%destroy()
    deallocate (sys)
    deallocate (sys_c)
end subroutine cmbd_destroy_system

type(c_ptr) function cmbd_init_damping(n_atoms, version_c, r_vdw, sigma, beta, a) bind(c)
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
    cmbd_init_damping = c_loc(damping)
end function cmbd_init_damping

subroutine cmbd_destroy_damping(damping_p) bind(c)
    type(c_ptr), value :: damping_p

    type(mbd_damping), pointer :: damping

    call c_f_pointer(damping_p, damping)
    if (allocated(damping%r_vdw)) deallocate (damping%r_vdw)
    if (allocated(damping%sigma)) deallocate (damping%sigma)
    deallocate (damping)
end subroutine cmbd_destroy_damping

function get_mbd_system(sys_cp)
    type(c_ptr), intent(in), value :: sys_cp
    type(mbd_system), pointer :: get_mbd_system

    type(cmbd_system), pointer :: sys_c

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, get_mbd_system)
end function

function get_mbd_calc(calc_cp)
    type(c_ptr), intent(in), value :: calc_cp
    type(mbd_calc), pointer :: get_mbd_calc

    type(cmbd_calc), pointer :: calc_c

    call c_f_pointer(calc_cp, calc_c)
    call c_f_pointer(calc_c%mbd_calc_f, get_mbd_calc)
end function

real(c_double) function cmbd_ts_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)

    type(cmbd_system), pointer :: sys_c
    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, sys)
    call c_f_pointer(damping_p, damping)
    cmbd_ts_energy = ts_energy(sys, alpha_0, C6, damping)
end function cmbd_ts_energy

real(c_double) function cmbd_mbd_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)

    type(cmbd_system), pointer :: sys_c
    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping
    type(mbd_result) :: res
    type(mbd_gradients) :: dene

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, sys)
    call c_f_pointer(damping_p, damping)
    res = mbd_energy( &
        sys, alpha_0, C6, damping, dene, mbd_grad(dcoords=present(gradients)) &
    )
    cmbd_mbd_energy = res%energy
    if (present(gradients)) gradients = transpose(dene%dcoords)
end function cmbd_mbd_energy

real(c_double) function cmbd_rpa_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients) bind(c)
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
    sys%do_rpa = .true.
    res = mbd_energy(sys, alpha_0, C6, damping, dene, mbd_grad())
    sys%do_rpa = .false.
    cmbd_rpa_energy = res%energy
end function cmbd_rpa_energy

real(c_double) function cmbd_mbd_rsscs_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients, eigvals, eigvecs) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: eigvals(3*n_atoms)
    real(c_double), intent(out), optional :: eigvecs(3*n_atoms, 3*n_atoms)

    type(cmbd_system), pointer :: sys_c
    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping
    type(mbd_result) :: res
    type(mbd_gradients) :: dene

    call c_f_pointer(sys_cp, sys_c)
    call c_f_pointer(sys_c%mbd_system_f, sys)
    call c_f_pointer(damping_p, damping)
    sys%get_eigs = present(eigvals)
    sys%get_modes = present(eigvecs)
    res = mbd_scs_energy(sys, 'rsscs', alpha_0, C6, damping, &
        dene, mbd_grad(dcoords=present(gradients)))
    if (sys%has_exc()) return
    cmbd_mbd_rsscs_energy = res%energy
    if (present(gradients)) gradients = transpose(dene%dcoords)
    if (present(eigvals)) eigvals = res%mode_eigs
    if (present(eigvecs)) eigvecs = res%modes
end function cmbd_mbd_rsscs_energy

real(c_double) function cmbd_mbd_scs_energy(sys_cp, n_atoms, alpha_0, C6, damping_p, gradients) bind(c)
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
    res = mbd_scs_energy(sys, 'scs', alpha_0, C6, damping, dene, mbd_grad())
    cmbd_mbd_scs_energy = res%energy
end function cmbd_mbd_scs_energy

subroutine cmbd_dipole_matrix(sys_cp, damping_p, k_point, dipmat_p) bind(c)
    type(c_ptr), intent(in), value :: sys_cp
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(in), optional :: k_point(3)
    type(c_ptr), intent(in), value :: dipmat_p

    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damp
    type(mbd_matrix_real) :: dipmat
    type(mbd_matrix_complex) :: dipmat_c
    real(dp), pointer :: dipmat_re(:, :)
    complex(dp), pointer :: dipmat_cplx(:, :)
    integer :: n_atoms

    sys => get_mbd_system(sys_cp)
    n_atoms = size(sys%coords, 2)
    call c_f_pointer(damping_p, damp)
    if (present(k_point)) then
        dipmat_c = dipole_matrix(sys, damp, k_point=k_point)
        call c_f_pointer(dipmat_p, dipmat_cplx, [3*n_atoms, 3*n_atoms])
        dipmat_cplx = transpose(dipmat_c%val)
    else
        dipmat = dipole_matrix(sys, damp)
        call c_f_pointer(dipmat_p, dipmat_re, [3*n_atoms, 3*n_atoms])
        dipmat_re = transpose(dipmat%val)
    end if
end subroutine cmbd_dipole_matrix

real(c_double) function cmbd_coulomb_energy( &
        sys_cp, n_atoms, q, m, w_t, version, r_vdw, beta, a, C) bind(c)
    type(c_ptr), value :: sys_cp
    integer(c_int), value, intent(in) :: n_atoms
    real(c_double), value, intent(in) :: a, beta
    real(c_double), intent(in) ::  C(3*n_atoms, 3*n_atoms), &
        w_t(3*n_atoms), q(n_atoms), m(n_atoms), r_vdw(n_atoms)
    character(c_char), intent(in) :: version(20)

    type(mbd_system), pointer :: sys
    type(mbd_damping) :: damp

    damp%version = f_string(version)
    damp%r_vdw = r_vdw
    damp%ts_d = a
    damp%ts_sr = beta
    sys => get_mbd_system(sys_cp)
    cmbd_coulomb_energy = coulomb_energy(sys, q, m, w_t, C, damp)
end function cmbd_coulomb_energy

real(c_double) function cmbd_dipole_energy( &
        sys_cp, n_atoms, a0, w, w_t, version, r_vdw, beta, a, C) bind(c)
    type(c_ptr), value :: sys_cp
    integer(c_int), value, intent(in) :: n_atoms
    real(c_double), intent(in) :: C(3*n_atoms, 3*n_atoms), &
        w_t(3*n_atoms), w(n_atoms), a0(n_atoms), r_vdw(n_atoms)
    real(c_double), value, intent(in) :: a, beta
    character(c_char), intent(in) :: version(20)

    type(mbd_system), pointer :: sys
    type(mbd_damping) :: damp

    damp%version = f_string(version)
    damp%r_vdw = r_vdw
    damp%beta = beta
    damp%a = a
    sys => get_mbd_system(sys_cp)
    cmbd_dipole_energy = dipole_energy(sys, a0, w, w_t, C, damp)
end function cmbd_dipole_energy

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

subroutine f_c_string(str_f, str_c)
    character(len=*), intent(in) :: str_f
    character(kind=c_char), intent(out) :: str_c(:)

    integer :: i

    do i = 1, min(len(trim(str_f)), size(str_c)-1)
        str_c(i) = str_f(i:i)
    end do
    str_c(i) = c_null_char
end subroutine

end module mbd_c_api
