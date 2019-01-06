! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_c_api
!! Implementation of C API.

use iso_c_binding
use mbd_constants
use mbd_calc, only: calc_t
use mbd_coulomb, only: dipole_energy, coulomb_energy
use mbd_damping, only: damping_t
use mbd_dipole, only: dipole_matrix
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_t, grad_request_t
use mbd_matrix, only: matrix_re_t, matrix_cplx_t
use mbd_methods, only: get_mbd_energy, get_mbd_scs_energy
use mbd_ts, only: ts_energy
use mbd_utils, only: result_t

implicit none

private
public :: cmbd_with_scalapack, cmbd_with_mpi
public :: cmbd_calc, cmbd_geom
public :: cmbd_init_calc, cmbd_destroy_calc, cmbd_init_geom, &
    cmbd_destroy_geom, cmbd_init_damping, cmbd_destroy_damping, cmbd_get_exception
public :: cmbd_ts_energy, cmbd_mbd_energy, cmbd_rpa_energy, cmbd_mbd_rsscs_energy, &
    cmbd_mbd_scs_energy, cmbd_dipole_matrix, cmbd_coulomb_energy, cmbd_dipole_energy, &
    cmbd_toggle_muted

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

type, bind(c) :: cmbd_geom
    type(c_ptr) :: mbd_geom_f = c_null_ptr
end type

contains

type(c_ptr) function cmbd_init_calc(n_freq) bind(c)
    integer(c_int), intent(in), value :: n_freq

    type(calc_t), pointer :: calc
    type(cmbd_calc), pointer :: calc_c

    allocate (calc)
    calc%param%n_frequency_grid = n_freq
    call calc%init()
    allocate (calc_c)
    calc_c%n_freq = ubound(calc%omega_grid, 1)
    calc_c%omega_grid = c_loc(calc%omega_grid)
    calc_c%omega_grid_w = c_loc(calc%omega_grid_w)
    calc_c%mbd_calc_f = c_loc(calc)
    cmbd_init_calc = c_loc(calc_c)
end function

subroutine cmbd_destroy_calc(calc_cp) bind(c)
    type(c_ptr), value :: calc_cp

    type(calc_t), pointer :: calc
    type(cmbd_calc), pointer :: calc_c

    call c_f_pointer(calc_cp, calc_c)
    call c_f_pointer(calc_c%mbd_calc_f, calc)
    deallocate (calc)
    deallocate (calc_c)
end subroutine

subroutine cmbd_get_exception(calc_cp, code, origin, msg) bind(c)
    type(c_ptr), value :: calc_cp
    integer(c_int), intent(out) :: code
    character(kind=c_char), intent(out) :: origin(50), msg(150)

    type(calc_t), pointer :: calc

    calc => get_mbd_calc(calc_cp)
    code = calc%exc%code
    call f_c_string(calc%exc%origin, origin)
    call f_c_string(calc%exc%msg, msg)
    calc%exc%code = 0
    calc%exc%origin = ''
    calc%exc%msg = ''
end subroutine

subroutine cmbd_toggle_muted(calc_cp) bind(c)
    type(c_ptr), value :: calc_cp

    type(calc_t), pointer :: calc

    calc => get_mbd_calc(calc_cp)
    calc%muted = .not. calc%muted
end subroutine

type(c_ptr) function cmbd_init_geom( &
        calc_cp, n_atoms, coords, lattice, k_grid) bind(c)
    type(c_ptr), value :: calc_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: coords(3, n_atoms)
    real(c_double), intent(in), optional :: lattice(3, 3)
    integer(c_int), intent(in), optional :: k_grid(3)

    type(calc_t), pointer :: calc
    type(geom_t), pointer :: geom
    type(cmbd_geom), pointer :: geom_c

    calc => get_mbd_calc(calc_cp)
    allocate (geom)
    geom%coords = coords
    if (present(lattice)) geom%lattice = lattice
    if (present(k_grid)) geom%k_grid = k_grid
    call geom%init(calc)
    allocate (geom_c)
    geom_c%mbd_geom_f = c_loc(geom)
    cmbd_init_geom = c_loc(geom_c)
end function

subroutine cmbd_destroy_geom(geom_cp) bind(c)
    type(c_ptr), value :: geom_cp

    type(cmbd_geom), pointer :: geom_c
    type(geom_t), pointer :: geom

    call c_f_pointer(geom_cp, geom_c)
    call c_f_pointer(geom_c%mbd_geom_f, geom)
    call geom%destroy()
    deallocate (geom)
    deallocate (geom_c)
end subroutine

type(c_ptr) function cmbd_init_damping(n_atoms, version_c, r_vdw, sigma, beta, a) bind(c)
    integer(c_int), intent(in), value :: n_atoms
    character(kind=c_char), intent(in) :: version_c(*)
    real(c_double), intent(in), optional :: r_vdw(n_atoms)
    real(c_double), intent(in), optional :: sigma(n_atoms)
    real(c_double), intent(in), value :: beta
    real(c_double), intent(in), value :: a

    type(damping_t), pointer :: damping

    allocate (damping)
    damping%version = f_string(version_c)
    if (present(r_vdw)) damping%r_vdw = r_vdw
    if (present(sigma)) damping%sigma = sigma
    damping%beta = beta
    damping%a = a
    damping%ts_sr = beta
    damping%ts_d = a
    cmbd_init_damping = c_loc(damping)
end function

subroutine cmbd_destroy_damping(damping_p) bind(c)
    type(c_ptr), value :: damping_p

    type(damping_t), pointer :: damping

    call c_f_pointer(damping_p, damping)
    if (allocated(damping%r_vdw)) deallocate (damping%r_vdw)
    if (allocated(damping%sigma)) deallocate (damping%sigma)
    deallocate (damping)
end subroutine

function get_mbd_geom(geom_cp)
    type(c_ptr), intent(in), value :: geom_cp
    type(geom_t), pointer :: get_mbd_geom

    type(cmbd_geom), pointer :: geom_c

    call c_f_pointer(geom_cp, geom_c)
    call c_f_pointer(geom_c%mbd_geom_f, get_mbd_geom)
end function

function get_mbd_calc(calc_cp)
    type(c_ptr), intent(in), value :: calc_cp
    type(calc_t), pointer :: get_mbd_calc

    type(cmbd_calc), pointer :: calc_c

    call c_f_pointer(calc_cp, calc_c)
    call c_f_pointer(calc_c%mbd_calc_f, get_mbd_calc)
end function

real(c_double) function cmbd_ts_energy(geom_cp, n_atoms, alpha_0, C6, damping_p) bind(c)
    type(c_ptr), intent(in), value :: geom_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p

    type(cmbd_geom), pointer :: geom_c
    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping

    call c_f_pointer(geom_cp, geom_c)
    call c_f_pointer(geom_c%mbd_geom_f, geom)
    call c_f_pointer(damping_p, damping)
    cmbd_ts_energy = ts_energy(geom, alpha_0, C6, damping)
end function

real(c_double) function cmbd_mbd_energy( &
    geom_cp, n_atoms, alpha_0, C6, damping_p, gradients, latt_gradients &
) bind(c)
    type(c_ptr), intent(in), value :: geom_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: latt_gradients(3, 3)

    type(cmbd_geom), pointer :: geom_c
    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping
    type(result_t) :: res
    type(grad_t) :: dene

    call c_f_pointer(geom_cp, geom_c)
    call c_f_pointer(geom_c%mbd_geom_f, geom)
    call c_f_pointer(damping_p, damping)
    res = get_mbd_energy( &
        geom, alpha_0, C6, damping, dene, grad_request_t( &
            dcoords=present(gradients), &
            dlattice=present(latt_gradients) &
        ) &
    )
    cmbd_mbd_energy = res%energy
    if (present(gradients)) gradients = transpose(dene%dcoords)
    if (present(latt_gradients)) latt_gradients = transpose(dene%dlattice)
end function

real(c_double) function cmbd_rpa_energy( &
    geom_cp, n_atoms, alpha_0, C6, damping_p, gradients, latt_gradients &
) bind(c)
    type(c_ptr), intent(in), value :: geom_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: latt_gradients(3, 3)

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping
    type(result_t) :: res
    type(grad_t) :: dene

    geom => get_mbd_geom(geom_cp)
    call c_f_pointer(damping_p, damping)
    geom%calc%do_rpa = .true.
    res = get_mbd_energy(geom, alpha_0, C6, damping, dene, grad_request_t())
    geom%calc%do_rpa = .false.
    cmbd_rpa_energy = res%energy
end function

real(c_double) function cmbd_mbd_rsscs_energy( &
    geom_cp, n_atoms, alpha_0, C6, damping_p, &
    gradients, latt_gradients, eigvals, eigvecs &
) bind(c)
    type(c_ptr), intent(in), value :: geom_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: latt_gradients(3, 3)
    real(c_double), intent(out), optional :: eigvals(3*n_atoms)
    real(c_double), intent(out), optional :: eigvecs(3*n_atoms, 3*n_atoms)

    type(cmbd_geom), pointer :: geom_c
    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping
    type(result_t) :: res
    type(grad_t) :: dene

    call c_f_pointer(geom_cp, geom_c)
    call c_f_pointer(geom_c%mbd_geom_f, geom)
    call c_f_pointer(damping_p, damping)
    geom%calc%get_eigs = present(eigvals)
    geom%calc%get_modes = present(eigvecs)
    res = get_mbd_scs_energy( &
        geom, 'rsscs', alpha_0, C6, damping, dene, grad_request_t( &
            dcoords=present(gradients), &
            dlattice=present(latt_gradients) &
        ) &
    )
    if (geom%has_exc()) return
    cmbd_mbd_rsscs_energy = res%energy
    if (present(gradients)) gradients = transpose(dene%dcoords)
    if (present(latt_gradients)) latt_gradients = transpose(dene%dlattice)
    if (present(eigvals)) eigvals = res%mode_eigs
    if (present(eigvecs)) eigvecs = res%modes
end function

real(c_double) function cmbd_mbd_scs_energy( &
    geom_cp, n_atoms, alpha_0, C6, damping_p, gradients, latt_gradients &
) bind(c)
    type(c_ptr), intent(in), value :: geom_cp
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: latt_gradients(3, 3)

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping
    type(result_t) :: res
    type(grad_t) :: dene

    geom => get_mbd_geom(geom_cp)
    call c_f_pointer(damping_p, damping)
    res = get_mbd_scs_energy( &
        geom, 'scs', alpha_0, C6, damping, dene, grad_request_t( &
            dcoords=present(gradients), &
            dlattice=present(latt_gradients) &
        ) &
    )
    cmbd_mbd_scs_energy = res%energy
    if (present(gradients)) gradients = transpose(dene%dcoords)
    if (present(latt_gradients)) latt_gradients = transpose(dene%dlattice)
end function

subroutine cmbd_dipole_matrix(geom_cp, damping_p, q_point, dipmat_p) bind(c)
    type(c_ptr), intent(in), value :: geom_cp
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(in), optional :: q_point(3)
    type(c_ptr), intent(in), value :: dipmat_p

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damp
    type(matrix_re_t) :: dipmat
    type(matrix_cplx_t) :: dipmat_c
    real(dp), pointer :: dipmat_re(:, :)
    complex(dp), pointer :: dipmat_cplx(:, :)
    integer :: n_atoms

    geom => get_mbd_geom(geom_cp)
    n_atoms = size(geom%coords, 2)
    call c_f_pointer(damping_p, damp)
    if (present(q_point)) then
        dipmat_c = dipole_matrix(geom, damp, q=q_point)
        call c_f_pointer(dipmat_p, dipmat_cplx, [3*n_atoms, 3*n_atoms])
        dipmat_cplx = transpose(dipmat_c%val)
    else
        dipmat = dipole_matrix(geom, damp)
        call c_f_pointer(dipmat_p, dipmat_re, [3*n_atoms, 3*n_atoms])
        dipmat_re = transpose(dipmat%val)
    end if
end subroutine

real(c_double) function cmbd_coulomb_energy( &
        geom_cp, n_atoms, q, m, w_t, version, r_vdw, beta, a, C) bind(c)
    type(c_ptr), value :: geom_cp
    integer(c_int), value, intent(in) :: n_atoms
    real(c_double), value, intent(in) :: a, beta
    real(c_double), intent(in) ::  C(3*n_atoms, 3*n_atoms), &
        w_t(3*n_atoms), q(n_atoms), m(n_atoms), r_vdw(n_atoms)
    character(c_char), intent(in) :: version(20)

    type(geom_t), pointer :: geom
    type(damping_t) :: damp

    damp%version = f_string(version)
    damp%r_vdw = r_vdw
    damp%ts_d = a
    damp%ts_sr = beta
    geom => get_mbd_geom(geom_cp)
    cmbd_coulomb_energy = coulomb_energy(geom, q, m, w_t, C, damp)
end function

real(c_double) function cmbd_dipole_energy( &
        geom_cp, n_atoms, a0, w, w_t, version, r_vdw, beta, a, C) bind(c)
    type(c_ptr), value :: geom_cp
    integer(c_int), value, intent(in) :: n_atoms
    real(c_double), intent(in) :: C(3*n_atoms, 3*n_atoms), &
        w_t(3*n_atoms), w(n_atoms), a0(n_atoms), r_vdw(n_atoms)
    real(c_double), value, intent(in) :: a, beta
    character(c_char), intent(in) :: version(20)

    type(geom_t), pointer :: geom
    type(damping_t) :: damp

    damp%version = f_string(version)
    damp%r_vdw = r_vdw
    damp%beta = beta
    damp%a = a
    geom => get_mbd_geom(geom_cp)
    cmbd_dipole_energy = dipole_energy(geom, a0, w, w_t, C, damp)
end function

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
end function

subroutine f_c_string(str_f, str_c)
    character(len=*), intent(in) :: str_f
    character(kind=c_char), intent(out) :: str_c(:)

    integer :: i

    do i = 1, min(len(trim(str_f)), size(str_c)-1)
        str_c(i) = str_f(i:i)
    end do
    str_c(i) = c_null_char
end subroutine

end module
