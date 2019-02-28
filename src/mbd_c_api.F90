! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_c_api
!! Implementation of C API.

use iso_c_binding
use mbd_constants
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
public :: cmbd_init_geom, cmbd_destroy_geom, cmbd_init_damping, &
    cmbd_destroy_damping, cmbd_get_exception
public :: cmbd_ts_energy, cmbd_mbd_energy, cmbd_rpa_energy, cmbd_mbd_rsscs_energy, &
    cmbd_mbd_scs_energy, cmbd_dipole_matrix, cmbd_coulomb_energy, cmbd_dipole_energy

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

contains

type(c_ptr) function cmbd_init_geom(n_atoms, coords, lattice, k_grid, n_freq) bind(c)
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: coords(3, n_atoms)
    real(c_double), intent(in), optional :: lattice(3, 3)
    integer(c_int), intent(in), optional :: k_grid(3)
    integer(c_int), intent(in), value :: n_freq

    type(geom_t), pointer :: geom

    allocate (geom)
    if (n_freq > 0) geom%param%n_freq = n_freq
    geom%coords = coords
    if (present(lattice)) geom%lattice = lattice
    if (present(k_grid)) geom%k_grid = k_grid
    call geom%init()
    cmbd_init_geom = c_loc(geom)
end function

subroutine cmbd_destroy_geom(geom_c) bind(c)
    type(c_ptr), value :: geom_c

    type(geom_t), pointer :: geom

    call c_f_pointer(geom_c, geom)
    call geom%destroy()
    deallocate (geom)
end subroutine

subroutine cmbd_get_exception(geom_c, code, origin, msg) bind(c)
    type(c_ptr), value :: geom_c
    integer(c_int), intent(out) :: code
    character(kind=c_char), intent(out) :: origin(50), msg(150)

    type(geom_t), pointer :: geom

    call c_f_pointer(geom_c, geom)
    code = geom%exc%code
    call f_c_string(geom%exc%origin, origin)
    call f_c_string(geom%exc%msg, msg)
    geom%exc%code = 0
    geom%exc%origin = ''
    geom%exc%msg = ''
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

subroutine cmbd_destroy_damping(damping_c) bind(c)
    type(c_ptr), value :: damping_c

    type(damping_t), pointer :: damping

    call c_f_pointer(damping_c, damping)
    if (allocated(damping%r_vdw)) deallocate (damping%r_vdw)
    if (allocated(damping%sigma)) deallocate (damping%sigma)
    deallocate (damping)
end subroutine

real(c_double) function cmbd_ts_energy(geom_c, n_atoms, alpha_0, C6, damping_c) bind(c)
    type(c_ptr), intent(in), value :: geom_c
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_c

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping

    call c_f_pointer(geom_c, geom)
    call c_f_pointer(damping_c, damping)
    cmbd_ts_energy = ts_energy(geom, alpha_0, C6, damping)
end function

real(c_double) function cmbd_mbd_energy( &
    geom_c, n_atoms, alpha_0, C6, damping_c, gradients, latt_gradients &
) bind(c)
    type(c_ptr), intent(in), value :: geom_c
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_c
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: latt_gradients(3, 3)

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping
    type(result_t) :: res
    type(grad_t) :: dene

    call c_f_pointer(geom_c, geom)
    call c_f_pointer(damping_c, damping)
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
    geom_c, n_atoms, alpha_0, C6, damping_c, gradients, latt_gradients, &
    rpa_orders &
) bind(c)
    type(c_ptr), intent(in), value :: geom_c
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_c
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: latt_gradients(3, 3)
    real(c_double), intent(out), optional :: rpa_orders(10)

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping
    type(result_t) :: res
    type(grad_t) :: dene

    call c_f_pointer(geom_c, geom)
    call c_f_pointer(damping_c, damping)
    geom%do_rpa = .true.
    geom%get_rpa_orders = present(rpa_orders)
    res = get_mbd_energy(geom, alpha_0, C6, damping, dene, grad_request_t())
    geom%do_rpa = .false.
    geom%get_rpa_orders = .false.
    cmbd_rpa_energy = res%energy
    if (present(rpa_orders)) rpa_orders = res%rpa_orders
end function

real(c_double) function cmbd_mbd_rsscs_energy( &
    geom_c, n_atoms, alpha_0, C6, damping_c, &
    gradients, latt_gradients, eigvals, eigvecs &
) bind(c)
    type(c_ptr), intent(in), value :: geom_c
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_c
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: latt_gradients(3, 3)
    real(c_double), intent(out), optional :: eigvals(3*n_atoms)
    real(c_double), intent(out), optional :: eigvecs(3*n_atoms, 3*n_atoms)

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping
    type(result_t) :: res
    type(grad_t) :: dene

    call c_f_pointer(geom_c, geom)
    call c_f_pointer(damping_c, damping)
    geom%get_eigs = present(eigvals)
    geom%get_modes = present(eigvecs)
    res = get_mbd_scs_energy( &
        geom, 'rsscs', alpha_0, C6, damping, dene, grad_request_t( &
            dcoords=present(gradients), &
            dlattice=present(latt_gradients) &
        ) &
    )
    geom%get_eigs = .false.
    geom%get_modes = .false.
    if (geom%has_exc()) return
    cmbd_mbd_rsscs_energy = res%energy
    if (present(gradients)) gradients = transpose(dene%dcoords)
    if (present(latt_gradients)) latt_gradients = transpose(dene%dlattice)
    if (present(eigvals)) eigvals = res%mode_eigs
    if (present(eigvecs)) eigvecs = res%modes
end function

real(c_double) function cmbd_mbd_scs_energy( &
    geom_c, n_atoms, alpha_0, C6, damping_c, gradients, latt_gradients &
) bind(c)
    type(c_ptr), intent(in), value :: geom_c
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: C6(n_atoms)
    type(c_ptr), intent(in), value :: damping_c
    real(c_double), intent(out), optional :: gradients(3, n_atoms)
    real(c_double), intent(out), optional :: latt_gradients(3, 3)

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damping
    type(result_t) :: res
    type(grad_t) :: dene

    call c_f_pointer(geom_c, geom)
    call c_f_pointer(damping_c, damping)
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

subroutine cmbd_dipole_matrix(geom_c, damping_c, q_point, dipmat_c) bind(c)
    type(c_ptr), intent(in), value :: geom_c
    type(c_ptr), intent(in), value :: damping_c
    real(c_double), intent(in), optional :: q_point(3)
    type(c_ptr), intent(in), value :: dipmat_c

    type(geom_t), pointer :: geom
    type(damping_t), pointer :: damp
    type(matrix_re_t) :: dipmat_re
    type(matrix_cplx_t) :: dipmat_cplx
    real(dp), pointer :: dipmat_re_p(:, :)
    complex(dp), pointer :: dipmat_cplx_p(:, :)
    integer :: n_atoms

    call c_f_pointer(geom_c, geom)
    n_atoms = size(geom%coords, 2)
    call c_f_pointer(damping_c, damp)
    if (present(q_point)) then
        dipmat_cplx = dipole_matrix(geom, damp, q=q_point)
        call c_f_pointer(dipmat_c, dipmat_cplx_p, [3*n_atoms, 3*n_atoms])
        dipmat_cplx_p = transpose(dipmat_cplx%val)
    else
        dipmat_re = dipole_matrix(geom, damp)
        call c_f_pointer(dipmat_c, dipmat_re_p, [3*n_atoms, 3*n_atoms])
        dipmat_re_p = transpose(dipmat_re%val)
    end if
end subroutine

real(c_double) function cmbd_coulomb_energy( &
        geom_c, n_atoms, q, m, w_t, version, r_vdw, beta, a, C) bind(c)
    type(c_ptr), value :: geom_c
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
    call c_f_pointer(geom_c, geom)
    cmbd_coulomb_energy = coulomb_energy(geom, q, m, w_t, C, damp)
end function

real(c_double) function cmbd_dipole_energy( &
        geom_c, n_atoms, a0, w, w_t, version, r_vdw, beta, a, C) bind(c)
    type(c_ptr), value :: geom_c
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
    call c_f_pointer(geom_c, geom)
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
