! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "version.h"

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
use mbd_ts, only: get_ts_energy
use mbd_utils, only: result_t

implicit none

private
public :: cmbd_with_scalapack, cmbd_with_mpi, cmbd_version_major, &
    cmbd_version_minor, cmbd_version_patch
public :: cmbd_init_geom, cmbd_destroy_geom, cmbd_init_damping, &
    cmbd_destroy_damping, cmbd_get_exception, cmbd_update_coords, cmbd_update_lattice, &
    cmbd_get_results, cmbd_destroy_result
public :: cmbd_ts_energy, cmbd_mbd_energy, cmbd_mbd_scs_energy, &
    cmbd_dipole_matrix, cmbd_coulomb_energy, cmbd_dipole_energy

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

integer(c_int), bind(c) :: cmbd_version_major = MBD_VERSION_MAJOR
integer(c_int), bind(c) :: cmbd_version_minor = MBD_VERSION_MINOR
integer(c_int), bind(c) :: cmbd_version_patch = MBD_VERSION_PATCH

contains

type(c_ptr) function cmbd_init_geom( &
    n_atoms, coords, lattice, k_grid, n_kpts, custom_k_pts, &
    n_freq, do_rpa, get_spectrum, get_rpa_orders, rpa_rescale_eigs, &
    max_atoms_per_block &
) bind(c)
    integer(c_int), value, intent(in) :: n_atoms
    real(c_double), intent(in) :: coords(3, n_atoms)
    real(c_double), optional, intent(in) :: lattice(3, 3)
    integer(c_int), optional, intent(in) :: k_grid(3)
    integer(c_int), value, intent(in) :: n_kpts
    real(c_double), optional, intent(in) :: custom_k_pts(3, n_kpts)
    integer(c_int), value, intent(in) :: n_freq
    logical(c_bool), value, intent(in) :: do_rpa
    logical(c_bool), value, intent(in) :: get_spectrum
    logical(c_bool), value, intent(in) :: get_rpa_orders
    logical(c_bool), value, intent(in) :: rpa_rescale_eigs
    integer(c_int), value, intent(in) :: max_atoms_per_block

    type(geom_t), pointer :: geom

    allocate (geom)
    geom%coords = coords
    if (present(lattice)) geom%lattice = lattice
    if (present(k_grid)) geom%k_grid = k_grid
    if (present(custom_k_pts)) geom%custom_k_pts = custom_k_pts
    if (n_freq > 0) geom%param%n_freq = n_freq
#ifdef WITH_SCALAPACK
    if (max_atoms_per_block > 0) geom%max_atoms_per_block = max_atoms_per_block
#endif
    geom%do_rpa = do_rpa
    geom%get_eigs = get_spectrum
    geom%get_modes = get_spectrum
    geom%get_rpa_orders = get_rpa_orders
    geom%param%rpa_rescale_eigs = rpa_rescale_eigs
    call geom%init()
    cmbd_init_geom = c_loc(geom)
end function

subroutine cmbd_update_coords(geom_c, coords_c) bind(c)
    type(c_ptr), value, intent(in) :: geom_c
    type(c_ptr), value, intent(in) :: coords_c

    type(geom_t), pointer :: geom
    real(c_double), pointer :: coords(:, :)

    call c_f_pointer(geom_c, geom)
    call c_f_pointer(coords_c, coords, [3, geom%siz()])
    geom%coords = coords
end subroutine

subroutine cmbd_update_lattice(geom_c, lattice) bind(c)
    type(c_ptr), value, intent(in) :: geom_c
    real(c_double), intent(in) :: lattice(3, 3)

    type(geom_t), pointer :: geom

    call c_f_pointer(geom_c, geom)
    geom%lattice = lattice
end subroutine

subroutine cmbd_destroy_geom(geom_c) bind(c)
    type(c_ptr), value, intent(in) :: geom_c

    type(geom_t), pointer :: geom

    call c_f_pointer(geom_c, geom)
    call geom%destroy()
    deallocate (geom)
end subroutine

subroutine cmbd_get_exception(geom_c, code, origin, msg) bind(c)
    type(c_ptr), value, intent(in) :: geom_c
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
    integer(c_int), value, intent(in) :: n_atoms
    character(kind=c_char), intent(in) :: version_c(*)
    real(c_double), optional, intent(in) :: r_vdw(n_atoms)
    real(c_double), optional, intent(in) :: sigma(n_atoms)
    real(c_double), value, intent(in) :: beta
    real(c_double), value, intent(in) :: a

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
    type(c_ptr), value, intent(in) :: damping_c

    type(damping_t), pointer :: damping

    call c_f_pointer(damping_c, damping)
    deallocate (damping)
end subroutine

real(c_double) function cmbd_ts_energy(geom_c, alpha_0_c, C6_c, damping_c) bind(c)
    type(c_ptr), value, intent(in) :: geom_c
    type(c_ptr), value, intent(in) :: alpha_0_c
    type(c_ptr), value, intent(in) :: C6_c
    type(c_ptr), value, intent(in) :: damping_c

    type(geom_t), pointer :: geom
    real(c_double), pointer :: alpha_0(:)
    real(c_double), pointer :: C6(:)
    type(damping_t), pointer :: damping

    call c_f_pointer(geom_c, geom)
    call c_f_pointer(alpha_0_c, alpha_0, [geom%siz()])
    call c_f_pointer(C6_c, C6, [geom%siz()])
    call c_f_pointer(damping_c, damping)
    cmbd_ts_energy = get_ts_energy(geom, alpha_0, C6, damping)
end function

type(c_ptr) function cmbd_mbd_energy(geom_c, alpha_0_c, C6_c, damping_c, grad) bind(c)
    type(c_ptr), value, intent(in) :: geom_c
    type(c_ptr), value, intent(in) :: alpha_0_c
    type(c_ptr), value, intent(in) :: C6_c
    type(c_ptr), value, intent(in) :: damping_c
    logical(c_bool), value, intent(in) :: grad

    type(geom_t), pointer :: geom
    real(c_double), pointer :: alpha_0(:)
    real(c_double), pointer :: C6(:)
    type(damping_t), pointer :: damping
    type(result_t), pointer :: res

    call c_f_pointer(geom_c, geom)
    call c_f_pointer(alpha_0_c, alpha_0, [geom%siz()])
    call c_f_pointer(C6_c, C6, [geom%siz()])
    call c_f_pointer(damping_c, damping)
    allocate (res)
    res = get_mbd_energy( &
        geom, alpha_0, C6, damping, grad_request_t( &
            dcoords=grad, dlattice=grad .and. allocated(geom%lattice) &
        ) &
    )
    cmbd_mbd_energy = c_loc(res)
end function

type(c_ptr) function cmbd_mbd_scs_energy( &
        geom_c, variant_c, alpha_0_c, C6_c, damping_c, grad) bind(c)
    type(c_ptr), value, intent(in) :: geom_c
    character(kind=c_char), intent(in) :: variant_c(*)
    type(c_ptr), value, intent(in) :: alpha_0_c
    type(c_ptr), value, intent(in) :: C6_c
    type(c_ptr), value, intent(in) :: damping_c
    logical(c_bool), value, intent(in) :: grad

    type(geom_t), pointer :: geom
    character(len=20) :: variant
    real(c_double), pointer :: alpha_0(:)
    real(c_double), pointer :: C6(:)
    type(damping_t), pointer :: damping
    type(result_t), pointer :: res

    call c_f_pointer(geom_c, geom)
    variant = f_string(variant_c)
    call c_f_pointer(alpha_0_c, alpha_0, [geom%siz()])
    call c_f_pointer(C6_c, C6, [geom%siz()])
    call c_f_pointer(damping_c, damping)
    allocate (res)
    res = get_mbd_scs_energy( &
        geom, variant, alpha_0, C6, damping, grad_request_t( &
            dcoords=grad, dlattice=grad .and. allocated(geom%lattice) &
        ) &
    )
    cmbd_mbd_scs_energy = c_loc(res)
end function

subroutine cmbd_get_results( &
    res_c, energy, gradients_c, lattice_gradients_c, eigvals_c, eigvecs_c, rpa_orders_c, &
    eigvals_k_c, eigvecs_k_c &
) bind(c)
    type(c_ptr), value, intent(in) :: res_c
    real(c_double), intent(out) :: energy
    type(c_ptr), value, intent(in) :: gradients_c
    type(c_ptr), value, intent(in) :: lattice_gradients_c
    type(c_ptr), value, intent(in) :: eigvals_c
    type(c_ptr), value, intent(in) :: eigvecs_c
    type(c_ptr), value, intent(in) :: rpa_orders_c
    type(c_ptr), value, intent(in) :: eigvals_k_c
    type(c_ptr), value, intent(in) :: eigvecs_k_c

    type(result_t), pointer :: res
    real(c_double), pointer :: gradients(:, :)
    real(c_double), pointer :: lattice_gradients(:, :)
    real(c_double), pointer :: eigvals(:)
    real(c_double), pointer :: eigvecs(:, :)
    real(c_double), pointer :: rpa_orders(:)
    real(c_double), pointer :: eigvals_k(:, :)
    complex(c_double_complex), pointer :: eigvecs_k(:, :, :)

    call c_f_pointer(res_c, res)
    energy = res%energy
    if (c_associated(gradients_c) .and. allocated(res%dE%dcoords)) then
        call c_f_pointer(gradients_c, gradients, [3, size(res%dE%dcoords, 1)])
        gradients = transpose(res%dE%dcoords)
    end if
    if (c_associated(lattice_gradients_c) .and. allocated(res%dE%dlattice)) then
        call c_f_pointer(lattice_gradients_c, lattice_gradients, [3, 3])
        lattice_gradients = transpose(res%dE%dlattice)
    end if
    if (c_associated(eigvals_c) .and. allocated(res%mode_eigs)) then
        call c_f_pointer(eigvals_c, eigvals, [size(res%mode_eigs)])
        eigvals = res%mode_eigs
    end if
    if (c_associated(eigvecs_c) .and. allocated(res%modes)) then
        call c_f_pointer(eigvecs_c, eigvecs, [size(res%modes, 1), size(res%modes, 2)])
        eigvecs = res%modes
    end if
    if (c_associated(rpa_orders_c) .and. allocated(res%rpa_orders)) then
        call c_f_pointer(rpa_orders_c, rpa_orders, [size(res%rpa_orders)])
        rpa_orders = res%rpa_orders
    end if
    if (c_associated(eigvals_k_c) .and. allocated(res%mode_eigs_k)) then
        call c_f_pointer( &
            eigvals_k_c, eigvals_k, &
            [size(res%mode_eigs_k, 1), size(res%mode_eigs_k, 2)] &
        )
        eigvals_k = res%mode_eigs_k
    end if
    if (c_associated(eigvecs_k_c) .and. allocated(res%modes_k)) then
        call c_f_pointer( &
            eigvecs_k_c, eigvecs_k, &
            [size(res%modes_k, 1), size(res%modes_k, 2), size(res%modes_k, 3)] &
        )
        eigvecs_k = res%modes_k
    end if
end subroutine

subroutine cmbd_destroy_result(res_c) bind(c)
    type(c_ptr), value, intent(in) :: res_c

    type(result_t), pointer :: res

    call c_f_pointer(res_c, res)
    deallocate (res)
end subroutine

subroutine cmbd_dipole_matrix(geom_c, damping_c, q_point, dipmat_c) bind(c)
    type(c_ptr), value, intent(in) :: geom_c
    type(c_ptr), value, intent(in) :: damping_c
    real(c_double), optional, intent(in) :: q_point(3)
    type(c_ptr), value, intent(in) :: dipmat_c

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
