! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_c_api

use iso_c_binding, only: c_ptr, c_int, c_double, c_f_pointer, c_loc, c_bool
use mbd, only: mbd_system, mbd_calc, mbd_damping, get_mbd_energy, init_grid, &
    mbd_rsscs_energy, mbd_scs_energy

implicit none

contains

type(c_ptr) function mbd_init_calc() bind(c)
    type(mbd_calc), pointer :: calc

    allocate (calc)
    call init_grid(calc)
    mbd_init_calc = c_loc(calc)
end function mbd_init_calc

subroutine mbd_destroy_calc(calc_p) bind(c)
    type(c_ptr), value :: calc_p

    type(mbd_calc), pointer :: calc

    call c_f_pointer(calc_p, calc)
    deallocate (calc)
end subroutine mbd_destroy_calc

type(c_ptr) function mbd_init_system(calc_p, n_atoms, coords, periodic, lattice, k_grid) bind(c)
    type(c_ptr), value :: calc_p
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: coords(n_atoms, 3)
    logical(c_bool), intent(in), value :: periodic
    real(c_double), intent(in) :: lattice(3, 3)
    integer(c_int), intent(in) :: k_grid(3)

    type(mbd_calc), pointer :: calc
    type(mbd_system), pointer :: sys

    call c_f_pointer(calc_p, calc)
    allocate (sys)
    sys%calc => calc
    sys%coords = coords
    if (periodic) then
        sys%periodic = .true.
        sys%lattice = lattice
        sys%k_grid = k_grid
    end if
    mbd_init_system = c_loc(sys)
end function mbd_init_system

subroutine mbd_destroy_system(sys_p) bind(c)
    type(c_ptr), value :: sys_p

    type(mbd_system), pointer :: sys

    call c_f_pointer(sys_p, sys)
    deallocate (sys)
end subroutine mbd_destroy_system

type(c_ptr) function mbd_init_damping(n_atoms, r_vdw, beta, a) bind(c)
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: r_vdw(n_atoms)
    real(c_double), intent(in), value :: beta
    real(c_double), intent(in), value :: a

    type(mbd_damping), pointer :: damping

    allocate (damping)
    damping%version = 'fermi,dip'
    damping%r_vdw = r_vdw
    damping%beta = beta
    damping%a = a
    mbd_init_damping = c_loc(damping)
end function mbd_init_damping

subroutine mbd_destroy_damping(damping_p) bind(c)
    type(c_ptr), value :: damping_p

    type(mbd_damping), pointer :: damping

    call c_f_pointer(damping_p, damping)
    deallocate (damping%r_vdw)
    deallocate (damping)
end subroutine mbd_destroy_damping

real(c_double) function calc_mbd_energy(sys_p, n_atoms, alpha_0, omega, damping_p) bind(c)
    type(c_ptr), intent(in), value :: sys_p
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: omega(n_atoms)
    type(c_ptr), intent(in), value :: damping_p

    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping

    call c_f_pointer(sys_p, sys)
    call c_f_pointer(damping_p, damping)
    calc_mbd_energy = get_mbd_energy(sys, alpha_0, omega, damping)
end function calc_mbd_energy

real(c_double) function calc_rpa_energy(sys_p, n_atoms, alpha_0, omega, damping_p) bind(c)
    type(c_ptr), intent(in), value :: sys_p
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: omega(n_atoms)
    type(c_ptr), intent(in), value :: damping_p

    type(mbd_system), pointer :: sys
    type(mbd_system) :: sys2
    type(mbd_damping), pointer :: damping

    call c_f_pointer(sys_p, sys)
    call c_f_pointer(damping_p, damping)
    sys2 = sys
    sys2%do_rpa = .true.
    calc_rpa_energy = get_mbd_energy(sys2, alpha_0, omega, damping)
end function calc_rpa_energy

real(c_double) function calc_mbd_rsscs_energy(sys_p, n_atoms, alpha_0, omega, damping_p) bind(c)
    type(c_ptr), intent(in), value :: sys_p
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: omega(n_atoms)
    type(c_ptr), intent(in), value :: damping_p

    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping

    call c_f_pointer(sys_p, sys)
    call c_f_pointer(damping_p, damping)
    calc_mbd_rsscs_energy = mbd_rsscs_energy(sys, alpha_0, omega, damping)
end function calc_mbd_rsscs_energy

real(c_double) function calc_mbd_scs_energy(sys_p, n_atoms, alpha_0, omega, damping_p) bind(c)
    type(c_ptr), intent(in), value :: sys_p
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: omega(n_atoms)
    type(c_ptr), intent(in), value :: damping_p

    type(mbd_system), pointer :: sys
    type(mbd_damping), pointer :: damping

    call c_f_pointer(sys_p, sys)
    call c_f_pointer(damping_p, damping)
    calc_mbd_scs_energy = mbd_scs_energy(sys, alpha_0, omega, damping)
end function calc_mbd_scs_energy

end module mbd_c_api
