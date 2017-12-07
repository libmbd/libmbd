module mbd_c_api

use iso_c_binding, only: c_ptr, c_int, c_double, c_f_pointer, c_loc
use mbd, only: mbd_system, mbd_calc, mbd_damping, get_mbd_energy, init_grid, &
    destroy_grid, mbd_rsscs_energy

implicit none

contains

type(c_ptr) function mbd_init_calc(n_grid) bind(c)
    integer(c_int), intent(in), value :: n_grid

    type(mbd_calc), pointer :: calc

    allocate (calc)
    call init_grid(calc, n_grid)
    mbd_init_calc = c_loc(calc)
end function mbd_init_calc

subroutine mbd_destroy_calc(calc_p) bind(c)
    type(c_ptr), value :: calc_p

    type(mbd_calc), pointer :: calc

    call c_f_pointer(calc_p, calc)
    call destroy_grid(calc)
    deallocate (calc)
end subroutine mbd_destroy_calc

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

subroutine calc_mbd_energy(calc_p, n_atoms, coords, alpha_0, omega, damping_p, energy) bind(c)
    type(c_ptr), intent(in), value :: calc_p
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: coords(n_atoms, 3)
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: omega(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out) :: energy

    type(mbd_system) :: sys
    type(mbd_damping), pointer :: damping

    call c_f_pointer(calc_p, sys%calc)
    call c_f_pointer(damping_p, damping)
    sys%coords = coords
    energy = get_mbd_energy(sys, alpha_0, omega, damping)
end subroutine calc_mbd_energy

subroutine calc_mbd_rsscs_energy(calc_p, n_atoms, coords, alpha_0, omega, damping_p, energy) bind(c)
    type(c_ptr), intent(in), value :: calc_p
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: coords(n_atoms, 3)
    real(c_double), intent(in) :: alpha_0(n_atoms)
    real(c_double), intent(in) :: omega(n_atoms)
    type(c_ptr), intent(in), value :: damping_p
    real(c_double), intent(out) :: energy

    type(mbd_system) :: sys
    type(mbd_damping), pointer :: damping

    call c_f_pointer(calc_p, sys%calc)
    call c_f_pointer(damping_p, damping)
    sys%coords = coords
    energy = mbd_rsscs_energy(sys, alpha_0, omega, damping)
end subroutine calc_mbd_rsscs_energy

end module mbd_c_api
