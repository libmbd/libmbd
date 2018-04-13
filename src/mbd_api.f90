! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_api

use mbd, only: mbd_system, mbd_calc_inner => mbd_calc, mbd_damping, &
    mbd_rsscs_energy, get_ts_energy, get_damping_parameters, init_grid, &
    mbd_result, scale_TS
use mbd_common, only: dp
use mbd_types, only: vecn
use mbd_vdw_param, only: default_vdw_params, species_index
use mbd_defaults

implicit none

private
public :: mbd_input, mbd_calc  ! types
public :: &  ! subroutines
    mbd_init, mbd_update_coords, mbd_update_lattice_vectors, &
    mbd_update_vdw_params_custom, mbd_update_vdw_params_from_ratios, &
    mbd_get_energy, mbd_get_gradients, mbd_get_lattice_derivs, &
    mbd_get_spectrum_modes, mbd_get_damping_parameters, &
    mbd_get_free_vdw_params

type :: mbd_input
    integer :: comm  ! MPI communicator

    ! which calculation will be done (mbd|ts)
    character(len=30) :: dispersion_type = 'mbd'
    logical :: calculate_forces = .true.
    logical :: calculate_spectrum = .false.

    real(dp) :: ts_ene_acc = TS_ENERGY_ACCURACY  ! accuracy of TS energy
    real(dp) :: ts_f_acc = TS_FORCES_ACCURACY  ! accuracy of TS gradients
    integer :: n_omega_grid = N_FREQUENCY_GRID  ! number of frequency grid points
    ! off-gamma shift of k-points in units of inter-k-point distance
    real(dp) :: k_grid_shift = K_GRID_SHIFT

    ! TS damping parameters
    real(dp) :: ts_d = TS_DAMPING_D
    real(dp) :: ts_sr
    ! MBD damping parameters
    real(dp) :: mbd_a = MBD_DAMPING_A
    real(dp) :: mbd_beta

    ! lattice vectors as column vectors, unallocated when not periodic
    real(dp), allocatable :: lattice_vectors(:, :)
    integer  :: k_grid(3)  ! number of k-points along reciprocal axes
    ! is there vacuum along some axes in a periodic calculation
    logical  :: vacuum_axis(3) = [.false., .false., .false.]
end type

type mbd_calc
    private
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    type(vecn) :: alpha_0
    type(vecn) :: C6
    character(len=30) :: dispersion_type
    type(mbd_calc_inner) :: calc
    type(mbd_result) :: results
end type

contains


subroutine mbd_init(calc, input)
    type(mbd_calc), target, intent(out) :: calc
    type(mbd_input), intent(in) :: input

    calc%sys%calc => calc%calc
    calc%sys%calc%comm = input%comm
    calc%dispersion_type = input%dispersion_type
    calc%sys%do_force = input%calculate_forces
    if (input%calculate_spectrum) then
        calc%sys%get_eigs = .true.
        calc%sys%get_modes = .true.
    end if
    calc%sys%calc%param%ts_energy_accuracy = input%ts_ene_acc
    ! TODO ... = input%ts_f_acc
    calc%sys%calc%param%n_frequency_grid = input%n_omega_grid
    calc%sys%calc%param%k_grid_shift = input%k_grid_shift
    calc%damp%beta = input%mbd_beta
    calc%damp%a = input%mbd_a
    calc%damp%ts_d = input%ts_d
    calc%damp%ts_sr = input%ts_sr
    calc%sys%k_grid = input%k_grid
    calc%sys%vacuum_axis = input%vacuum_axis
    call init_grid(calc%calc)
end subroutine


subroutine mbd_update_coords(calc, coords)
    type(mbd_calc), intent(inout) :: calc
    real(dp), intent(in) :: coords(:, :)

    allocate (calc%sys%coords(3, size(coords, 2)))
    calc%sys%coords = coords
end subroutine


subroutine mbd_update_lattice_vectors(calc, latt_vecs)
    type(mbd_calc), intent(inout) :: calc
    real(dp), intent(in) :: latt_vecs(:, :)

    calc%sys%lattice = latt_vecs
end subroutine


subroutine mbd_update_vdw_params_custom(calc, alpha_0, C6, r_vdw, dalpha_0, dC6, dr_vdw)
    type(mbd_calc), intent(inout) :: calc
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    real(dp), intent(in) :: r_vdw(:)
    real(dp), intent(in), optional :: dalpha_0(:, :, :), dC6(:, :, :), dr_vdw(:, :, :)

    calc%alpha_0%val = alpha_0
    calc%C6%val = C6
    calc%damp%r_vdw%val = r_vdw
    if (present(dalpha_0)) calc%alpha_0%dr = dalpha_0
    if (present(dC6)) calc%C6%dr = dC6
    if (present(dr_vdw)) calc%damp%r_vdw%dr = dr_vdw
end subroutine


subroutine mbd_update_vdw_params_from_ratios(calc, ratios, free_values, dratios)
    type(mbd_calc), intent(inout) :: calc
    real(dp), intent(in) :: ratios(:)
    real(dp), intent(in) :: free_values(:, :)
    real(dp), intent(in), optional :: dratios(:, :, :)

    type(vecn) :: vols, ones

    allocate (ones%val(size(ratios)), source=1d0)
    vols%val = ratios
    if (present(dratios)) vols%dr = dratios
    calc%alpha_0 = scale_TS(vecn(free_values(1, :)), vols, ones, 1d0)
    calc%C6 = scale_TS(vecn(free_values(2, :)), vols, ones, 2d0)
    calc%damp%r_vdw = scale_TS(vecn(free_values(3, :)), vols, ones, 1d0/3)
end subroutine


subroutine mbd_get_energy(calc, energy)
    type(mbd_calc), intent(inout) :: calc
    real(dp), intent(out) :: energy

    select case (calc%dispersion_type)
    case ('mbd')
        calc%results = mbd_rsscs_energy(calc%sys, calc%alpha_0, calc%C6, calc%damp)
        energy = calc%results%energy
    case ('ts')
        energy = get_ts_energy(calc%sys, calc%alpha_0%val, calc%C6%val, calc%damp)
    end select
end subroutine


subroutine mbd_get_gradients(calc, gradients)  ! 3 by N  dE/dR
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(out) :: gradients(:, :)

    gradients = transpose(calc%results%forces)
end subroutine


subroutine mbd_get_lattice_derivs(calc, latt_derivs)  ! 3 by 3  (dE/d{abc}_i)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(out) :: latt_derivs(:, :)

    ! TODO
end subroutine


subroutine mbd_get_spectrum_modes(calc, spectrum, modes)
    type(mbd_calc), intent(inout) :: calc
    real(dp), intent(out) :: spectrum(:)
    real(dp), intent(out), optional :: modes(:, :)
    ! TODO document that this can be called only once

    spectrum = calc%results%mode_enes
    if (present(modes)) then
        modes = calc%results%modes
    end if
end subroutine


subroutine mbd_get_damping_parameters(xc, mbd_beta, ts_sr)
    character(len=*), intent(in) :: xc
    real(dp), intent(out) :: mbd_beta, ts_sr

    real(dp) :: d1, d2, d3, d4, d5, d6

    call get_damping_parameters(xc, d1, ts_sr, d2, d3, d4, d5, d6, mbd_beta)
end subroutine


subroutine mbd_get_free_vdw_params(atom_types, table_type, free_values)
    character(len=*), intent(in) :: atom_types(:)  ! e.g. ['Ar', 'Ar']
    character(len=*), intent(in) :: table_type  ! either "ts" or "ts_surf"
    real(dp), intent(out) :: free_values(:, :)  ! 3 by N (alpha_0, C6, R_vdw)

    select case (table_type)
    case ('ts')
        free_values = default_vdw_params(:, species_index(atom_types))
    end select
end subroutine

end module
