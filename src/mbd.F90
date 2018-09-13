! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd

use mbd_constants
use mbd_system_type, only: mbd_system, mbd_calc
use mbd_core, only: mbd_scs_energy, mbd_result, scale_TS
use mbd_damping_type, only: mbd_damping
use mbd_gradients_type, only: mbd_gradients, mbd_grad_switch
use mbd_ts, only: ts_energy
use mbd_common, only: printer
use mbd_vdw_param, only: ts_vdw_params, tssurf_vdw_params, species_index

implicit none

private
public :: mbd_input, mbd_calculation  ! types
public :: mbd_get_free_vdw_params  ! subroutines

type :: mbd_input
    integer :: comm = -1  ! MPI communicator

    ! which calculation will be done (mbd|ts)
    character(len=30) :: dispersion_type = 'mbd'
    logical :: calculate_forces = .true.
    logical :: calculate_spectrum = .false.

    real(dp) :: ts_ene_acc = TS_ENERGY_ACCURACY  ! accuracy of TS energy
    real(dp) :: ts_f_acc = TS_FORCES_ACCURACY  ! accuracy of TS gradients
    integer :: n_omega_grid = N_FREQUENCY_GRID  ! number of frequency grid points
    ! off-gamma shift of k-points in units of inter-k-point distance
    real(dp) :: k_grid_shift = K_GRID_SHIFT

    character(len=20) :: xc = ''
    ! TS damping parameters
    real(dp) :: ts_d = TS_DAMPING_D
    real(dp) :: ts_sr = -1
    ! MBD damping parameters
    real(dp) :: mbd_a = MBD_DAMPING_A
    real(dp) :: mbd_beta = -1

    integer :: k_grid(3)  ! number of k-points along reciprocal axes
    ! is there vacuum along some axes in a periodic calculation
    logical :: vacuum_axis(3) = [.false., .false., .false.]
    character(len=3), allocatable :: atom_types(:)
    character(len=10) :: vdw_params_kind = 'ts'
    real(dp), allocatable :: free_values(:, :)
    logical :: zero_negative_eigvals = .false.

    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: lattice_vectors(:, :)

    character(len=10) :: parallel_mode = 'auto'
end type

type mbd_calculation
    private
    type(mbd_system) :: sys
    type(mbd_damping) :: damp
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
    character(len=30) :: dispersion_type
    type(mbd_calc) :: calc
    type(mbd_result) :: results
    type(mbd_gradients) :: denergy
    logical :: do_gradients
    real(dp), allocatable :: free_values(:, :)
contains
    procedure :: init => mbd_calc_init
    procedure :: destroy => mbd_calc_destroy
    procedure :: update_coords => mbd_calc_update_coords
    procedure :: update_lattice_vectors => mbd_calc_update_lattice_vectors
    procedure :: update_vdw_params_custom => mbd_calc_update_vdw_params_custom
    procedure :: update_vdw_params_from_ratios => &
        mbd_calc_update_vdw_params_from_ratios
    procedure :: get_energy => mbd_calc_get_energy
    procedure :: get_gradients => mbd_calc_get_gradients
    procedure :: get_lattice_derivs => mbd_calc_get_lattice_derivs
    procedure :: get_spectrum_modes => mbd_calc_get_spectrum_modes
    procedure :: get_exception => mbd_calc_get_exception
    procedure :: print_info => mbd_calc_print_info
end type

contains


subroutine mbd_calc_init(this, input)
    class(mbd_calculation), target, intent(inout) :: this
    type(mbd_input), intent(in) :: input

#ifdef WITH_MPI
    if (input%comm /= -1) this%sys%comm = input%comm
#endif
    this%dispersion_type = input%dispersion_type
    this%do_gradients = input%calculate_forces
    if (input%calculate_spectrum) then
        this%sys%get_eigs = .true.
        this%sys%get_modes = .true.
    end if
    this%calc%param%ts_energy_accuracy = input%ts_ene_acc
    ! TODO ... = input%ts_f_acc
    this%calc%param%n_frequency_grid = input%n_omega_grid
    this%calc%param%k_grid_shift = input%k_grid_shift
    this%calc%param%zero_negative_eigs = input%zero_negative_eigvals
    this%sys%k_grid = input%k_grid
    this%sys%vacuum_axis = input%vacuum_axis
    this%sys%coords = input%coords
    if (allocated(input%lattice_vectors)) this%sys%lattice = input%lattice_vectors
    this%sys%parallel_mode = input%parallel_mode
    call this%calc%init_grid()
    call this%sys%init(this%calc)
    if (allocated(input%free_values)) then
        this%free_values = input%free_values
    else
        this%free_values = &
            mbd_get_free_vdw_params(input%atom_types, input%vdw_params_kind)
    end if
    if (input%xc == '') then
        this%damp%beta = input%mbd_beta
        this%damp%a = input%mbd_a
        this%damp%ts_d = input%ts_d
        this%damp%ts_sr = input%ts_sr
    else
        this%calc%exc = this%damp%set_params_from_xc(input%xc, 'MBD@rsSCS')
        if (this%sys%has_exc()) return
        this%calc%exc = this%damp%set_params_from_xc(input%xc, 'TS')
        if (this%sys%has_exc()) return
    end if
end subroutine


subroutine mbd_calc_destroy(this)
    class(mbd_calculation), target, intent(inout) :: this

    deallocate (this%calc%omega_grid, this%calc%omega_grid_w)
    call this%sys%destroy()
end subroutine


subroutine mbd_calc_update_coords(this, coords)
    class(mbd_calculation), intent(inout) :: this
    real(dp), intent(in) :: coords(:, :)

    this%sys%coords = coords
end subroutine


subroutine mbd_calc_update_lattice_vectors(this, latt_vecs)
    class(mbd_calculation), intent(inout) :: this
    real(dp), intent(in) :: latt_vecs(:, :)

    this%sys%lattice = latt_vecs
end subroutine


subroutine mbd_calc_update_vdw_params_custom(this, alpha_0, C6, r_vdw)
    class(mbd_calculation), intent(inout) :: this
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    real(dp), intent(in) :: r_vdw(:)

    this%alpha_0 = alpha_0
    this%C6 = C6
    this%damp%r_vdw = r_vdw
end subroutine


subroutine mbd_calc_update_vdw_params_from_ratios(this, ratios)
    class(mbd_calculation), intent(inout) :: this
    real(dp), intent(in) :: ratios(:)

    real(dp), allocatable :: ones(:)

    allocate (ones(size(ratios)), source=1d0)
    this%alpha_0 = scale_TS(this%free_values(1, :), ratios, ones, 1d0)
    this%C6 = scale_TS(this%free_values(2, :), ratios, ones, 2d0)
    this%damp%r_vdw = scale_TS(this%free_values(3, :), ratios, ones, 1d0/3)
end subroutine


subroutine mbd_calc_get_energy(this, energy)
    class(mbd_calculation), intent(inout) :: this
    real(dp), intent(out) :: energy

    select case (this%dispersion_type)
    case ('mbd')
        this%results = mbd_scs_energy( &
            this%sys, 'rsscs', this%alpha_0, this%C6, this%damp, &
            this%denergy, mbd_grad_switch(dcoords=this%do_gradients) &
        )
        energy = this%results%energy
    case ('ts')
        energy = ts_energy(this%sys, this%alpha_0, this%C6, this%damp)
    end select
end subroutine


subroutine mbd_calc_get_gradients(this, gradients)  ! 3 by N  dE/dR
    class(mbd_calculation), intent(in) :: this
    real(dp), intent(out) :: gradients(:, :)

    gradients = transpose(this%denergy%dcoords)
end subroutine


subroutine mbd_calc_get_lattice_derivs(this, latt_derivs)  ! 3 by 3  (dE/d{abc}_i)
    class(mbd_calculation), intent(in) :: this
    real(dp), intent(out) :: latt_derivs(:, :)

    ! TODO
end subroutine


subroutine mbd_calc_get_spectrum_modes(this, spectrum, modes)
    class(mbd_calculation), intent(inout) :: this
    real(dp), intent(out) :: spectrum(:)
    real(dp), intent(out), optional :: modes(:, :)
    ! TODO document that this can be called only once

    spectrum = this%results%mode_eigs
    if (present(modes)) then
        modes = this%results%modes
    end if
end subroutine


subroutine mbd_calc_get_exception(this, code, origin, msg)
    class(mbd_calculation), intent(inout) :: this
    integer, intent(out) :: code
    character(*), intent(out) :: origin
    character(*), intent(out) :: msg

    code = this%calc%exc%code
    if (code == 0) return
    origin = this%calc%exc%origin
    msg = this%calc%exc%msg
    this%calc%exc%code = 0
    this%calc%exc%origin = ''
    this%calc%exc%msg = ''
end subroutine


subroutine mbd_calc_print_info(this, info)
    class(mbd_calculation), intent(inout) :: this
    procedure(printer) :: info

    call this%calc%info%print(info)
end subroutine


function mbd_get_free_vdw_params(atom_types, table_type) result(free_values)
    character(len=*), intent(in) :: atom_types(:)  ! e.g. ['Ar', 'Ar']
    character(len=*), intent(in) :: table_type  ! either "ts" or "ts_surf"
    ! 3 by N (alpha_0, C6, R_vdw)
    real(dp) :: free_values(3, size(atom_types))

    select case (table_type)
    case ('ts')
        free_values = ts_vdw_params(:, species_index(atom_types))
    case ('tssurf')
        free_values = tssurf_vdw_params(:, species_index(atom_types))
    end select
end function

end module
