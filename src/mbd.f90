! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd

use mbd_build_flags, only: WITH_MPI
use mbd_mpi, only: sync_sum, broadcast, MPI_COMM_WORLD
use mbd_common, only: tostr, nan, print_matrix, dp, pi, printer, exception
use mbd_linalg, only: &
    operator(.cprod.), diag, invert, diagonalize, sdiagonalize, diagonalized, &
    sdiagonalized, inverted, sinvert, add_diag, repeatn, mult_cprod
use mbd_types, only: mat3n3n, mat33, scalar, vecn

implicit none

private
public :: mbd_param, mbd_calc, mbd_damping, mbd_work, mbd_system, &
    init_grid, get_mbd_energy, dipole_matrix, mbd_rsscs_energy, mbd_scs_energy, &
    run_tests, get_sigma_selfint
public :: get_ts_energy, init_eqi_grid, eval_mbd_nonint_density, &
    eval_mbd_int_density, nbody_coeffs, get_damping_parameters, v_to_r, &
    clock_rate

interface operator(.prod.)
    module procedure T_damped__
end interface

real(dp), parameter :: ang = 1.8897259886d0
integer, parameter :: n_timestamps = 100

type :: mbd_param
    real(dp) :: ts_energy_accuracy = 1.d-10
    real(dp) :: ts_cutoff_radius = 50.d0*ang
    real(dp) :: dipole_low_dim_cutoff = 100.d0*ang
    real(dp) :: dipole_cutoff = 400.d0*ang  ! used only when Ewald is off
    real(dp) :: ewald_real_cutoff_scaling = 1.d0
    real(dp) :: ewald_rec_cutoff_scaling = 1.d0
    real(dp) :: k_grid_shift = 0.5d0
    logical :: ewald_on = .true.
    logical :: zero_negative_eigs = .false.
    integer :: mbd_nbody_max = 3
    integer :: rpa_order_max = 10
    integer :: n_frequency_grid = 15
end type

type :: mbd_timing
    logical :: measure_time = .true.
    integer :: timestamps(n_timestamps), ts_counts(n_timestamps)
    integer :: ts_cnt, ts_rate, ts_cnt_max, ts_aid
end type mbd_timing

type :: mbd_calc
    type(mbd_param) :: param
    type(mbd_timing) :: tm
    real(dp), allocatable :: omega_grid(:)
    real(dp), allocatable :: omega_grid_w(:)
    logical :: parallel = .false.
    integer :: comm = MPI_COMM_WORLD
    integer :: io = -1
    integer :: my_task = 0
    integer :: n_tasks = 1
    logical :: mute = .false.
end type mbd_calc

type :: mbd_damping
    character(len=20) :: version
    real(dp) :: beta = 0.d0
    real(dp) :: a = 6.d0
    real(dp) :: ts_d = 20.d0
    real(dp) :: ts_sr = 0.d0
    real(dp) :: mayer_scaling = 1.d0
    real(dp), allocatable :: r_vdw(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: damping_custom(:, :)
    real(dp), allocatable :: potential_custom(:, :, :, :)
end type mbd_damping

type :: mbd_work
    logical :: get_eigs = .false.
    logical :: get_modes = .false.
    logical :: get_rpa_orders = .false.
    integer :: i_kpt = 0
    real(dp), allocatable :: k_pts(:, :)
    real(dp), allocatable :: mode_enes(:)
    real(dp), allocatable :: modes(:, :)
    real(dp), allocatable :: rpa_orders(:)
    real(dp), allocatable :: mode_enes_k(:, :)
    complex(dp), allocatable :: modes_k(:, :, :)
    real(dp), allocatable :: rpa_orders_k(:, :)
    real(dp), allocatable :: forces(:, :)
    type(exception) :: exc
end type

type :: mbd_system
    type(mbd_calc), pointer :: calc
    type(mbd_work) :: work
    real(dp), allocatable :: coords(:, :)
    logical :: periodic = .false.
    logical :: vacuum_axis(3) = (/ .false., .false., .false. /)
    real(dp) :: lattice(3, 3)
    integer :: k_grid(3)
    integer :: supercell(3)
    logical :: do_rpa = .false.
    logical :: do_reciprocal = .true.
    logical :: do_force = .false.
end type mbd_system

contains


real(dp) function mbd_rsscs_energy(sys, alpha_0, C6, damp)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp

    real(dp), allocatable :: alpha_dyn(:, :)
    type(vecn) :: alpha_dyn_rsscs(0:ubound(sys%calc%omega_grid, 1))
    real(dp), allocatable :: C6_rsscs(:)
    real(dp), allocatable :: R_vdw_rsscs(:)
    type(mbd_damping) :: damp_rsscs, damp_mbd
    integer :: n_freq

    n_freq = ubound(sys%calc%omega_grid, 1)
    allocate (alpha_dyn(0:n_freq, size(sys%coords, 1)))
    alpha_dyn = alpha_dynamic_ts(sys%calc, alpha_0, C6)
    damp_rsscs = damp
    damp_rsscs%version = 'fermi,dip,gg'
    alpha_dyn_rsscs = run_scs(sys, alpha_dyn, damp_rsscs)
    C6_rsscs = get_C6_from_alpha(sys%calc, alpha_dyn_rsscs)
    R_vdw_rsscs = damp%R_vdw*(alpha_dyn_rsscs(0)%val/alpha_dyn(0, :))**(1.d0/3)
    damp_mbd%version = 'fermi,dip'
    damp_mbd%r_vdw = R_vdw_rsscs
    damp_mbd%beta = damp%beta
    mbd_rsscs_energy = get_mbd_energy(sys, alpha_dyn_rsscs(0)%val, C6_rsscs, damp_mbd)
    if (has_exc(sys)) call print_exc(sys)
end function mbd_rsscs_energy


real(dp) function mbd_scs_energy(sys, alpha_0, C6, damp)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp

    real(dp), allocatable :: alpha_dyn(:, :)
    type(vecn) :: alpha_dyn_scs(0:ubound(sys%calc%omega_grid, 1))
    real(dp), allocatable :: C6_scs(:)
    real(dp), allocatable :: R_vdw_scs(:)
    type(mbd_damping) :: damp_scs, damp_mbd
    integer :: n_freq

    n_freq = ubound(sys%calc%omega_grid, 1)
    allocate (alpha_dyn(0:n_freq, size(sys%coords, 1)))
    alpha_dyn = alpha_dynamic_ts(sys%calc, alpha_0, C6)
    damp_scs = damp
    damp_scs%version = 'dip,gg'
    alpha_dyn_scs = run_scs(sys, alpha_dyn, damp_scs)
    C6_scs = get_C6_from_alpha(sys%calc, alpha_dyn_scs)
    R_vdw_scs = damp%R_vdw*(alpha_dyn_scs(0)%val(:)/alpha_dyn(0, :))**(1.d0/3)
    damp_mbd%version = 'dip,1mexp'
    damp_mbd%r_vdw = R_vdw_scs
    damp_mbd%beta = 1.d0
    damp_mbd%a = damp%a
    mbd_scs_energy = get_mbd_energy(sys, alpha_dyn_scs(0)%val, C6_scs, damp_mbd)
    if (has_exc(sys)) call print_exc(sys)
end function mbd_scs_energy


function get_ts_energy(sys, alpha_0, C6, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp
    real(dp) :: ene

    real(dp) :: C6_ij, r(3), r_norm, R_vdw_ij, ene_shell, ene_pair, R_cell(3)
    type(scalar) :: f_damp
    integer :: i_shell, i_cell, i_atom, j_atom, range_cell(3), idx_cell(3)
    real(dp), parameter :: shell_thickness = 10.d0
    logical :: is_crystal, is_parallel

    is_crystal = sys%periodic
    is_parallel = sys%calc%parallel

    ene = 0.d0
    i_shell = 0
    do
        i_shell = i_shell+1
        ene_shell = 0.d0
        if (is_crystal) then
            range_cell = supercell_circum(sys, sys%lattice, i_shell*shell_thickness)
        else
            range_cell = (/ 0, 0, 0 /)
        end if
        idx_cell = (/ 0, 0, -1 /)
        do i_cell = 1, product(1+2*range_cell)
            call shift_cell(idx_cell, -range_cell, range_cell)
            ! MPI code begin
            if (is_parallel .and. is_crystal) then
                if (sys%calc%my_task /= modulo(i_cell, sys%calc%n_tasks)) cycle
            end if
            ! MPI code end
            if (is_crystal) then
                R_cell = matmul(idx_cell, sys%lattice)
            else
                R_cell = (/ 0.d0, 0.d0, 0.d0 /)
            end if
            do i_atom = 1, size(sys%coords, 1)
                ! MPI code begin
                if (is_parallel .and. .not. is_crystal) then
                    if (sys%calc%my_task /= modulo(i_atom, sys%calc%n_tasks)) cycle
                end if
                ! MPI code end
                do j_atom = 1, i_atom
                    if (i_cell == 1) then
                        if (i_atom == j_atom) cycle
                    end if
                    r = sys%coords(i_atom, :)-sys%coords(j_atom, :)-R_cell
                    r_norm = sqrt(sum(r**2))
                    if (r_norm > sys%calc%param%ts_cutoff_radius) cycle
                    if (is_crystal) then
                        if (r_norm >= i_shell*shell_thickness &
                            .or. r_norm < (i_shell-1)*shell_thickness) then
                            cycle
                        end if
                    end if
                    C6_ij = combine_C6( &
                        C6(i_atom), C6(j_atom), &
                        alpha_0(i_atom), alpha_0(j_atom))
                    if (allocated(damp%r_vdw)) then
                        R_vdw_ij = damp%r_vdw(i_atom)+damp%r_vdw(j_atom)
                    end if
                    select case (damp%version)
                        case ("fermi")
                            f_damp = damping_fermi(r, damp%ts_sr*R_vdw_ij, damp%ts_d, .false.)
                        case ("fermi2")
                            f_damp = damping_fermi(r, damp%ts_sr*R_vdw_ij, damp%ts_d, .false.)
                            f_damp%val = f_damp%val**2
                        case ("custom")
                            f_damp%val = damp%damping_custom(i_atom, j_atom)
                    end select
                    ene_pair = -C6_ij*f_damp%val/r_norm**6
                    if (i_atom == j_atom) then
                        ene_shell = ene_shell+ene_pair/2
                    else
                        ene_shell = ene_shell+ene_pair
                    endif
                end do ! j_atom
            end do ! i_atom
        end do ! i_cell
        ! MPI code begin
        if (is_parallel) then
            call sync_sum(ene_shell, sys%calc%comm)
        end if
        ! MPI code end
        ene = ene+ene_shell
        if (.not. is_crystal) exit
        if (i_shell > 1 .and. abs(ene_shell) < sys%calc%param%ts_energy_accuracy) then
            call printer(sys%calc%io, "Periodic TS converged in " &
                //trim(tostr(i_shell))//" shells, " &
                //trim(tostr(i_shell*shell_thickness/ang))//" angstroms")
            exit
        endif
    end do ! i_shell
end function get_ts_energy


type(mat3n3n) function dipole_matrix(sys, damp, k_point) result(dipmat)
    type(mbd_system), intent(inout) :: sys
    type(mbd_damping), intent(in) :: damp
    real(dp), intent(in), optional :: k_point(3)

    real(dp) :: R_cell(3), r(3), r_norm, R_vdw_ij, &
        sigma_ij, volume, ewald_alpha, real_space_cutoff, f_ij
    type(mat33) :: Tpp
    complex(dp) :: Tpp_c(3, 3)
    character(len=1) :: parallel_mode
    integer :: i_atom, j_atom, i_cell, idx_cell(3), range_cell(3), i, j, n_atoms
    logical :: mute, do_ewald

    do_ewald = .false.
    mute = sys%calc%mute
    n_atoms = size(sys%coords, 1)
    if (sys%calc%parallel) then
        parallel_mode = 'A' ! atoms
        if (sys%periodic .and. n_atoms < sys%calc%n_tasks) then
            parallel_mode = 'C' ! cells
        end if
    else
        parallel_mode = ''
    end if

    if (present(k_point)) then
        allocate (dipmat%cplx(3*n_atoms, 3*n_atoms), source=(0.d0, 0.d0))
    else
        allocate (dipmat%re(3*n_atoms, 3*n_atoms), source=0.d0)
        if (sys%do_force) then
            allocate (dipmat%re_dr(3*n_atoms, 3*n_atoms, 3), source=0.d0)
        end if
    end if
    ! MPI code end
    if (sys%periodic) then
        if (any(sys%vacuum_axis)) then
            real_space_cutoff = sys%calc%param%dipole_low_dim_cutoff
        else if (sys%calc%param%ewald_on) then
            do_ewald = .true.
            volume = max(abs(dble(product(diagonalized(sys%lattice)))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1.d0/3)
            real_space_cutoff = 6.d0/ewald_alpha*sys%calc%param%ewald_real_cutoff_scaling
            if (.not. mute) call printer(sys%calc%io, &
                'Ewald: using alpha = '//trim(tostr(ewald_alpha)) &
                //', real cutoff = '//trim(tostr(real_space_cutoff)))
        else
            real_space_cutoff = sys%calc%param%dipole_cutoff
        end if
        range_cell = supercell_circum(sys, sys%lattice, real_space_cutoff)
    else
        range_cell(:) = 0
    end if
    if (sys%periodic) then
        if (.not. mute) call printer(sys%calc%io, &
            'Ewald: summing real part in cell vector range of ' &
            //trim(tostr(1+2*range_cell(1)))//'x' &
            //trim(tostr(1+2*range_cell(2)))//'x' &
            //trim(tostr(1+2*range_cell(3))))
    end if
    call ts(sys%calc, 11)
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, product(1+2*range_cell)
        call shift_cell(idx_cell, -range_cell, range_cell)
        ! MPI code begin
        if (parallel_mode == 'C') then
            if (sys%calc%my_task /= modulo(i_cell, sys%calc%n_tasks)) cycle
        end if
        ! MPI code end
        if (sys%periodic) then
            R_cell = matmul(idx_cell, sys%lattice)
        else
            R_cell(:) = 0.d0
        end if
        do i_atom = 1, n_atoms
            ! MPI code begin
            if (parallel_mode == 'A') then
                if (sys%calc%my_task /= modulo(i_atom, sys%calc%n_tasks)) cycle
            end if
            ! MPI code end
            !$omp parallel do private(r, r_norm, R_vdw_ij, sigma_ij, overlap_ij, C6_ij, &
            !$omp    Tpp, i, j, Tpp_c)
            do j_atom = i_atom, n_atoms
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                r = sys%coords(i_atom, :)-sys%coords(j_atom, :)-R_cell
                r_norm = sqrt(sum(r**2))
                if (sys%periodic .and. r_norm > real_space_cutoff) cycle
                if (allocated(damp%R_vdw)) then
                    R_vdw_ij = damp%R_vdw(i_atom)+damp%R_vdw(j_atom)
                end if
                if (allocated(damp%sigma)) then
                    sigma_ij = damp%mayer_scaling*sqrt(sum(damp%sigma([i_atom, j_atom])**2))
                end if
                select case (damp%version)
                    case ("bare")
                        Tpp = T_bare_v2(r, sys%do_force)
                    case ("dip,1mexp")
                        Tpp%val = T_1mexp_coulomb(r, damp%beta*R_vdw_ij, damp%a)
                    case ("fermi,dip")
                        Tpp = damping_fermi( &
                            r, damp%beta*R_vdw_ij, damp%a, sys%do_force &
                        ).prod.T_bare_v2(r, sys%do_force)
                    case ("sqrtfermi,dip")
                        Tpp = damping_sqrtfermi( &
                            r, damp%beta*R_vdw_ij, damp%a, sys%do_force &
                        ).prod.T_bare_v2(r, sys%do_force)
                    case ("custom,dip")
                        Tpp%val = damp%damping_custom(i_atom, j_atom)*T_bare(r)
                    case ("dip,custom")
                        Tpp%val = damp%potential_custom(i_atom, j_atom, :, :)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma_ij, sys%do_force)
                    case ("fermi,dip,gg")
                        Tpp = op1minus(damping_fermi( &
                            r, damp%beta*R_vdw_ij, damp%a, sys%do_force &
                        )).prod.T_erf_coulomb(r, sigma_ij, sys%do_force)
                        do_ewald = .false.
                    case ("sqrtfermi,dip,gg")
                        Tpp = op1minus(damping_sqrtfermi( &
                            r, damp%beta*R_vdw_ij, damp%a, sys%do_force &
                        )).prod.T_erf_coulomb(r, sigma_ij, sys%do_force)
                        do_ewald = .false.
                    case ("custom,dip,gg")
                        f_ij = 1.d0-damp%damping_custom(i_atom, j_atom)
                        Tpp = T_erf_coulomb(r, sigma_ij, sys%do_force)
                        Tpp%val = f_ij*Tpp%val
                        do_ewald = .false.
                end select
                if (do_ewald) then
                    Tpp%val = Tpp%val+T_erfc(r, ewald_alpha)-T_bare(r)
                end if
                if (present(k_point)) then
                    Tpp_c = Tpp%val*exp(-cmplx(0.d0, 1.d0, 8)*( &
                        dot_product(k_point, r)))
                end if
                i = 3*(i_atom-1)
                j = 3*(j_atom-1)
                if (present(k_point)) then
                    associate (T => dipmat%cplx(i+1:i+3, j+1:j+3))
                        T = T + Tpp_c
                    end associate
                else
                    associate (T => dipmat%re(i+1:i+3, j+1:j+3))
                        T = T + Tpp%val
                    end associate
                    if (sys%do_force) then
                        associate (T => dipmat%re_dr(i+1:i+3, j+1:j+3, :))
                            T = T + Tpp%dr
                        end associate
                    end if
                end if
            end do ! j_atom
            !$omp end parallel do
        end do ! i_atom
    end do ! i_cell
    call ts(sys%calc, -11)
    ! MPI code begin
    if (sys%calc%parallel) then
        if (present(k_point)) then
            call sync_sum(dipmat%cplx, sys%calc%comm)
        else
            call sync_sum(dipmat%re, sys%calc%comm)
        end if
    end if
    ! MPI code end
    if (do_ewald) then
        call add_ewald_dipole_parts(sys, ewald_alpha, dipmat, k_point)
    end if
    if (present(k_point)) then
        do i_atom = 1, 3*n_atoms
            do j_atom = i_atom+1, 3*n_atoms
                dipmat%cplx(j_atom, i_atom) = conjg(dipmat%cplx(i_atom, j_atom))
            end do
        end do
    else
        do i_atom = 1, 3*n_atoms
            do j_atom = i_atom+1, 3*n_atoms
                dipmat%re(j_atom, i_atom) = dipmat%re(i_atom, j_atom)
            end do
        end do
    end if
end function dipole_matrix


subroutine add_ewald_dipole_parts(sys, alpha, dipmat, k_point)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha
    real(dp), intent(in), optional :: k_point(3)
    type(mat3n3n), intent(inout) :: dipmat

    logical :: is_parallel, mute, do_surface
    real(dp) :: rec_unit_cell(3, 3), volume, G_vector(3), r(3), k_total(3), &
        k_sq, rec_space_cutoff, Tpp(3, 3), k_prefactor(3, 3), elem
    complex(dp) :: Tpp_c(3, 3)
    integer :: &
        i_atom, j_atom, i, j, i_xyz, j_xyz, idx_G_vector(3), i_G_vector, &
        range_G_vector(3)
    character(len=1) :: parallel_mode

    is_parallel = sys%calc%parallel
    mute = sys%calc%mute
    if (is_parallel) then
        parallel_mode = 'A' ! atoms
        if (size(sys%coords, 1) < sys%calc%n_tasks) then
            parallel_mode = 'G' ! G vectors
        end if
    else
        parallel_mode = ''
    end if

    ! MPI code begin
    if (is_parallel) then
        ! will be restored by syncing at the end
        if (present(k_point)) then
            dipmat%cplx = dipmat%cplx/sys%calc%n_tasks
        else
            dipmat%re = dipmat%re/sys%calc%n_tasks
        end if
    end if
    ! MPI code end
    rec_unit_cell = 2*pi*inverted(transpose(sys%lattice))
    volume = abs(dble(product(diagonalized(sys%lattice))))
    rec_space_cutoff = 10.d0*alpha*sys%calc%param%ewald_rec_cutoff_scaling
    range_G_vector = supercell_circum(sys, rec_unit_cell, rec_space_cutoff)
    if (.not. mute) then
        call printer(sys%calc%io, 'Ewald: using reciprocal cutoff = ' &
            //trim(tostr(rec_space_cutoff)))
        call printer(sys%calc%io, 'Ewald: summing reciprocal part in G vector range of ' &
            //trim(tostr(1+2*range_G_vector(1)))//'x' &
            //trim(tostr(1+2*range_G_vector(2)))//'x' &
            //trim(tostr(1+2*range_G_vector(3))))
    end if
    call ts(sys%calc, 12)
    idx_G_vector = (/ 0, 0, -1 /)
    do i_G_vector = 1, product(1+2*range_G_vector)
        call shift_cell(idx_G_vector, -range_G_vector, range_G_vector)
        if (i_G_vector == 1) cycle
        ! MPI code begin
        if (parallel_mode == 'G') then
            if (sys%calc%my_task /= modulo(i_G_vector, sys%calc%n_tasks)) cycle
        end if
        ! MPI code end
        G_vector = matmul(idx_G_vector, rec_unit_cell)
        if (present(k_point)) then
            k_total = k_point+G_vector
        else
            k_total = G_vector
        end if
        k_sq = sum(k_total**2)
        if (sqrt(k_sq) > rec_space_cutoff) cycle
        k_prefactor(:, :) = 4*pi/volume*exp(-k_sq/(4*alpha**2))
        forall (i_xyz = 1:3, j_xyz = 1:3) &
                k_prefactor(i_xyz, j_xyz) = k_prefactor(i_xyz, j_xyz) &
                *k_total(i_xyz)*k_total(j_xyz)/k_sq
        do i_atom = 1, size(sys%coords, 1)
            ! MPI code begin
            if (parallel_mode == 'A') then
                if (sys%calc%my_task /= modulo(i_atom, sys%calc%n_tasks)) cycle
            end if
            ! MPI code end
            !$omp parallel do private(r, Tpp, i, j, Tpp_c)
            do j_atom = i_atom, size(sys%coords, 1)
                r = sys%coords(i_atom, :)-sys%coords(j_atom, :)
                if (present(k_point)) then
                    Tpp_c = k_prefactor*exp(cmplx(0.d0, 1.d0, 8) &
                        *dot_product(G_vector, r))
                else
                    Tpp = k_prefactor*cos(dot_product(G_vector, r))
                end if
                i = 3*(i_atom-1)
                j = 3*(j_atom-1)
                if (present(k_point)) then
                    associate (T => dipmat%cplx(i+1:i+3, j+1:j+3))
                        T = T + Tpp_c
                    end associate
                else
                    associate (T => dipmat%re(i+1:i+3, j+1:j+3))
                        T = T + Tpp
                    end associate
                end if
            end do ! j_atom
            !$omp end parallel do
        end do ! i_atom
    end do ! i_G_vector
    ! MPI code begin
    if (is_parallel) then
        if (present(k_point)) then
            call sync_sum(dipmat%cplx, sys%calc%comm)
        else
            call sync_sum(dipmat%re, sys%calc%comm)
        end if
    end if
    ! MPI code end
    call add_diag(dipmat, -4*alpha**3/(3*sqrt(pi))) ! self energy
    do_surface = .true.
    if (present(k_point)) then
        k_sq = sum(k_point**2)
        if (sqrt(k_sq) > 1.d-15) then
            do_surface = .false.
            do i_atom = 1, size(sys%coords, 1)
            do j_atom = i_atom, size(sys%coords, 1)
                do i_xyz = 1, 3
                do j_xyz = 1, 3
                    i = 3*(i_atom-1)+i_xyz
                    j = 3*(j_atom-1)+j_xyz
                    elem = 4*pi/volume*k_point(i_xyz)*k_point(j_xyz)/k_sq &
                        *exp(-k_sq/(4*alpha**2))
                    if (present(k_point)) then
                        dipmat%cplx(i, j) = dipmat%cplx(i, j) + elem
                    else
                        dipmat%re(i, j) = dipmat%re(i, j) + elem
                    end if ! present(k_point)
                end do ! j_xyz
                end do ! i_xyz
            end do ! j_atom
            end do ! i_atom
        end if ! k_sq >
    end if ! k_point present
    if (do_surface) then ! surface energy
        do i_atom = 1, size(sys%coords, 1)
        do j_atom = i_atom, size(sys%coords, 1)
            do i_xyz = 1, 3
                i = 3*(i_atom-1)+i_xyz
                j = 3*(j_atom-1)+i_xyz
                if (present(k_point)) then
                    dipmat%cplx(i, j) = dipmat%cplx(i, j) + 4*pi/(3*volume)
                else
                    dipmat%re(i, j) = dipmat%re(i, j) + 4*pi/(3*volume)
                end if
            end do ! i_xyz
        end do ! j_atom
        end do ! i_atom
    end if
    call ts(sys%calc, -12)
end subroutine


subroutine init_grid(calc)
    type(mbd_calc), intent(inout) :: calc

    integer :: n

    n = calc%param%n_frequency_grid
    allocate (calc%omega_grid(0:n))
    allocate (calc%omega_grid_w(0:n))
    calc%omega_grid(0) = 0.d0
    calc%omega_grid_w(0) = 0.d0
    call get_omega_grid(n, 0.6d0, calc%omega_grid(1:n), calc%omega_grid_w(1:n))
    call printer(calc%io, &
        "Initialized a radial integration grid of "//trim(tostr(n))//" points." &
    )
    call printer(calc%io, &
        "Relative quadrature error in C6 of carbon atom: "// &
        trim(tostr(test_frequency_grid(calc))) &
    )
end subroutine


real(dp) function test_frequency_grid(calc) result(error)
    type(mbd_calc), intent(in) :: calc

    real(dp) :: alpha(0:ubound(calc%omega_grid, 1), 1)

    alpha = alpha_dynamic_ts(calc, (/ 21.d0 /), (/ 99.5d0 /))
    error = abs(get_total_C6_from_alpha(calc, alpha)/99.5d0-1.d0)
end function


subroutine get_omega_grid(n, L, x, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: L
    real(dp), intent(out) :: x(n), w(n)

    call gauss_legendre(n, x, w)
    w = 2*L/(1-x)**2*w
    x = L*(1+x)/(1-x)
    w = w(n:1:-1)
    x = x(n:1:-1)
end subroutine get_omega_grid


subroutine gauss_legendre(n, r, w)
    use mbd_build_flags, only: LEGENDRE_PREC

    integer, intent(in) :: n
    real(dp), intent(out) :: r(n), w(n)

    integer, parameter :: q = LEGENDRE_PREC
    integer, parameter :: n_iter = 1000
    real(q) :: x, f, df, dx
    integer :: k, iter, i
    real(q) :: Pk(0:n), Pk1(0:n-1), Pk2(0:n-2)

    if (n == 1) then
        r(1) = 0.d0
        w(1) = 2.d0
        return
    end if
    Pk2(0) = 1._q  ! k = 0
    Pk1(0:1) = (/ 0._q, 1._q /)  ! k = 1
    do k = 2, n
        Pk(0:k) = ((2*k-1)*(/ 0.0_q, Pk1(0:k-1) /)-(k-1)*(/ Pk2(0:k-2), 0._q, 0._q /))/k
        if (k < n) then
            Pk2(0:k-1) = Pk1(0:k-1)
            Pk1(0:k) = Pk(0:k)
        end if
    end do
    ! now Pk contains k-th Legendre polynomial
    do i = 1, n
        x = cos(pi*(i-0.25_q)/(n+0.5_q))
        do iter = 1, n_iter
            df = 0._q
            f = Pk(n)
            do k = n-1, 0, -1
                df = f + x*df
                f = Pk(k) + x*f
            end do
            dx = f/df
            x = x-dx
            if (abs(dx) < 10*epsilon(dx)) exit
        end do
        r(i) = dble(x)
        w(i) = dble(2/((1-x**2)*df**2))
    end do
end subroutine


subroutine init_eqi_grid(calc, n, a, b)
    type(mbd_calc), intent(inout) :: calc
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b

    real(dp) :: delta
    integer :: i

    if (allocated(calc%omega_grid)) deallocate(calc%omega_grid)
    if (allocated(calc%omega_grid_w)) deallocate(calc%omega_grid_w)
    allocate (calc%omega_grid(0:n))
    allocate (calc%omega_grid_w(0:n))
    calc%omega_grid(0) = 0.d0
    calc%omega_grid_w(0) = 0.d0
    delta = (b-a)/n
    calc%omega_grid(1:n) = (/ (a+delta/2+i*delta, i = 0, n-1) /)
    calc%omega_grid_w(1:n) = delta
end subroutine


function run_scs(sys, alpha, damp) result(alpha_scs)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha(0:, :)
    type(mbd_damping), intent(in) :: damp
    type(vecn) :: alpha_scs(0:ubound(alpha, 1))

    type(mat3n3n) :: alpha_full, T
    integer :: i_grid_omega
    logical :: mute
    type(mbd_damping) :: damp_local

    mute = sys%calc%mute

    do i_grid_omega = 0, ubound(sys%calc%omega_grid, 1)
        damp_local = damp
        damp_local%sigma = get_sigma_selfint(alpha(i_grid_omega, :))
        T = dipole_matrix(sys, damp_local)
        if (sys%do_force) then
            alpha_full = T
        else
            call move_alloc(T%re, alpha_full%re)
        end if
        call add_diag(alpha_full, repeatn(1.d0/alpha(i_grid_omega, :), 3))
        call ts(sys%calc, 32)
        call sinvert(alpha_full%re, sys%work%exc)
        if (has_exc(sys)) return
        call ts(sys%calc, -32)
        alpha_scs(i_grid_omega)%val = contract_polarizability(alpha_full%re)
        sys%calc%mute = .true.
    end do

    sys%calc%mute = mute
end function run_scs


function get_mbd_energy(sys, alpha_0, C6, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp
    real(dp) :: ene

    logical :: is_parallel, do_rpa, is_reciprocal, is_crystal
    real(dp), allocatable :: alpha(:, :)

    is_parallel = sys%calc%parallel
    is_crystal = sys%periodic
    do_rpa = sys%do_rpa
    is_reciprocal = sys%do_reciprocal
    if (.not. is_crystal) then
        if (.not. do_rpa) then
            ene = get_single_mbd_energy(sys, alpha_0, C6, damp)
        else
            allocate (alpha(0:ubound(sys%calc%omega_grid, 1), size(alpha_0)))
            alpha = alpha_dynamic_ts(sys%calc, alpha_0, C6)
            ene = get_single_rpa_energy(sys, alpha, damp)
            deallocate (alpha)
        end if
    else
        if (is_reciprocal) then
            ene = get_reciprocal_mbd_energy(sys, alpha_0, C6, damp)
        else
            ene = get_supercell_mbd_energy(sys, alpha_0, C6, damp)
        end if
    end if
end function get_mbd_energy


real(dp) function get_supercell_mbd_energy(sys, alpha_0, C6, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp

    logical :: do_rpa
    real(dp) :: R_cell(3)
    integer :: i_atom, i
    integer :: i_cell
    integer :: idx_cell(3), n_cells

    real(dp), allocatable :: &
        xyz_super(:, :), alpha_0_super(:), C6_super(:), &
        R_vdw_super(:), alpha_ts_super(:, :)
    type(mbd_system) :: sys_super
    type(mbd_damping) :: damp_super

    do_rpa = sys%do_rpa

    sys_super%calc = sys%calc
    sys_super%work = sys%work
    n_cells = product(sys%supercell)
    do i = 1, 3
        sys_super%lattice(i, :) = sys%lattice(i, :)*sys%supercell(i)
    end do
    allocate (sys_super%coords(n_cells*size(sys%coords, 1), 3))
    allocate (alpha_0_super(n_cells*size(alpha_0)))
    allocate (alpha_ts_super(0:ubound(sys%calc%omega_grid, 1), n_cells*size(alpha_0)))
    allocate (C6_super(n_cells*size(C6)))
    if (allocated(damp%r_vdw)) allocate (damp_super%r_vdw(n_cells*size(damp%r_vdw)))
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, n_cells
        call shift_cell(idx_cell, (/ 0, 0, 0 /), sys%supercell-1)
        R_cell = matmul(idx_cell, sys%lattice)
        do i_atom = 1, size(sys%coords, 1)
            i = (i_cell-1)*size(sys%coords, 1)+i_atom
            sys_super%coords(i, :) = sys%coords(i_atom, :)+R_cell
            alpha_0_super(i) = alpha_0(i_atom)
            C6_super(i) = C6(i_atom)
            if (allocated(damp%R_vdw)) then
                damp_super%R_vdw(i) = damp%R_vdw(i_atom)
            end if
        end do
    end do
    if (do_rpa) then
        alpha_ts_super = alpha_dynamic_ts(sys%calc, alpha_0_super, C6_super)
        ene = get_single_rpa_energy( &
            sys_super, alpha_ts_super, damp_super &
        )
    else
        ene = get_single_mbd_energy( &
            sys_super, alpha_0_super, C6_super, damp_super &
        )
    end if
    deallocate (xyz_super)
    deallocate (alpha_0_super)
    deallocate (alpha_ts_super)
    deallocate (C6_super)
    deallocate (R_vdw_super)
    ene = ene/n_cells
    if (sys%work%get_rpa_orders) then
        sys%work%rpa_orders =sys_super%work%rpa_orders/n_cells
    end if
end function get_supercell_mbd_energy
    

real(dp) function get_single_mbd_energy(sys, alpha_0, C6, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp

    type(mat3n3n) :: relay
    real(dp), allocatable :: eigs(:)
    real(dp), allocatable :: omega(:)
    integer :: i_xyz, i
    integer :: n_negative_eigs, n_atoms
    logical :: is_parallel
    real(dp), allocatable :: c_lambda12i_c(:, :)
    real(dp), allocatable :: c_lambda14i(:, :)

    is_parallel = sys%calc%parallel

    n_atoms = size(sys%coords, 1)
    allocate (eigs(3*n_atoms))
    ! relay%re = T
    relay = dipole_matrix(sys, damp)
    omega = omega_eff(C6, alpha_0)
    call form_mbd_matrix(relay, alpha_0, omega)
    call ts(sys%calc, 21)
    if (.not. is_parallel .or. sys%calc%my_task == 0) then
        if (sys%work%get_modes .or. sys%do_force) then
            call sdiagonalize('V', relay%re, eigs, sys%work%exc)
            call move_alloc(relay%re, sys%work%modes)
        else
            call sdiagonalize('N', relay%re, eigs, sys%work%exc)
        end if
        if (has_exc(sys)) return
    end if
    ! MPI code begin
    if (is_parallel) then
        call broadcast(relay%re, sys%calc%comm)
        call broadcast(eigs, sys%calc%comm)
    end if
    ! MPI code end
    call ts(sys%calc, -21)
    if (sys%work%get_eigs) then
        sys%work%mode_enes = sqrt(eigs)
        where (eigs < 0) sys%work%mode_enes = 0.d0
    end if
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        call printer(sys%calc%io, &
            "CDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues" &
        )
        if (sys%calc%param%zero_negative_eigs) where (eigs < 0) eigs = 0.d0
    end if
    ene = 1.d0/2*sum(sqrt(eigs))-3.d0/2*sum(omega)
    if (sys%do_force) then
        allocate (c_lambda14i(3*n_atoms, 3*n_atoms))
        allocate (sys%work%forces(n_atoms, 3))
        forall (i = 1:3*n_atoms)
            c_lambda14i(:, i) = eigs(i)**(-1.d0/4)*sys%work%modes(:, i)
        end forall
        c_lambda12i_c = matmul(c_lambda14i, transpose(c_lambda14i))
        do i_xyz = 1, 3
            relay%re = -relay%re_dr(:, :, i_xyz)
            call form_mbd_matrix(relay, alpha_0, omega)
            relay%re = relay%re-transpose(relay%re)
            sys%work%forces(:, i_xyz) = 1.d0/2*contract_forces(c_lambda12i_c*relay%re)
        end do
    end if
end function get_single_mbd_energy


subroutine form_mbd_matrix(T, alpha_0, omega)
    type(mat3n3n), intent(inout) :: T
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: omega(:)

    real(dp), allocatable :: aw(:)

    aw = repeatn(omega*sqrt(alpha_0), 3)
    call mult_cprod(T, aw, aw)
    call add_diag(T, repeatn(omega**2, 3))
end subroutine form_mbd_matrix


real(dp) function get_reciprocal_mbd_energy(sys, alpha_0, C6, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp

    logical :: &
        is_parallel, do_rpa, mute
    integer :: i_kpt, n_kpts, n_atoms
    real(dp) :: k_point(3), alpha_ts(0:ubound(sys%calc%omega_grid, 1), size(sys%coords, 1))

    n_atoms = size(sys%coords, 1)
    sys%work%k_pts = make_k_grid( &
        make_g_grid(sys%calc, sys%k_grid(1), sys%k_grid(2), sys%k_grid(3)), sys%lattice &
    )
    n_kpts = size(sys%work%k_pts, 1)
    is_parallel = sys%calc%parallel
    do_rpa = sys%do_rpa
    mute = sys%calc%mute

    sys%calc%parallel = .false.

    alpha_ts = alpha_dynamic_ts(sys%calc, alpha_0, C6)
    ene = 0.d0
    if (sys%work%get_eigs) &
        allocate (sys%work%mode_enes_k(n_kpts, 3*n_atoms), source=0.d0)
    if (sys%work%get_modes) &
        allocate (sys%work%modes_k(n_kpts, 3*n_atoms, 3*n_atoms), source=(0.d0, 0.d0))
    if (sys%work%get_rpa_orders) &
        allocate (sys%work%rpa_orders_k(n_kpts, sys%calc%param%rpa_order_max), source=0.d0)
    do i_kpt = 1, n_kpts
        ! MPI code begin
        if (is_parallel) then
            if (sys%calc%my_task /= modulo(i_kpt, sys%calc%n_tasks)) cycle
        end if
        ! MPI code end
        k_point = sys%work%k_pts(i_kpt, :)
        sys%work%i_kpt = i_kpt
        if (do_rpa) then
            ene = ene + get_single_reciprocal_rpa_ene(sys, alpha_ts, k_point, damp)
        else
            ene = ene + get_single_reciprocal_mbd_ene(sys, alpha_0, C6, k_point, damp)
        end if
        sys%calc%mute = .true.
    end do ! k_point loop
    ! MPI code begin
    if (is_parallel) then
        call sync_sum(ene, sys%calc%comm)
        if (sys%work%get_eigs) call sync_sum(sys%work%mode_enes_k, sys%calc%comm)
        if (sys%work%get_modes) call sync_sum(sys%work%modes_k, sys%calc%comm)
        if (sys%work%get_rpa_orders) call sync_sum(sys%work%rpa_orders_k, sys%calc%comm)
    end if
    ! MPI code end
    ene = ene/size(sys%work%k_pts, 1)
    if (sys%work%get_rpa_orders) sys%work%rpa_orders = sys%work%rpa_orders/n_kpts

    sys%calc%parallel = is_parallel
    sys%calc%mute = mute
end function get_reciprocal_mbd_energy


real(dp) function get_single_reciprocal_mbd_ene(sys, alpha_0, C6, k_point, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    real(dp), intent(in) :: k_point(3)
    type(mbd_damping), intent(in) :: damp


    type(mat3n3n) :: relay
    real(dp), allocatable :: eigs(:)
    real(dp), allocatable :: omega(:)
    integer :: i_atom, j_atom, i, j
    integer :: n_negative_eigs
    logical :: is_parallel

    is_parallel = sys%calc%parallel

    allocate (eigs(3*size(sys%coords, 1)))
    omega = omega_eff(C6, alpha_0)
    ! relay = T
    relay = dipole_matrix(sys, damp, k_point)
    do i_atom = 1, size(sys%coords, 1)
        do j_atom = i_atom, size(sys%coords, 1)
            i = 3*(i_atom-1)
            j = 3*(j_atom-1)
            relay%cplx(i+1:i+3, j+1:j+3) = & ! relay = sqrt(a*a)*w*w*T
                omega(i_atom)*omega(j_atom) &
                *sqrt(alpha_0(i_atom)*alpha_0(j_atom))* &
                relay%cplx(i+1:i+3, j+1:j+3)
        end do
    end do
    ! relay = w^2+sqrt(a*a)*w*w*T
    call add_diag(relay, repeatn(omega**2, 3))
    call ts(sys%calc, 22)
    if (.not. is_parallel .or. sys%calc%my_task == 0) then
        if (sys%work%get_modes) then
            call sdiagonalize('V', relay%cplx, eigs, sys%work%exc)
            sys%work%modes_k(sys%work%i_kpt, :, :) = relay%cplx
        else
            call sdiagonalize('N', relay%cplx, eigs, sys%work%exc)
        end if
        if (has_exc(sys)) return
    end if
    ! MPI code begin
    if (is_parallel) then
        call broadcast(relay%cplx, sys%calc%comm)
        call broadcast(eigs, sys%calc%comm)
    end if
    ! MPI code end
    call ts(sys%calc, -22)
    if (sys%work%get_eigs) then
        sys%work%mode_enes = sqrt(eigs)
        where (eigs < 0) sys%work%mode_enes = 0.d0
    end if
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        call printer(sys%calc%io, &
            "CDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues" &
        )
        if (sys%calc%param%zero_negative_eigs) where (eigs < 0) eigs = 0.d0
    end if
    ene = 1.d0/2*sum(sqrt(eigs))-3.d0/2*sum(omega)
end function get_single_reciprocal_mbd_ene


real(dp) function get_single_rpa_energy(sys, alpha, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha(0:, :)
    type(mbd_damping), intent(in) :: damp

    type(mat3n3n) :: relay, AT
    complex(dp), allocatable :: eigs(:)
    integer :: i_atom, i_grid_omega, i
    integer :: n_order, n_negative_eigs
    logical :: is_parallel, mute
    type(mbd_damping) :: damp_alpha

    is_parallel = sys%calc%parallel
    mute = sys%calc%mute

    sys%calc%parallel = .false.

    ene = 0.d0
    damp_alpha = damp
    allocate (eigs(3*size(sys%coords, 1)))
    do i_grid_omega = 0, ubound(sys%calc%omega_grid, 1)
        ! MPI code begin
        if (is_parallel) then
            if (sys%calc%my_task /= modulo(i_grid_omega, sys%calc%n_tasks)) cycle
        end if
        ! MPI code end
        damp_alpha%sigma = get_sigma_selfint(alpha(i_grid_omega, :))
        ! relay = T
        relay = dipole_matrix(sys, damp_alpha)
        do i_atom = 1, size(sys%coords, 1)
            i = 3*(i_atom-1)
            relay%re(i+1:i+3, :i) = &
                alpha(i_grid_omega, i_atom)*transpose(relay%re(:i, i+1:i+3))
        end do
        do i_atom = 1, size(sys%coords, 1)
            i = 3*(i_atom-1)
            relay%re(i+1:i+3, i+1:) = &
                alpha(i_grid_omega, i_atom)*relay%re(i+1:i+3, i+1:)
        end do
        ! relay = alpha*T
        if (sys%work%get_rpa_orders) AT = relay
        ! relay = 1+alpha*T
        call add_diag(relay, 1.d0)
        call ts(sys%calc, 23)
        call diagonalize('N', relay%re, eigs, sys%work%exc)
        call ts(sys%calc, -23)
        if (has_exc(sys)) return
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            call printer(sys%calc%io, "1+AT matrix has " &
                //trim(tostr(n_negative_eigs))//" negative eigenvalues")
        end if
        ene = ene+1.d0/(2*pi)*sum(log(dble(eigs)))*sys%calc%omega_grid_w(i_grid_omega)
        if (sys%work%get_rpa_orders) then
            call ts(sys%calc, 24)
            call diagonalize('N', AT%re, eigs, sys%work%exc)
            call ts(sys%calc, -24)
            if (has_exc(sys)) return
            allocate (sys%work%rpa_orders(sys%calc%param%rpa_order_max))
            do n_order = 2, sys%calc%param%rpa_order_max
                sys%work%rpa_orders(n_order) = sys%work%rpa_orders(n_order) &
                    +(-1.d0/(2*pi)*(-1)**n_order &
                    *sum(dble(eigs)**n_order)/n_order) &
                    *sys%calc%omega_grid_w(i_grid_omega)
            end do
        end if
        sys%calc%mute = .true.
    end do
    if (is_parallel) then
        call sync_sum(ene, sys%calc%comm)
        if (sys%work%get_rpa_orders) then
            call sync_sum(sys%work%rpa_orders, sys%calc%comm)
        end if
    end if
end function get_single_rpa_energy


real(dp) function get_single_reciprocal_rpa_ene(sys, alpha, k_point, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha(0:, :)
    real(dp), intent(in) :: k_point(3)
    type(mbd_damping), intent(in) :: damp

    type(mat3n3n) :: relay, AT
    complex(dp), allocatable :: eigs(:)
    integer :: i_atom, i_grid_omega, i
    integer :: n_order, n_negative_eigs
    logical :: is_parallel, mute
    type(mbd_damping) :: damp_alpha

    is_parallel = sys%calc%parallel
    mute = sys%calc%mute

    sys%calc%parallel = .false.

    ene = 0.d0
    damp_alpha = damp
    allocate (eigs(3*size(sys%coords, 1)))
    do i_grid_omega = 0, ubound(sys%calc%omega_grid, 1)
        ! MPI code begin
        if (is_parallel) then
            if (sys%calc%my_task /= modulo(i_grid_omega, sys%calc%n_tasks)) cycle
        end if
        ! MPI code end
        damp_alpha%sigma = get_sigma_selfint(alpha(i_grid_omega, :))
        ! relay = T
        relay = dipole_matrix(sys, damp_alpha, k_point)
        do i_atom = 1, size(sys%coords, 1)
            i = 3*(i_atom-1)
            relay%cplx(i+1:i+3, :i) = &
                alpha(i_grid_omega, i_atom)*conjg(transpose(relay%cplx(:i, i+1:i+3)))
        end do
        do i_atom = 1, size(sys%coords, 1)
            i = 3*(i_atom-1)
            relay%cplx(i+1:i+3, i+1:) = &
                alpha(i_grid_omega, i_atom)*relay%cplx(i+1:i+3, i+1:)
        end do
        ! relay = alpha*T
        if (sys%work%get_rpa_orders) AT = relay
        do i = 1, 3*size(sys%coords, 1)
            relay%cplx(i, i) = 1.d0+relay%cplx(i, i) ! relay = 1+alpha*T
        end do
        call ts(sys%calc, 25)
        call diagonalize('N', relay%cplx, eigs, sys%work%exc)
        if (has_exc(sys)) return
        call ts(sys%calc, -25)
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            call printer(sys%calc%io, "1+AT matrix has " &
                //trim(tostr(n_negative_eigs))//" negative eigenvalues")
        end if
        ene = ene+1.d0/(2*pi)*dble(sum(log(eigs)))*sys%calc%omega_grid_w(i_grid_omega)
        if (sys%work%get_rpa_orders) then
            call ts(sys%calc, 26)
            call diagonalize('N', AT%cplx, eigs, sys%work%exc)
            if (has_exc(sys)) return
            call ts(sys%calc, -26)
            do n_order = 2, sys%calc%param%rpa_order_max
                sys%work%rpa_orders_k(sys%work%i_kpt, n_order) = &
                    sys%work%rpa_orders_k(sys%work%i_kpt, n_order) &
                    +(-1.d0)/(2*pi)*(-1)**n_order &
                    *dble(sum(eigs**n_order))/n_order &
                    *sys%calc%omega_grid_w(i_grid_omega)
            end do
        end if
        sys%calc%mute = .true.
    end do
    if (is_parallel) then
        call sync_sum(ene, sys%calc%comm)
        if (sys%work%get_rpa_orders) then
            call sync_sum(sys%work%rpa_orders_k(sys%work%i_kpt, :), sys%calc%comm)
        end if
    end if

    sys%calc%parallel = is_parallel
    sys%calc%mute = mute
end function get_single_reciprocal_rpa_ene


! function mbd_nbody( &
!         xyz, &
!         alpha_0, &
!         omega, &
!         version, &
!         R_vdw, beta, a, &
!         calc%my_task, calc%n_tasks) &
!         result(ene_orders)
!     real(dp), intent(in) :: &
!         xyz(:, :), &
!         alpha_0(size(xyz, 1)), &
!         omega(size(xyz, 1)), &
!         R_vdw(size(xyz, 1)), &
!         beta, a
!     character(len=*), intent(in) :: version
!     integer, intent(in), optional :: calc%my_task, calc%n_tasks
!     real(dp) :: ene_orders(20)
!
!     integer :: &
!         multi_index(calc%param%mbd_nbody_max), i_body, j_body, i_tuple, &
!         i_atom_ind, j_atom_ind, i_index
!     real(dp) :: ene
!     logical :: is_parallel
!     
!     is_parallel = .false.
!     if (present(calc%n_tasks)) then
!         if (calc%n_tasks > 0) then
!             is_parallel = .true.
!         end if
!     end if
!     ene_orders(:) = 0.d0
!     do i_body = 2, calc%param%mbd_nbody_max
!         i_tuple = 0
!         multi_index(1:i_body-1) = 1
!         multi_index(i_body:calc%param%mbd_nbody_max) = 0
!         do
!             multi_index(i_body) = multi_index(i_body)+1
!             do i_index = i_body, 2, -1
!                 if (multi_index(i_index) > size(xyz, 1)) then
!                     multi_index(i_index) = 1
!                     multi_index(i_index-1) = multi_index(i_index-1)+1
!                 end if
!             end do
!             if (multi_index(1) > size(xyz, 1)) exit
!             if (any(multi_index(1:i_body-1)-multi_index(2:i_body) >= 0)) cycle
!             i_tuple = i_tuple+1
!             if (is_parallel) then
!                 if (calc%my_task /= modulo(i_tuple, calc%n_tasks)) cycle
!             end if
!             ene = get_mbd_energy( &
!                 xyz(multi_index(1:i_body), :), &
!                 alpha_0(multi_index(1:i_body)), &
!                 omega(multi_index(1:i_body)), &
!                 version, &
!                 R_vdw(multi_index(1:i_body)), &
!                 beta, a)
!             ene_orders(i_body) = ene_orders(i_body) &
!                 +ene+3.d0/2*sum(omega(multi_index(1:i_body)))
!         end do ! i_tuple
!     end do ! i_body
!     if (is_parallel) then
!         call sync_sum(ene_orders, size(ene_orders))
!     end if
!     ene_orders(1) = 3.d0/2*sum(omega)
!     do i_body = 2, min(calc%param%mbd_nbody_max, size(xyz, 1))
!         do j_body = 1, i_body-1
!             ene_orders(i_body) = ene_orders(i_body) &
!                 -nbody_coeffs(j_body, i_body, size(xyz, 1))*ene_orders(j_body)
!         end do
!     end do
!     ene_orders(1) = sum(ene_orders(2:calc%param%mbd_nbody_max))
! end function mbd_nbody


function eval_mbd_nonint_density(calc, pts, xyz, charges, masses, omegas) result(rho)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(in) :: &
        pts(:, :), &
        xyz(:, :), &
        charges(:), &
        masses(:), &
        omegas(:)
    real(dp) :: rho(size(pts, 1))

    integer :: i_pt, i_atom, n_atoms
    real(dp), dimension(size(xyz, 1)) :: pre, kernel, rsq

    pre = charges*(masses*omegas/pi)**(3.d0/2)
    kernel = masses*omegas
    n_atoms = size(xyz, 1)
    rho(:) = 0.d0
    do i_pt = 1, size(pts, 1)
        if (calc%my_task /= modulo(i_pt, calc%n_tasks)) cycle
        forall (i_atom = 1:n_atoms)
            rsq(i_atom) = sum((pts(i_pt, :)-xyz(i_atom, :))**2)
        end forall
        rho(i_pt) = sum(pre*exp(-kernel*rsq))
    end do
    if (calc%parallel) then
        call sync_sum(rho, calc%comm)
    end if
end function


function eval_mbd_int_density(calc, pts, xyz, charges, masses, omegas, modes) result(rho)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(in) :: &
        pts(:, :), &
        xyz(:, :), &
        charges(:), &
        masses(:), &
        omegas(:), &
        modes(:, :)
    real(dp) :: rho(size(pts, 1))

    integer :: i_pt, i_atom, n_atoms, i, i_xyz, j_xyz
    integer :: self(3), other(3*(size(xyz, 1)-1))
    real(dp) :: &
        pre(size(xyz, 1)), &
        factor(size(xyz, 1)), &
        rdiffsq(3, 3), &
        omegas_p(3*size(xyz, 1), 3*size(xyz, 1)), &
        kernel(3, 3, size(xyz, 1)), &
        rdiff(3)

    omegas_p = matmul(matmul(modes, diag(omegas)), transpose(modes))
    n_atoms = size(xyz, 1)
    kernel(:, :, :) = 0.d0
    pre(:) = 0.d0
    do i_atom = 1, n_atoms
        if (calc%my_task /= modulo(i_atom, calc%n_tasks)) cycle
        self(:) = (/ (3*(i_atom-1)+i, i = 1, 3) /)
        other(:) = (/ (i, i = 1, 3*(i_atom-1)),  (i, i = 3*i_atom+1, 3*n_atoms) /)
        kernel(:, :, i_atom) = masses(i_atom) &
            *(omegas_p(self, self) &
                -matmul(matmul(omegas_p(self, other), inverted(omegas_p(other, other))), &
                    omegas_p(other, self)))
        pre(i_atom) = charges(i_atom)*(masses(i_atom)/pi)**(3.d0/2) &
            *sqrt(product(omegas)/product(sdiagonalized(omegas_p(other, other))))
    end do
    if (calc%parallel) then
        call sync_sum(kernel, calc%comm)
        call sync_sum(pre, calc%comm)
    end if
    rho(:) = 0.d0
    do i_pt = 1, size(pts, 1)
        if (calc%my_task /= modulo(i_pt, calc%n_tasks)) cycle
        do i_atom = 1, n_atoms
            rdiff(:) = pts(i_pt, :)-xyz(i_atom, :)
            forall (i_xyz = 1:3, j_xyz = 1:3)
                rdiffsq(i_xyz, j_xyz) = rdiff(i_xyz)*rdiff(j_xyz)
            end forall
            factor(i_atom) = sum(kernel(:, :, i_atom)*rdiffsq(:, :))
        end do
        rho(i_pt) = sum(pre*exp(-factor))
    end do
    if (calc%parallel) then
        call sync_sum(rho, calc%comm)
    end if
end function


function nbody_coeffs(k, m, N) result(a)
    integer, intent(in) :: k, m, N
    integer :: a

    integer :: i

    a = 1
    do i = N-m+1, N-k
        a = a*i
    end do
    do i = 1, m-k
        a = a/i
    end do
end function nbody_coeffs


function contract_polarizability(alpha_3n_3n) result(alpha_n)
    real(dp), intent(in) :: alpha_3n_3n(:, :)
    real(dp) :: alpha_n(size(alpha_3n_3n, 1)/3)

    integer :: i_atom, i_xyz

    alpha_n(:) = 0.d0
    do i_atom = 1, size(alpha_n)
        associate (A => alpha_n(i_atom))
            do i_xyz = 1, 3
                ! this convoluted contraction is necessary because alpha_3n_3n is
                ! calucated as upper triangular
                A = A + sum(alpha_3n_3n(i_xyz:3*i_atom:3, 3*(i_atom-1)+i_xyz))
                A = A + sum(alpha_3n_3n(3*(i_atom-1)+i_xyz, 3*i_atom+i_xyz::3))
            end do
        end associate
    end do
    alpha_n = alpha_n/3
end function contract_polarizability


function contract_forces(relay) result(atomvec)
    real(dp), intent(in) :: relay(:, :)
    real(dp) :: atomvec(size(relay, 1)/3)

    integer :: i_atom

    atomvec(:) = 0.d0
    do i_atom = 1, size(atomvec)
        associate (A => atomvec(i_atom), i => 3*(i_atom-1))
            A = A + sum(relay(:, i+1:i+3))
        end associate
    end do
    atomvec = atomvec
end function contract_forces


function T_bare(rxyz) result(T)
    real(dp), intent(in) :: rxyz(3)
    real(dp) :: T(3, 3)

    integer :: i, j
    real(dp) :: r_sq, r_5

    r_sq = sum(rxyz(:)**2)
    r_5 = sqrt(r_sq)**5
    do i = 1, 3
        T(i, i) = (3.d0*rxyz(i)**2-r_sq)/r_5
        do j = i+1, 3
            T(i, j) = 3.d0*rxyz(i)*rxyz(j)/r_5
            T(j, i) = T(i, j)
        end do
    end do
    T = -T
end function


type(mat33) function T_bare_v2(r, deriv) result(T)
    real(dp), intent(in) :: r(3)
    logical, intent(in) :: deriv

    integer :: a, b, c
    real(dp) :: r_1, r_2, r_5, r_7

    r_2 = sum(r**2)
    r_1 = sqrt(r_2)
    r_5 = r_1**5
    forall (a = 1:3)
        T%val(a, a) = (-3*r(a)**2+r_2)/r_5
        forall (b = a+1:3)
            T%val(a, b) = -3*r(a)*r(b)/r_5
            T%val(b, a) = T%val(a, b)
        end forall
    end forall
    if (deriv) then
        allocate (T%dr(3, 3, 3))
        r_7 = r_1**7
        forall (a = 1:3)
            T%dr(a, a, a) = -3*(3*r(a)/r_5-5*r(a)**3/r_7)
            forall (b = a+1:3)
                T%dr(a, a, b) = -3*(r(b)/r_5-5*r(a)**2*r(b)/r_7)
                T%dr(a, b, a) = T%dr(a, a, b)
                T%dr(b, a, a) = T%dr(a, a, b)
                T%dr(b, b, a) = -3*(r(a)/r_5-5*r(b)**2*r(a)/r_7)
                T%dr(b, a, b) = T%dr(b, b, a)
                T%dr(a, b, b) = T%dr(b, b, a)
                forall (c = b+1:3)
                    T%dr(a, b, c) = 15*r(a)*r(b)*r(c)/r_7
                    T%dr(a, c, b) = T%dr(a, b, c)
                    T%dr(b, a, c) = T%dr(a, b, c)
                    T%dr(b, c, a) = T%dr(a, b, c)
                    T%dr(c, a, b) = T%dr(a, b, c)
                    T%dr(c, b, a) = T%dr(a, b, c)
                end forall
            end forall
        end forall
    end if
end function


real(dp) function B_erfc(r, a) result(B)
    real(dp), intent(in) :: r, a

    B = (erfc(a*r)+(2*a*r/sqrt(pi))*exp(-(a*r)**2))/r**3
end function


real(dp) elemental function C_erfc(r, a) result(C)
    real(dp), intent(in) :: r, a

    C = (3*erfc(a*r)+(2*a*r/sqrt(pi))*(3.d0+2*(a*r)**2)*exp(-(a*r)**2))/r**5
end function


function T_erfc(rxyz, alpha) result(T)
    real(dp), intent(in) :: rxyz(3), alpha
    real(dp) :: T(3, 3)

    integer :: i, j
    real(dp) :: r, B, C

    r = sqrt(sum(rxyz(:)**2))
    B = B_erfc(r, alpha)
    C = C_erfc(r, alpha)
    do i = 1, 3
        do j = i, 3
            T(i, j) = -C*rxyz(i)*rxyz(j)
            if (i /= j) T(j, i) = T(i, j)
        end do
        T(i, i) = T(i, i)+B
    end do
end function


type(scalar) function damping_fermi(r, s_vdw, d, deriv) result(f)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: s_vdw
    real(dp), intent(in) :: d
    logical, intent(in) :: deriv

    real(dp) :: pre, eta, r_1

    r_1 = sqrt(sum(r**2))
    eta = r_1/s_vdw
    f%val = 1.d0/(1+exp(-d*(eta-1)))
    pre = d/(2+2*cosh(d-d*eta))
    if (deriv) then
        f%dr = pre*r/(r_1*s_vdw)
        f%dvdw = -pre*r_1/s_vdw**2
    end if
end function


type(scalar) function damping_sqrtfermi(r, s_vdw, d, deriv) result(f)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: s_vdw
    real(dp), intent(in) :: d
    logical, intent(in) :: deriv

    f = damping_fermi(r, s_vdw, d, deriv)
    f%val = sqrt(f%val)
end function


type(scalar) function op1minus(f)
    type(scalar), intent(in) :: f

    op1minus%val = 1-f%val
    if (allocated(f%dr)) op1minus%dr = -f%dr
    if (allocated(f%dvdw)) op1minus%dvdw = -f%dvdw
end function


type(mat33) function T_damped__(f, T) result(fT)
    type(scalar), intent(in) :: f
    type(mat33), intent(in) :: T

    integer :: c

    fT%val = f%val*T%val
    if (allocated(f%dr) .or. allocated(T%dr)) &
        allocate (fT%dr(3, 3, 3), source=0.d0)
    if (allocated(f%dvdw) .or. allocated(T%dvdw)) &
        allocate (fT%dvdw(3, 3), source=0.d0)
    if (allocated(f%dr)) forall (c = 1:3) fT%dr(:, :, c) = f%dr(c)*T%val
    if (allocated(T%dr)) fT%dr = fT%dr + f%val*T%dr
    if (allocated(f%dvdw)) fT%dvdw = f%dvdw*T%val
    if (allocated(T%dvdw)) fT%dvdw = fT%dvdw + f%val*T%dvdw
    if (allocated(T%dsigma)) fT%dsigma = f%val*T%dsigma
end function


type(mat33) function T_erf_coulomb(r, sigma, deriv) result(T)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: sigma
    logical, intent(in) :: deriv

    real(dp) :: theta, erf_theta, r_5, r_1, zeta
    type(mat33) :: bare
    real(dp) :: tmp33(3, 3), tmp333(3, 3, 3), rr_r5(3, 3)
    integer :: a, c

    bare = T_bare_v2(r, deriv)
    r_1 = sqrt(sum(r**2))
    r_5 = r_1**5
    rr_r5 = (r.cprod.r)/r_5
    zeta = r_1/sigma
    theta = 2*zeta/sqrt(pi)*exp(-zeta**2)
    erf_theta = erf(zeta)-theta
    T%val = erf_theta*bare%val+2*(zeta**2)*theta*rr_r5
    if (deriv) then
        allocate (T%dr(3, 3, 3))
        tmp33 = 2*zeta*theta*(bare%val+(3-2*zeta**2)*rr_r5)
        forall (c = 1:3) T%dr(:, :, c) = tmp33*r(c)/(r_1*sigma)
        tmp333 = bare%dr/3
        forall (a = 1:3, c = 1:3) tmp333(a, a, c) = tmp333(a, a, c) + r(c)/r_5
        T%dr = T%dr + erf_theta*bare%dr-2*(zeta**2)*theta*tmp333
        T%dsigma = -tmp33*r_1/sigma**2
    end if
end function


function T_1mexp_coulomb(rxyz, sigma, a) result(T)
    real(dp), intent(in) :: rxyz(3), sigma, a
    real(dp) :: T(3, 3)

    real(dp) :: r_sigma, zeta_1, zeta_2

    r_sigma = (sqrt(sum(rxyz**2))/sigma)**a
    zeta_1 = 1.d0-exp(-r_sigma)-a*r_sigma*exp(-r_sigma)
    zeta_2 = -r_sigma*a*exp(-r_sigma)*(1+a*(-1+r_sigma))
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function


subroutine get_damping_parameters(xc, ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, &
        mbd_ts_erf_beta, mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta)
    character(len=*), intent(in) :: xc
    real(dp), intent(out) :: &
        ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, mbd_ts_erf_beta, &
        mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta

    ts_d = 20.d0
    ts_s_r = 1.d0
    mbd_scs_a = 2.d0
    mbd_ts_a = 6.d0
    mbd_ts_erf_beta = 1.d0
    mbd_ts_fermi_beta = 1.d0
    mbd_rsscs_a = 6.d0
    mbd_rsscs_beta = 1.d0
    select case (xc)
        case ("pbe")
            ts_s_r = 0.94d0
            mbd_scs_a = 2.56d0
            mbd_ts_erf_beta = 1.07d0
            mbd_ts_fermi_beta = 0.81d0
            mbd_rsscs_beta = 0.83d0
        case ("pbe0")
            ts_s_r = 0.96d0
            mbd_scs_a = 2.53d0
            mbd_ts_erf_beta = 1.08d0
            mbd_ts_fermi_beta = 0.83d0
            mbd_rsscs_beta = 0.85d0
        case ("hse")
            ts_s_r = 0.96d0
            mbd_scs_a = 2.53d0
            mbd_ts_erf_beta = 1.08d0
            mbd_ts_fermi_beta = 0.83d0
            mbd_rsscs_beta = 0.85d0
        case ("blyp")
            ts_s_r = 0.62d0
        case ("b3lyp")
            ts_s_r = 0.84d0
        case ("revpbe")
            ts_s_r = 0.60d0
        case ("am05")
            ts_s_r = 0.84d0
    endselect
end subroutine get_damping_parameters


elemental function terf(r, r0, a)
    real(dp), intent(in) :: r, r0, a
    real(dp) :: terf

    terf = 0.5d0*(erf(a*(r+r0))+erf(a*(r-r0)))
end function


function alpha_dynamic_ts(calc, alpha_0, C6) result(alpha)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    real(dp) :: alpha(0:ubound(calc%omega_grid, 1), size(alpha_0))

    integer :: i_freq
    real(dp), allocatable :: omega(:)

    omega = omega_eff(C6, alpha_0)
    forall (i_freq = 0:ubound(calc%omega_grid, 1))
        alpha(i_freq, :) = alpha_osc(alpha_0, omega, calc%omega_grid(i_freq))
    end forall
end function


elemental function alpha_osc(alpha_0, omega, u) result(alpha)
    real(dp), intent(in) :: alpha_0, omega, u
    real(dp) :: alpha

    alpha = alpha_0/(1+(u/omega)**2)
end function


elemental function combine_C6 (C6_i, C6_j, alpha_0_i, alpha_0_j) result(C6_ij)
    real(dp), intent(in) :: C6_i, C6_j, alpha_0_i, alpha_0_j
    real(dp) :: C6_ij

    C6_ij = 2*C6_i*C6_j/(alpha_0_j/alpha_0_i*C6_i+alpha_0_i/alpha_0_j*C6_j)
end function


elemental function V_to_R(V) result(R)
    real(dp), intent(in) :: V
    real(dp) :: R

    R = (3.d0*V/(4.d0*pi))**(1.d0/3)
end function


elemental function omega_eff(C6, alpha) result(omega)
    real(dp), intent(in) :: C6, alpha
    real(dp) :: omega

    omega = 4.d0/3*C6/alpha**2
end function


elemental function get_sigma_selfint(alpha) result(sigma)
    real(dp), intent(in) :: alpha
    real(dp) :: sigma

    sigma = (sqrt(2.d0/pi)*alpha/3.d0)**(1.d0/3)
end function


function get_C6_from_alpha(calc, alpha) result(C6)
    type(mbd_calc), intent(in) :: calc
    type(vecn), intent(in) :: alpha(0:)
    real(dp) :: C6(size(alpha(0)%val))

    integer :: i_atom, i

    do i_atom = 1, size(alpha(0)%val)
        C6(i_atom) = 3.d0/pi*sum( &
            ([(alpha(i)%val(i_atom), i = 0, ubound(alpha, 1))]**2) * &
            calc%omega_grid_w(:) &
        )
    end do
end function


function get_total_C6_from_alpha(calc, alpha) result(C6)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(in) :: alpha(:, :)
    real(dp) :: C6

    C6 = 3.d0/pi*sum((sum(alpha, 2)**2)*calc%omega_grid_w(:))
end function


function supercell_circum(sys, uc, radius) result(sc)
    type(mbd_system), intent(in) :: sys
    real(dp), intent(in) :: uc(3, 3), radius
    integer :: sc(3)

    real(dp) :: ruc(3, 3), layer_sep(3)
    integer :: i

    ruc = 2*pi*inverted(transpose(uc))
    forall (i = 1:3) layer_sep(i) = sum(uc(i, :)*ruc(i, :)/sqrt(sum(ruc(i, :)**2)))
    sc = ceiling(radius/layer_sep+0.5d0)
    where (sys%vacuum_axis) sc = 0
end function


subroutine shift_cell(ijk, first_cell, last_cell)
    integer, intent(inout) :: ijk(3)
    integer, intent(in) :: first_cell(3), last_cell(3)

    integer :: i_dim, i

    do i_dim = 3, 1, -1
        i = ijk(i_dim)+1
        if (i <= last_cell(i_dim)) then
            ijk(i_dim) = i
            return
        else
            ijk(i_dim) = first_cell(i_dim)
        end if
    end do
end subroutine


function make_g_grid(calc, n1, n2, n3) result(g_grid)
    type(mbd_calc), intent(in) :: calc
    integer, intent(in) :: n1, n2, n3
    real(dp) :: g_grid(n1*n2*n3, 3)

    integer :: g_kpt(3), i_kpt, kpt_range(3)
    real(dp) :: g_kpt_shifted(3)

    g_kpt = (/ 0, 0, -1 /)
    kpt_range = (/ n1, n2, n3 /)
    do i_kpt = 1, n1*n2*n3
        call shift_cell (g_kpt, (/ 0, 0, 0 /), kpt_range-1)
        g_kpt_shifted = dble(g_kpt)+calc%param%k_grid_shift
        where (2*g_kpt_shifted > kpt_range)
            g_kpt_shifted = g_kpt_shifted-dble(kpt_range)
        end where
        g_grid(i_kpt, :) = g_kpt_shifted/kpt_range
    end do
end function make_g_grid


logical function has_exc(sys)
    type(mbd_system), intent(in) :: sys

    has_exc = sys%work%exc%label /= ''
end function


subroutine print_exc(sys)
    type(mbd_system), intent(in) :: sys

    call printer(sys%calc%io, &
        'Error "' // trim(sys%work%exc%label) // '" in "' // &
        trim(sys%work%exc%origin) // ': ' // trim(sys%work%exc%msg))
end subroutine


function make_k_grid(g_grid, uc) result(k_grid)
    real(dp), intent(in) :: g_grid(:, :), uc(3, 3)
    real(dp) :: k_grid(size(g_grid, 1), 3)

    integer :: i_kpt
    real(dp) :: ruc(3, 3)

    ruc = 2*pi*inverted(transpose(uc))
    do i_kpt = 1, size(g_grid, 1)
        k_grid(i_kpt, :) = matmul(g_grid(i_kpt, :), ruc)
    end do
end function make_k_grid


subroutine ts(calc, id, always)
    type(mbd_calc), intent(inout) :: calc
    integer, intent(in) :: id
    logical, intent(in), optional :: always

    if (calc%tm%measure_time .or. present(always)) then
        call system_clock(calc%tm%ts_cnt, calc%tm%ts_rate, calc%tm%ts_cnt_max)
        if (id > 0) then
            calc%tm%timestamps(id) = calc%tm%timestamps(id)-calc%tm%ts_cnt
        else
            calc%tm%ts_aid = abs(id)
            calc%tm%timestamps(calc%tm%ts_aid) = calc%tm%timestamps(calc%tm%ts_aid)+calc%tm%ts_cnt
            calc%tm%ts_counts(calc%tm%ts_aid) = calc%tm%ts_counts(calc%tm%ts_aid)+1
        end if
    end if
end subroutine ts


function clock_rate() result(rate)
    integer :: cnt, rate, cnt_max

    call system_clock(cnt, rate, cnt_max) 
end function clock_rate


!!! tests !!!

subroutine run_tests()
    use mbd_common, only: diff3, tostr, diff5

    integer :: n_failed, n_all
    type(mbd_calc), target :: calc

    call init_grid(calc)
    n_failed = 0
    n_all = 0
    call exec_test('T_bare derivative')
    call exec_test('T_GG derivative explicit')
    call exec_test('T_GG derivative implicit')
    call exec_test('MBD derivative explicit')
    write (6, *) &
        trim(tostr(n_failed)) // '/' // trim(tostr(n_all)) // ' tests failed'
    if (n_failed /= 0) stop 1

    contains

    subroutine exec_test(test_name)
        character(len=*), intent(in) :: test_name

        integer :: n_failed_in

        write (6, '(A,A,A)', advance='no') 'Executing test "', test_name, '"... '
        n_failed_in = n_failed
        select case (test_name)
        case ('T_bare derivative'); call test_T_bare_deriv()
        case ('T_GG derivative explicit'); call test_T_GG_deriv_expl()
        case ('T_GG derivative implicit'); call test_T_GG_deriv_impl()
        case ('MBD derivative explicit'); call test_mbd_deriv_expl()
        end select
        n_all = n_all + 1
        if (n_failed == n_failed_in) write (6, *) 'OK'
    end subroutine

    subroutine failed()
        n_failed = n_failed + 1
        write (6, *) 'FAILED!'
    end subroutine

    subroutine test_T_bare_deriv()
        real(dp) :: r(3), r_diff(3)
        type(mat33) :: T
        real(dp) :: diff(3, 3)
        real(dp) :: T_diff_anl(3, 3, 3)
        real(dp) :: T_diff_num(3, 3, -2:2)
        integer :: a, b, c, i_step
        real(dp) :: delta

        delta = 1d-3
        r = [1.12d0, -2.12d0, 0.12d0]
        T = T_bare_v2(r, deriv=.true.)
        T_diff_anl = T%dr(:, :, :)
        do c = 1, 3
            do i_step = -2, 2
                if (i_step == 0) continue
                r_diff = r
                r_diff(c) = r_diff(c)+i_step*delta
                T = T_bare_v2(r_diff, deriv=.false.)
                T_diff_num(:, :, i_step) = T%val
            end do
            forall (a = 1:3, b = 1:3)
                T_diff_num(a, b, 0) = diff5(T_diff_num(a, b, :), delta)
            end forall
            diff = T_diff_num(:, :, 0)-T_diff_anl(:, :, c)
            if (any(abs(diff) > 1d-12)) then
                call failed()
                call print_matrix('delta dT(:, :, ' // trim(tostr(c)) // ')', diff)
            end if
        end do
    end subroutine test_T_bare_deriv

    subroutine test_T_GG_deriv_expl()
        real(dp) :: r(3), r_diff(3)
        type(mat33) :: T
        real(dp) :: diff(3, 3)
        real(dp) :: T_diff_anl(3, 3, 3)
        real(dp) :: T_diff_num(3, 3, -2:2)
        integer :: a, b, c, i_step
        real(dp) :: delta
        real(dp) :: sigma

        delta = 1d-3
        r = [1.02d0, -2.22d0, 0.15d0]
        sigma = 1.2d0
        T = T_erf_coulomb(r, sigma, deriv=.true.)
        T_diff_anl = T%dr
        do c = 1, 3
            do i_step = -2, 2
                if (i_step == 0) continue
                r_diff = r
                r_diff(c) = r_diff(c)+i_step*delta
                T = T_erf_coulomb(r_diff, sigma, deriv=.false.)
                T_diff_num(:, :, i_step) = T%val
            end do
            forall (a = 1:3, b = 1:3)
                T_diff_num(a, b, 0) = diff5(T_diff_num(a, b, :), delta)
            end forall
            diff = T_diff_num(:, :, 0)-T_diff_anl(:, :, c)
            if (any(abs(diff) > 1d-12)) then
                call failed()
                call print_matrix('delta dTGG_{ab,' // trim(tostr(c)) // '}', diff)
            end if
        end do
    end subroutine test_T_GG_deriv_expl

    subroutine test_T_GG_deriv_impl()
        real(dp) :: r(3)
        type(mat33) :: T
        real(dp) :: diff(3, 3)
        real(dp) :: T_diff_anl(3, 3)
        real(dp) :: T_diff_num(3, 3, -2:2)
        integer :: a, b, i_step
        real(dp) :: delta
        real(dp) :: sigma, dsigma_dr, sigma_diff

        delta = 1d-3
        r = [1.02d0, -2.22d0, 0.15d0]
        sigma = 1.2d0
        dsigma_dr = -0.3d0
        T = T_erf_coulomb(r, sigma, deriv=.true.)
        T_diff_anl = T%dsigma(:, :)*dsigma_dr
        do i_step = -2, 2
            if (i_step == 0) continue
            sigma_diff = sigma+i_step*delta*dsigma_dr
            T = T_erf_coulomb(r, sigma_diff, deriv=.false.)
            T_diff_num(:, :, i_step) = T%val
        end do
        forall (a = 1:3, b = 1:3)
            T_diff_num(a, b, 0) = diff5(T_diff_num(a, b, :), delta)
        end forall
        diff = T_diff_num(:, :, 0)-T_diff_anl
        if (any(abs(diff) > 1d-12)) then
            call failed()
            call print_matrix('delta dTGG', diff)
        end if
    end subroutine test_T_GG_deriv_impl

    subroutine test_mbd_deriv_expl()
        real(dp) :: delta
        type(mbd_system) :: sys
        type(mbd_damping) :: damp
        real(dp), allocatable :: coords(:, :)
        real(dp), allocatable :: forces(:, :)
        real(dp), allocatable :: diff(:, :)
        real(dp), allocatable :: alpha_0(:)
        real(dp), allocatable :: omega(:)
        real(dp) :: ene(-2:2)
        integer :: i_atom, n_atoms, i_xyz, i_step

        delta = 1d-3
        n_atoms = 3
        allocate (coords(n_atoms, 3), source=0.d0)
        allocate (forces(n_atoms, 3))
        coords(2, 1) = 4.d0*ang
        coords(3, 2) = 4.d0*ang
        sys%calc => calc
        sys%coords = coords
        sys%do_force = .true.
        damp%version = 'fermi,dip'
        damp%r_vdw = [3.55d0, 3.55d0, 3.55d0]
        damp%beta = 0.83
        alpha_0 = [11.d0, 11.d0, 11.d0]
        omega = [0.7d0, 0.7d0, 0.7d0]
        ene(0) = get_single_mbd_energy(sys, alpha_0, omega, damp)
        sys%do_force = .false.
        do i_atom = 1, n_atoms
            do i_xyz = 1, 3
                do i_step = -2, 2
                    if (i_step == 0) continue
                    sys%coords = coords
                    sys%coords(i_atom, i_xyz) = sys%coords(i_atom, i_xyz)+i_step*delta
                    ene(i_step) = get_single_mbd_energy(sys, alpha_0, omega, damp)
                end do
                forces(i_atom, i_xyz) = diff5(ene, delta)
            end do
        end do
        diff = forces-sys%work%forces
        if (any(abs(diff) > 1d-12)) then
            call failed()
            call print_matrix('delta forces', diff)
        end if
    end subroutine test_mbd_deriv_expl
end subroutine run_tests

end module mbd
