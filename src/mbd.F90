! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LEGENDRE_PREC
#define LEGENDRE_PREC 8
#endif
module mbd

use mbd_common, only: tostr, nan, print_matrix, dp, pi, exception
use mbd_linalg, only: inv, invh, inverse, eig, eigh, eigvals, eigvalsh
use mbd_types, only: mat3n3n, mat33, scalar, vecn, operator(.cprod.)
use mbd_parallel, only: mbd_blacs_grid, mbd_blacs
use mbd_defaults

implicit none

#ifndef MODULE_UNIT_TESTS
private
public :: mbd_param, mbd_calc, mbd_damping, mbd_result, mbd_system, &
    init_grid, get_mbd_energy, dipole_matrix, mbd_rsscs_energy, mbd_scs_energy, &
    get_sigma_selfint, scale_TS
public :: get_ts_energy, init_eqi_grid, get_damping_parameters, clock_rate
#endif

interface operator(.prod.)
    module procedure T_damped__
end interface

real(dp), parameter :: ang = 1.8897259886d0
integer, parameter :: n_timestamps = 100

type :: mbd_param
    real(dp) :: ts_energy_accuracy = TS_ENERGY_ACCURACY
    real(dp) :: ts_cutoff_radius = 50.d0*ang
    real(dp) :: dipole_low_dim_cutoff = 100.d0*ang
    real(dp) :: dipole_cutoff = 400.d0*ang  ! used only when Ewald is off
    real(dp) :: ewald_real_cutoff_scaling = 1.d0
    real(dp) :: ewald_rec_cutoff_scaling = 1.d0
    real(dp) :: k_grid_shift = K_GRID_SHIFT
    logical :: ewald_on = .true.
    logical :: zero_negative_eigs = .false.
    integer :: rpa_order_max = 10
    integer :: n_frequency_grid = N_FREQUENCY_GRID
end type

type :: mbd_timing
    logical :: measure_time = .true.
    integer :: timestamps(n_timestamps), ts_counts(n_timestamps)
    integer :: ts_cnt, ts_rate, ts_cnt_max, ts_aid
end type mbd_timing

type :: mbd_info
    character(len=120) :: ewald_alpha = '', ewald_rsum = '', ts_conv = '', &
        ewald_cutoff = '', ewald_recsum = '', freq_n = '', freq_error = '', &
        neg_eig = ''
end type mbd_info

type :: mbd_calc
    type(mbd_param) :: param
    type(mbd_timing) :: tm
    real(dp), allocatable :: omega_grid(:)
    real(dp), allocatable :: omega_grid_w(:)
    integer :: comm = -1
    type(exception) :: exc
    type(mbd_info) :: info
end type mbd_calc

type :: mbd_damping
    character(len=20) :: version
    real(dp) :: beta = 0.d0
    real(dp) :: a = MBD_DAMPING_A
    real(dp) :: ts_d = TS_DAMPING_D
    real(dp) :: ts_sr = 0.d0
    real(dp) :: mayer_scaling = 1.d0
    type(vecn) :: r_vdw
    type(vecn) :: sigma
    real(dp), allocatable :: damping_custom(:, :)
    real(dp), allocatable :: potential_custom(:, :, :, :)
end type mbd_damping

type :: mbd_result
    real(dp) :: energy
    real(dp), allocatable :: k_pts(:, :)
    real(dp), allocatable :: mode_enes(:)
    real(dp), allocatable :: modes(:, :)
    real(dp), allocatable :: rpa_orders(:)
    real(dp), allocatable :: mode_enes_k(:, :)
    complex(dp), allocatable :: modes_k(:, :, :)
    complex(dp), allocatable :: modes_k_single(:, :)
    real(dp), allocatable :: rpa_orders_k(:, :)
    real(dp), allocatable :: gradients(:, :)
end type

type :: mbd_system
    type(mbd_calc), pointer :: calc
    real(dp), allocatable :: coords(:, :)  ! 3 by n_atoms
    logical :: periodic = .false.
    logical :: vacuum_axis(3) = (/ .false., .false., .false. /)
    real(dp) :: lattice(3, 3)  ! vectors in columns
    integer :: k_grid(3)
    integer :: supercell(3)
    logical :: do_rpa = .false.
    logical :: do_reciprocal = .true.
    logical :: do_gradients = .false.
    logical :: get_eigs = .false.
    logical :: get_modes = .false.
    logical :: get_rpa_orders = .false.
    type(mbd_blacs_grid) :: blacs_grid
    contains
    procedure :: siz => system_siz
    procedure :: has_exc => mbd_system_has_exc
end type mbd_system

contains


type(mbd_result) function mbd_rsscs_energy(sys, alpha_0, C6, damp)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha_0
    type(vecn), intent(in) :: C6
    type(mbd_damping), intent(in) :: damp

    type(vecn), allocatable :: alpha_dyn(:), alpha_dyn_rsscs(:)
    type(vecn) :: C6_rsscs
    type(mbd_damping) :: damp_rsscs, damp_mbd
    integer :: n_freq, i_freq

    n_freq = ubound(sys%calc%omega_grid, 1)
    allocate (alpha_dyn(0:n_freq))
    allocate (alpha_dyn_rsscs(0:n_freq))
    alpha_dyn = alpha_dynamic_ts(sys%calc, alpha_0, C6)
    damp_rsscs = damp
    damp_rsscs%version = 'fermi,dip,gg'
    do i_freq = 0, n_freq
        alpha_dyn_rsscs(i_freq) = run_scs(sys, alpha_dyn(i_freq), damp_rsscs)
    end do
    if (sys%has_exc()) return
    C6_rsscs = get_C6_from_alpha(sys%calc, alpha_dyn_rsscs)
    damp_mbd%version = 'fermi,dip'
    damp_mbd%r_vdw = scale_TS(damp%R_vdw, alpha_dyn_rsscs(0), alpha_dyn(0), 1.d0/3)
    damp_mbd%beta = damp%beta
    mbd_rsscs_energy = get_mbd_energy(sys, alpha_dyn_rsscs(0), C6_rsscs, damp_mbd)
end function mbd_rsscs_energy


type(mbd_result) function mbd_scs_energy(sys, alpha_0, C6, damp)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha_0
    type(vecn), intent(in) :: C6
    type(mbd_damping), intent(in) :: damp

    type(vecn), allocatable :: alpha_dyn(:), alpha_dyn_scs(:)
    type(vecn) :: C6_scs
    type(mbd_damping) :: damp_scs, damp_mbd
    integer :: n_freq, i_freq

    n_freq = ubound(sys%calc%omega_grid, 1)
    allocate (alpha_dyn(0:n_freq))
    allocate (alpha_dyn_scs(0:n_freq))
    alpha_dyn = alpha_dynamic_ts(sys%calc, alpha_0, C6)
    damp_scs = damp
    damp_scs%version = 'dip,gg'
    do i_freq = 0, n_freq
        alpha_dyn_scs(i_freq) = run_scs(sys, alpha_dyn(i_freq), damp_scs)
    end do
    if (sys%has_exc()) return
    C6_scs = get_C6_from_alpha(sys%calc, alpha_dyn_scs)
    damp_mbd%r_vdw = scale_TS(damp%R_vdw, alpha_dyn_scs(0), alpha_dyn(0), 1.d0/3)
    damp_mbd%version = 'dip,1mexp'
    damp_mbd%beta = 1.d0
    damp_mbd%a = damp%a
    mbd_scs_energy = get_mbd_energy(sys, alpha_dyn_scs(0), C6_scs, damp_mbd)
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
    logical :: is_crystal

    is_crystal = sys%periodic
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
            if (is_crystal) then
                R_cell = matmul(sys%lattice, idx_cell)
            else
                R_cell = (/ 0.d0, 0.d0, 0.d0 /)
            end if
            do i_atom = 1, sys%siz()
                do j_atom = 1, i_atom
                    if (i_cell == 1) then
                        if (i_atom == j_atom) cycle
                    end if
                    r = sys%coords(:, i_atom)-sys%coords(:, j_atom)-R_cell
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
                    if (allocated(damp%r_vdw%val)) then
                        R_vdw_ij = damp%r_vdw%val(i_atom)+damp%r_vdw%val(j_atom)
                    end if
                    select case (damp%version)
                        case ("fermi")
                            f_damp = damping_fermi( &
                                r, damp%ts_sr*R_vdw_ij, damp%ts_d, .false. &
                            )
                        case ("fermi2")
                            f_damp = damping_fermi( &
                                r, damp%ts_sr*R_vdw_ij, damp%ts_d, .false. &
                            )
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
        ene = ene+ene_shell
        if (.not. is_crystal) exit
        if (i_shell > 1 .and. &
                abs(ene_shell) < sys%calc%param%ts_energy_accuracy) then
            sys%calc%info%ts_conv = "Periodic TS converged in " // &
                trim(tostr(i_shell)) // " shells, " // &
                trim(tostr(i_shell*shell_thickness/ang)) // " angstroms"
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
    integer :: i_atom, j_atom, i_cell, idx_cell(3), range_cell(3), i, j, &
        n_atoms, my_i_atom, my_j_atom
    logical :: do_ewald

    do_ewald = .false.
    n_atoms = sys%siz()
    call dipmat%init(n_atoms, sys%blacs_grid)
    if (present(k_point)) then
        allocate (dipmat%cplx(3*n_atoms, 3*n_atoms), source=(0.d0, 0.d0))
    else
        allocate (dipmat%re(3*n_atoms, 3*n_atoms), source=0.d0)
        if (sys%do_gradients) then
            allocate (dipmat%re_dr(3*n_atoms, 3*n_atoms, 3), source=0.d0)
            if (allocated(damp%r_vdw%dr)) then
                allocate (dipmat%re_dvdw(3*n_atoms, 3*n_atoms), source=0.d0)
            end if
            if (allocated(damp%sigma%dr)) then
                allocate (dipmat%re_dsigma(3*n_atoms, 3*n_atoms), source=0.d0)
            end if
        end if
    end if
    ! MPI code end
    if (sys%periodic) then
        if (any(sys%vacuum_axis)) then
            real_space_cutoff = sys%calc%param%dipole_low_dim_cutoff
        else if (sys%calc%param%ewald_on) then
            do_ewald = .true.
            volume = max(abs(dble(product(eigvals(sys%lattice)))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1.d0/3)
            real_space_cutoff = &
                6.d0/ewald_alpha*sys%calc%param%ewald_real_cutoff_scaling
            sys%calc%info%ewald_alpha = &
                'Ewald: using alpha = ' // trim(tostr(ewald_alpha)) // &
                ', real cutoff = ' // trim(tostr(real_space_cutoff))
        else
            real_space_cutoff = sys%calc%param%dipole_cutoff
        end if
        range_cell = supercell_circum(sys, sys%lattice, real_space_cutoff)
    else
        range_cell(:) = 0
    end if
    if (sys%periodic) then
        sys%calc%info%ewald_rsum = &
            'Ewald: summing real part in cell vector range of ' // &
            trim(tostr(1+2*range_cell(1))) // 'x' // &
            trim(tostr(1+2*range_cell(2))) // 'x' // &
            trim(tostr(1+2*range_cell(3)))
    end if
    call ts(sys%calc, 11)
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, product(1+2*range_cell)
        call shift_cell(idx_cell, -range_cell, range_cell)
        if (sys%periodic) then
            R_cell = matmul(sys%lattice, idx_cell)
        else
            R_cell(:) = 0.d0
        end if
        do my_i_atom = 1, size(dipmat%blacs%i_atom)
            i_atom = dipmat%blacs%i_atom(my_i_atom)
            do my_j_atom = 1, size(dipmat%blacs%j_atom)
                j_atom = dipmat%blacs%j_atom(my_j_atom)
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                r = sys%coords(:, i_atom)-sys%coords(:, j_atom)-R_cell
                r_norm = sqrt(sum(r**2))
                if (sys%periodic .and. r_norm > real_space_cutoff) cycle
                if (allocated(damp%R_vdw%val)) then
                    R_vdw_ij = sum(damp%R_vdw%val([i_atom, j_atom]))
                end if
                if (allocated(damp%sigma%val)) then
                    sigma_ij = damp%mayer_scaling * &
                        sqrt(sum(damp%sigma%val([i_atom, j_atom])**2))
                end if
                select case (damp%version)
                    case ("bare")
                        Tpp = T_bare_v2(r, sys%do_gradients)
                    case ("dip,1mexp")
                        Tpp%val = T_1mexp_coulomb(r, damp%beta*R_vdw_ij, damp%a)
                    case ("fermi,dip")
                        Tpp = damping_fermi( &
                            r, damp%beta*R_vdw_ij, damp%a, sys%do_gradients &
                        ).prod.T_bare_v2(r, sys%do_gradients)
                    case ("sqrtfermi,dip")
                        Tpp = damping_sqrtfermi( &
                            r, damp%beta*R_vdw_ij, damp%a, sys%do_gradients &
                        ).prod.T_bare_v2(r, sys%do_gradients)
                    case ("custom,dip")
                        Tpp%val = damp%damping_custom(i_atom, j_atom)*T_bare(r)
                    case ("dip,custom")
                        Tpp%val = damp%potential_custom(i_atom, j_atom, :, :)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma_ij, sys%do_gradients)
                    case ("fermi,dip,gg")
                        Tpp = op1minus(damping_fermi( &
                            r, damp%beta*R_vdw_ij, damp%a, sys%do_gradients &
                        )).prod.T_erf_coulomb(r, sigma_ij, sys%do_gradients)
                        do_ewald = .false.
                    case ("sqrtfermi,dip,gg")
                        Tpp = op1minus(damping_sqrtfermi( &
                            r, damp%beta*R_vdw_ij, damp%a, sys%do_gradients &
                        )).prod.T_erf_coulomb(r, sigma_ij, sys%do_gradients)
                        do_ewald = .false.
                    case ("custom,dip,gg")
                        f_ij = 1.d0-damp%damping_custom(i_atom, j_atom)
                        Tpp = T_erf_coulomb(r, sigma_ij, sys%do_gradients)
                        Tpp%val = f_ij*Tpp%val
                        do_ewald = .false.
                end select
                if (allocated(Tpp%dvdw)) then
                    Tpp%dvdw = damp%beta*Tpp%dvdw
                end if
                if (do_ewald) then
                    Tpp%val = Tpp%val+T_erfc(r, ewald_alpha)-T_bare(r)
                end if
                if (present(k_point)) then
                    Tpp_c = Tpp%val*exp(-cmplx(0.d0, 1.d0, 8)*( &
                        dot_product(k_point, r)))
                end if
                i = 3*(my_i_atom-1)
                j = 3*(my_j_atom-1)
                if (present(k_point)) then
                    associate (T => dipmat%cplx(i+1:i+3, j+1:j+3))
                        T = T + Tpp_c
                    end associate
                else
                    associate (T => dipmat%re(i+1:i+3, j+1:j+3))
                        T = T + Tpp%val
                    end associate
                    if (allocated(dipmat%re_dr)) then
                        associate (T => dipmat%re_dr(i+1:i+3, j+1:j+3, :))
                            T = T + Tpp%dr
                        end associate
                    end if
                    if (allocated(dipmat%re_dvdw)) then
                        associate (dTdRvdw => dipmat%re_dvdw(i+1:i+3, j+1:j+3))
                            dTdRvdw = dTdRvdw + Tpp%dvdw
                        end associate
                    end if
                    if (allocated(dipmat%re_dsigma)) then
                        associate (dTdsigma => dipmat%re_dsigma(i+1:i+3, j+1:j+3))
                            dTdsigma = dTdsigma + Tpp%dsigma
                        end associate
                    end if
                end if
            end do ! j_atom
        end do ! i_atom
    end do ! i_cell
    call ts(sys%calc, -11)
    if (do_ewald) then
        call add_ewald_dipole_parts(sys, ewald_alpha, dipmat, k_point)
    end if
end function dipole_matrix


subroutine add_ewald_dipole_parts(sys, alpha, dipmat, k_point)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha
    real(dp), intent(in), optional :: k_point(3)
    type(mat3n3n), intent(inout) :: dipmat

    logical :: do_surface
    real(dp) :: rec_unit_cell(3, 3), volume, G_vector(3), r(3), k_total(3), &
        k_sq, rec_space_cutoff, Tpp(3, 3), k_prefactor(3, 3), elem
    complex(dp) :: Tpp_c(3, 3)
    integer :: &
        i_atom, j_atom, i, j, i_xyz, j_xyz, idx_G_vector(3), i_G_vector, &
        range_G_vector(3), my_i_atom, my_j_atom

    rec_unit_cell = 2*pi*inverse(transpose(sys%lattice))
    volume = abs(dble(product(eigvals(sys%lattice))))
    rec_space_cutoff = 10.d0*alpha*sys%calc%param%ewald_rec_cutoff_scaling
    range_G_vector = supercell_circum(sys, rec_unit_cell, rec_space_cutoff)
    sys%calc%info%ewald_cutoff = 'Ewald: using reciprocal cutoff = ' // &
        trim(tostr(rec_space_cutoff))
    sys%calc%info%ewald_recsum = &
        'Ewald: summing reciprocal part in G vector range of ' // &
        trim(tostr(1+2*range_G_vector(1))) // 'x' // &
        trim(tostr(1+2*range_G_vector(2))) // 'x' // &
        trim(tostr(1+2*range_G_vector(3)))
    call ts(sys%calc, 12)
    idx_G_vector = (/ 0, 0, -1 /)
    do i_G_vector = 1, product(1+2*range_G_vector)
        call shift_cell(idx_G_vector, -range_G_vector, range_G_vector)
        if (i_G_vector == 1) cycle
        G_vector = matmul(rec_unit_cell, idx_G_vector)
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
        do my_i_atom = 1, size(dipmat%blacs%i_atom)
            i_atom = dipmat%blacs%i_atom(my_i_atom)
            do my_j_atom = 1, size(dipmat%blacs%j_atom)
                j_atom = dipmat%blacs%j_atom(my_j_atom)
                r = sys%coords(:, i_atom)-sys%coords(:, j_atom)
                if (present(k_point)) then
                    Tpp_c = k_prefactor*exp(cmplx(0.d0, 1.d0, 8) &
                        *dot_product(G_vector, r))
                else
                    Tpp = k_prefactor*cos(dot_product(G_vector, r))
                end if
                i = 3*(my_i_atom-1)
                j = 3*(my_j_atom-1)
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
        end do ! i_atom
    end do ! i_G_vector
    call dipmat%add_diag_scalar(-4*alpha**3/(3*sqrt(pi))) ! self energy
    do_surface = .true.
    if (present(k_point)) then
        k_sq = sum(k_point**2)
        if (sqrt(k_sq) > 1.d-15) then
            do_surface = .false.
            do my_i_atom = 1, size(dipmat%blacs%i_atom)
            do my_j_atom = 1, size(dipmat%blacs%j_atom)
                do i_xyz = 1, 3
                do j_xyz = 1, 3
                    i = 3*(my_i_atom-1)+i_xyz
                    j = 3*(my_j_atom-1)+j_xyz
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
        do my_i_atom = 1, size(dipmat%blacs%i_atom)
        do my_j_atom = 1, size(dipmat%blacs%j_atom)
            do i_xyz = 1, 3
                i = 3*(my_i_atom-1)+i_xyz
                j = 3*(my_j_atom-1)+i_xyz
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
    calc%info%freq_n = &
        "Initialized a radial integration grid of " // trim(tostr(n)) // &
        " points."
    calc%info%freq_error = &
        "Relative quadrature error in C6 of carbon atom: " // &
        trim(tostr(test_frequency_grid(calc)))
end subroutine


real(dp) function test_frequency_grid(calc) result(error)
    type(mbd_calc), intent(in) :: calc

    type(vecn) :: alpha(0:ubound(calc%omega_grid, 1)), C6

    alpha = alpha_dynamic_ts(calc, vecn([21.d0]), vecn([99.5d0]))
    C6 = get_C6_from_alpha(calc, alpha)
    error = abs(C6%val(1)/99.5d0-1.d0)
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
        Pk(0:k) = ((2*k-1) * &
            (/ 0.0_q, Pk1(0:k-1) /)-(k-1)*(/ Pk2(0:k-2), 0._q, 0._q /))/k
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


type(vecn) function run_scs(sys, alpha, damp) result(alpha_scs)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha
    type(mbd_damping), intent(in) :: damp

    type(mat3n3n) :: alpha_full, dQ_add, dQ
    integer :: n_atoms, i_xyz, i_atom, my_i_atom, my_j_atom
    type(mbd_damping) :: damp_local

    n_atoms = sys%siz()
    damp_local = damp
    damp_local%sigma = get_sigma_selfint(alpha)
    alpha_full = dipole_matrix(sys, damp_local)
    call alpha_full%add_diag(1.d0/alpha%val)
    call ts(sys%calc, 32)
    call invh(alpha_full%re, sys%calc%exc)
    if (sys%has_exc()) return
    call ts(sys%calc, -32)
    alpha_scs%val = contract_polarizability(alpha_full)
    if (.not. sys%do_gradients) return
    dQ = alpha_full
    dQ_add%blacs = alpha_full%blacs
    allocate (alpha_scs%dr(n_atoms, n_atoms, 3))
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            dQ%re(:, :) = 0.d0
            do my_i_atom = 1, size(dQ%blacs%i_atom)
                if (i_atom == dQ%blacs%i_atom(my_i_atom)) then
                    associate (i => (my_i_atom-1)*3)
                        dQ%re(i+1:i+3, :) = dQ%re_dr(i+1:i+3, :, i_xyz)
                    end associate
                end if
            end do
            do my_j_atom = 1, size(dQ%blacs%j_atom)
                if (i_atom == dQ%blacs%j_atom(my_j_atom)) then
                    associate (j => (my_j_atom-1)*3)
                        dQ%re(:, j+1:j+3) = -dQ%re_dr(:, j+1:j+3, i_xyz)
                    end associate
                end if
            end do
            if (allocated(alpha%dr)) then
                call dQ%add_diag(-alpha%dr(:, i_atom, i_xyz)/alpha%val(:)**2)
            end if
            if (allocated(damp%sigma%dr)) then
                dQ_add%re = dQ%re_dsigma
                call dQ_add%mult_dsigma( &
                    damp%sigma%val, damp%sigma%dr(:, i_atom, i_xyz) &
                )
                call dQ%add(dQ_add)
            end if
            if (allocated(damp%r_vdw%dr)) then
                dQ_add%re = dQ%re_dvdw
                call dQ_add%mult_cross_add(damp%r_vdw%dr(:, i_atom, i_xyz))
                call dQ%add(dQ_add)
            end if
            dQ%re = -matmul(alpha_full%re, matmul(dQ%re, alpha_full%re))
            alpha_scs%dr(:, i_atom, i_xyz) = contract_polarizability(dQ)
        end do
    end do
end function run_scs


type(mbd_result) function get_mbd_energy(sys, alpha_0, C6, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha_0
    type(vecn), intent(in) :: C6
    type(mbd_damping), intent(in) :: damp

    logical :: do_rpa, is_reciprocal, is_crystal
    type(vecn), allocatable :: alpha(:)

    is_crystal = sys%periodic
    do_rpa = sys%do_rpa
    is_reciprocal = sys%do_reciprocal
    if (.not. is_crystal) then
        if (.not. do_rpa) then
            ene = get_single_mbd_energy(sys, alpha_0, C6, damp)
        else
            allocate (alpha(0:ubound(sys%calc%omega_grid, 1)))
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


type(mbd_result) function get_supercell_mbd_energy(sys, alpha_0, C6, damp) &
        result(ene)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha_0
    type(vecn), intent(in) :: C6
    type(mbd_damping), intent(in) :: damp

    logical :: do_rpa
    real(dp) :: R_cell(3)
    integer :: idx_cell(3), n_cells, i_atom, i, i_cell, n_atoms

    type(vecn) :: alpha_0_super, C6_super
    type(vecn), allocatable :: alpha_ts_super(:)
    type(mbd_system) :: sys_super
    type(mbd_damping) :: damp_super
    type(mbd_result) :: ene_super

    do_rpa = sys%do_rpa

    sys_super%calc = sys%calc
    n_cells = product(sys%supercell)
    n_atoms = sys%siz()
    do i = 1, 3
        sys_super%lattice(:, i) = sys%lattice(:, i)*sys%supercell(i)
    end do
    allocate (sys_super%coords(3, n_cells*n_atoms))
    allocate (alpha_0_super%val(n_cells*n_atoms))
    allocate (C6_super%val(n_cells*n_atoms))
    allocate (alpha_ts_super(0:ubound(sys%calc%omega_grid, 1)))
    if (allocated(damp%r_vdw%val)) then
        allocate (damp_super%r_vdw%val(n_cells*n_atoms))
    end if
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, n_cells
        call shift_cell(idx_cell, (/ 0, 0, 0 /), sys%supercell-1)
        R_cell = matmul(sys%lattice, idx_cell)
        do i_atom = 1, n_atoms
            i = (i_cell-1)*n_atoms+i_atom
            sys_super%coords(:, i) = sys%coords(:, i_atom)+R_cell
            alpha_0_super%val(i) = alpha_0%val(i_atom)
            C6_super%val(i) = C6%val(i_atom)
            if (allocated(damp%R_vdw%val)) then
                damp_super%R_vdw%val(i) = damp%R_vdw%val(i_atom)
            end if
        end do
    end do
    if (do_rpa) then
        alpha_ts_super = alpha_dynamic_ts(sys%calc, alpha_0_super, C6_super)
        ene_super = get_single_rpa_energy(sys_super, alpha_ts_super, damp_super)
    else
        ene_super = get_single_mbd_energy( &
            sys_super, alpha_0_super, C6_super, damp_super &
        )
    end if
    ene%energy = ene_super%energy/n_cells
    if (sys%get_rpa_orders) then
        ene%rpa_orders = ene_super%rpa_orders/n_cells
    end if
end function get_supercell_mbd_energy


type(mbd_result) function get_single_mbd_energy( &
        sys, alpha_0, C6, damp, k_point) result(ene)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha_0
    type(vecn), intent(in) :: C6
    type(mbd_damping), intent(in) :: damp
    real(dp), intent(in), optional :: k_point(3)

    type(mat3n3n) :: relay, dQ, T, dQ_add, modes
    real(dp), allocatable :: eigs(:)
    type(vecn) :: omega
    integer :: i_xyz, i, i_atom
    integer :: n_negative_eigs, n_atoms
    logical :: do_impl_deriv
    real(dp), allocatable :: c_lambda12i_c(:, :)

    do_impl_deriv = allocated(alpha_0%dr) .or. allocated(C6%dr) .or. &
        allocated(damp%r_vdw%dr) .or. allocated(damp%sigma%dr)
    n_atoms = sys%siz()
    allocate (eigs(3*n_atoms))
    T = dipole_matrix(sys, damp, k_point)
    if (do_impl_deriv) then
        call relay%copy_from(T)
    else
        call relay%move_from(T)
    end if
    omega = omega_eff(C6, alpha_0)
    call relay%mult_cross(omega%val*sqrt(alpha_0%val))
    call relay%add_diag(omega%val**2)
    call ts(sys%calc, 21)
    if (sys%get_modes .or. sys%do_gradients) then
        call modes%alloc_from(relay)
        call eigh(modes, eigs, sys%calc%exc, src=relay)
        if (sys%get_modes) then
            if (allocated(modes%re)) then
                call move_alloc(modes%re, ene%modes)
            else
                call move_alloc(modes%cplx, ene%modes_k_single)
            end if
        end if
    else
        eigs = eigvalsh(relay, sys%calc%exc, destroy=.true.)
    end if
    if (sys%has_exc()) return
    call ts(sys%calc, -21)
    if (sys%get_eigs) then
        ene%mode_enes = sqrt(eigs)
        where (eigs < 0) ene%mode_enes = 0.d0
    end if
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        sys%calc%info%neg_eig = "CDM Hamiltonian has " // &
            trim(tostr(n_negative_eigs)) //  " negative eigenvalues"
        if (sys%calc%param%zero_negative_eigs) where (eigs < 0) eigs = 0.d0
    end if
    ene%energy = 1.d0/2*sum(sqrt(eigs))-3.d0/2*sum(omega%val)
    if (.not. sys%do_gradients) return
    allocate (c_lambda12i_c(3*n_atoms, 3*n_atoms))
    allocate (ene%gradients(n_atoms, 3), source=0.d0)
    forall (i = 1:3*n_atoms)
        c_lambda12i_c(:, i) = eigs(i)**(-1.d0/4)*modes%re(:, i)
    end forall
    c_lambda12i_c = matmul(c_lambda12i_c, transpose(c_lambda12i_c))
    dQ%blacs = T%blacs
    do i_xyz = 1, 3
        dQ%re = -T%re_dr(:, :, i_xyz)
        call dQ%mult_cross(omega%val*sqrt(alpha_0%val))
        dQ%re = c_lambda12i_c*dQ%re
        ene%gradients(:, i_xyz) = 1.d0/4*contract_gradients(dQ)
    end do
    if (.not. do_impl_deriv) return
    dQ_add%blacs = T%blacs
    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            dQ%re(:, :) = 0.d0
            if (allocated(omega%dr)) then
                call dQ%add_diag(2*omega%val*omega%dr(:, i_atom, i_xyz))
                call dQ%add(T%multed_cross( &
                    omega%val*sqrt(alpha_0%val), &
                    omega%dr(:, i_atom, i_xyz)*sqrt(alpha_0%val) &
                ))
            end if
            if (allocated(alpha_0%dr)) then
                call dQ%add(T%multed_cross( &
                    omega%val*sqrt(alpha_0%val), &
                    omega%val*alpha_0%dr(:, i_atom, i_xyz)/(2*sqrt(alpha_0%val)) &
                ))
            end if
            if (allocated(damp%r_vdw%dr)) then
                dQ_add%re = T%re_dvdw
                call dQ_add%mult_cross_add(damp%r_vdw%dr(:, i_atom, i_xyz))
                call dQ_add%mult_cross(omega%val*sqrt(alpha_0%val))
                call dQ%add(dQ_add)
            end if
            ene%gradients(i_atom, i_xyz) = ene%gradients(i_atom, i_xyz) + &
                1.d0/4*sum(c_lambda12i_c*dQ%re)
        end do
    end do
    if (allocated(omega%dr)) then
        ene%gradients = ene%gradients - 3.d0/2*sum(omega%dr, 1)
    end if
end function get_single_mbd_energy


type(mbd_result) function get_reciprocal_mbd_energy(sys, alpha_0, C6, damp) &
        result(ene)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha_0
    type(vecn), intent(in) :: C6
    type(mbd_damping), intent(in) :: damp

    logical :: do_rpa
    integer :: i_kpt, n_kpts, n_atoms
    real(dp) :: k_point(3)
    type(vecn), allocatable :: alpha_ts(:)
    type(mbd_result) :: ene_k

    n_atoms = sys%siz()
    ene%k_pts = make_k_grid(make_g_grid( &
        sys%calc, sys%k_grid(1), sys%k_grid(2), sys%k_grid(3) &
    ), sys%lattice)
    n_kpts = size(ene%k_pts, 2)
    do_rpa = sys%do_rpa

    allocate (alpha_ts(0:ubound(sys%calc%omega_grid, 1)))
    alpha_ts = alpha_dynamic_ts(sys%calc, alpha_0, C6)
    ene%energy = 0.d0
    if (sys%get_eigs) &
        allocate (ene%mode_enes_k(3*n_atoms, n_kpts), source=0.d0)
    if (sys%get_modes) &
        allocate (ene%modes_k(3*n_atoms, 3*n_atoms, n_kpts), source=(0.d0, 0.d0))
    if (sys%get_rpa_orders) allocate ( &
        ene%rpa_orders_k(sys%calc%param%rpa_order_max, n_kpts), source=0.d0 &
    )
    do i_kpt = 1, n_kpts
        k_point = ene%k_pts(:, i_kpt)
        if (do_rpa) then
            ene_k = get_single_reciprocal_rpa_ene(sys, alpha_ts, k_point, damp)
            if (sys%get_rpa_orders) then
                ene%rpa_orders_k(:, i_kpt) = ene_k%rpa_orders
            end if
        else
            ene_k = get_single_mbd_energy(sys, alpha_0, C6, damp, k_point)
            if (sys%get_eigs) ene%mode_enes_k(:, i_kpt) = ene_k%mode_enes
            if (sys%get_modes) ene%modes_k(:, :, i_kpt) = ene_k%modes_k_single
        end if
        ene%energy = ene%energy + ene_k%energy
    end do ! k_point loop
    ene%energy = ene%energy/size(ene%k_pts, 2)
    if (sys%get_rpa_orders) ene%rpa_orders = ene%rpa_orders/n_kpts
end function get_reciprocal_mbd_energy


type(mbd_result) function get_single_rpa_energy(sys, alpha, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha(0:)
    type(mbd_damping), intent(in) :: damp

    type(mat3n3n) :: relay, AT
    complex(dp), allocatable :: eigs(:)
    integer :: i_grid_omega, i, my_i_atom
    integer :: n_order, n_negative_eigs
    type(mbd_damping) :: damp_alpha

    ene%energy = 0.d0
    damp_alpha = damp
    allocate (eigs(3*sys%siz()))
    do i_grid_omega = 0, ubound(sys%calc%omega_grid, 1)
        damp_alpha%sigma = get_sigma_selfint(alpha(i_grid_omega))
        ! relay = T
        relay = dipole_matrix(sys, damp_alpha)
        do my_i_atom = 1, size(relay%blacs%i_atom)
            associate ( &
                    i_atom => relay%blacs%i_atom(my_i_atom), &
                    relay_sub => relay%re(3*(my_i_atom-1)+1:, :) &
            )
                relay_sub(:3, :) = relay_sub(:3, :) * &
                    alpha(i_grid_omega)%val(i_atom)
            end associate
        end do
        ! relay = alpha*T
        if (sys%get_rpa_orders) AT = relay
        ! relay = 1+alpha*T
        call relay%add_diag_scalar(1.d0)
        call ts(sys%calc, 23)
        eigs = eigvals(relay%re, sys%calc%exc, destroy=.true.)
        call ts(sys%calc, -23)
        if (sys%has_exc()) return
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            sys%calc%info%neg_eig = "1+AT matrix has " // &
                trim(tostr(n_negative_eigs)) // " negative eigenvalues"
        end if
        ene%energy = ene%energy + &
            1.d0/(2*pi)*sum(log(dble(eigs)))*sys%calc%omega_grid_w(i_grid_omega)
        if (sys%get_rpa_orders) then
            call ts(sys%calc, 24)
            eigs = eigvals(AT%re, sys%calc%exc, destroy=.true.)
            call ts(sys%calc, -24)
            if (sys%has_exc()) return
            allocate (ene%rpa_orders(sys%calc%param%rpa_order_max))
            do n_order = 2, sys%calc%param%rpa_order_max
                ene%rpa_orders(n_order) = ene%rpa_orders(n_order) &
                    +(-1.d0/(2*pi)*(-1)**n_order &
                    *sum(dble(eigs)**n_order)/n_order) &
                    *sys%calc%omega_grid_w(i_grid_omega)
            end do
        end if
    end do
end function get_single_rpa_energy


type(mbd_result) function get_single_reciprocal_rpa_ene( &
        sys, alpha, k_point, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    type(vecn), intent(in) :: alpha(0:)
    real(dp), intent(in) :: k_point(3)
    type(mbd_damping), intent(in) :: damp

    type(mat3n3n) :: relay, AT
    complex(dp), allocatable :: eigs(:)
    integer :: i_atom, i_grid_omega, i
    integer :: n_order, n_negative_eigs
    type(mbd_damping) :: damp_alpha

    ene%energy = 0d0
    damp_alpha = damp
    allocate (eigs(3*sys%siz()))
    do i_grid_omega = 0, ubound(sys%calc%omega_grid, 1)
        damp_alpha%sigma = get_sigma_selfint(alpha(i_grid_omega))
        ! relay = T
        relay = dipole_matrix(sys, damp_alpha, k_point)
        do i_atom = 1, sys%siz()
            i = 3*(i_atom-1)
            relay%cplx(i+1:i+3, :i) = alpha(i_grid_omega)%val(i_atom) * &
                conjg(transpose(relay%cplx(:i, i+1:i+3)))
        end do
        do i_atom = 1, sys%siz()
            i = 3*(i_atom-1)
            relay%cplx(i+1:i+3, i+1:) = &
                alpha(i_grid_omega)%val(i_atom)*relay%cplx(i+1:i+3, i+1:)
        end do
        ! relay = alpha*T
        if (sys%get_rpa_orders) AT = relay
        do i = 1, 3*sys%siz()
            relay%cplx(i, i) = 1.d0+relay%cplx(i, i) ! relay = 1+alpha*T
        end do
        call ts(sys%calc, 25)
        eigs = eigvals(relay%cplx, sys%calc%exc, destroy=.true.)
        if (sys%has_exc()) return
        call ts(sys%calc, -25)
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            sys%calc%info%neg_eig = "1+AT matrix has " // &
                trim(tostr(n_negative_eigs)) // " negative eigenvalues"
        end if
        ene%energy = ene%energy + &
            1.d0/(2*pi)*dble(sum(log(eigs)))*sys%calc%omega_grid_w(i_grid_omega)
        if (sys%get_rpa_orders) then
            call ts(sys%calc, 26)
            eigs = eigvals(AT%cplx, sys%calc%exc, destroy=.true.)
            if (sys%has_exc()) return
            call ts(sys%calc, -26)
            do n_order = 2, sys%calc%param%rpa_order_max
                ene%rpa_orders(n_order) = ene%rpa_orders(n_order) + &
                    (-1.d0)/(2*pi)*(-1)**n_order * &
                    dble(sum(eigs**n_order))/n_order * &
                    sys%calc%omega_grid_w(i_grid_omega)
            end do
        end if
    end do
end function get_single_reciprocal_rpa_ene


function contract_polarizability(alpha_full) result(alpha)
    type(mat3n3n), intent(in) :: alpha_full
    real(dp) :: alpha(alpha_full%blacs%n_atoms)

    integer :: i_xyz, my_j_atom

    alpha(:) = 0.d0
    do my_j_atom = 1, size(alpha_full%blacs%j_atom)
        associate (j_atom => alpha_full%blacs%j_atom(my_j_atom))
            do i_xyz = 1, 3
                alpha(j_atom) = alpha(j_atom) + &
                    sum(alpha_full%re(i_xyz::3, 3*(my_j_atom-1)+i_xyz))
            end do
        end associate
    end do
    alpha = alpha/3
end function contract_polarizability


function contract_gradients(relay) result(gradient)
    type(mat3n3n), intent(in) :: relay
    real(dp) :: gradient(relay%blacs%n_atoms)

    integer :: my_j_atom

    gradient(:) = 0.d0
    do my_j_atom = 1, size(relay%blacs%j_atom)
        associate ( &
                j_atom => relay%blacs%j_atom(my_j_atom), &
                relay_sub => relay%re(:, 3*(my_j_atom-1)+1:))
            gradient(j_atom) = gradient(j_atom) + 2*sum(relay_sub(:, :3))
        end associate
    end do
end function contract_gradients


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
    type(vecn), intent(in) :: alpha_0
    type(vecn), intent(in) :: C6
    type(vecn) :: alpha(0:ubound(calc%omega_grid, 1))

    integer :: i_freq
    type(vecn) :: omega

    omega = omega_eff(C6, alpha_0)
    do i_freq = 0, ubound(alpha, 1)
        alpha(i_freq) = alpha_osc(alpha_0, omega, calc%omega_grid(i_freq))
    end do
end function


! equation 14
type(vecn) function alpha_osc(alpha_0, omega, u) result(alpha)
    type(vecn), intent(in) :: alpha_0, omega
    real(dp), intent(in) :: u

    integer :: n_atoms, i

    alpha%val = alpha_0%val/(1+(u/omega%val)**2)
    n_atoms = size(alpha_0%val)
    if (allocated(alpha_0%dr) .or. allocated(omega%dr)) then
        allocate (alpha%dr(n_atoms, n_atoms, 3), source=0.d0)
    end if
    if (allocated(alpha_0%dr)) then
        do i = 1, n_atoms
            alpha%dr(i, :, :) = alpha%val(i)*alpha_0%dr(i, :, :)/alpha_0%val(i)
        end do
    end if
    if (allocated(omega%dr)) then
        do i = 1, n_atoms
            alpha%dr(i, :, :) = alpha%dr(i, :, :) + &
                alpha%val(i)*2.d0/omega%val(i)*omega%dr(i, :, :) / &
                (1.d0+(omega%val(i)/u)**2)
        end do
    end if
end function


! equation 13
type(vecn) function scale_TS(X, alpha_0_sc, alpha_0, q) result(X_sc)
    type(vecn), intent(in) :: X, alpha_0_sc, alpha_0
    real(dp), intent(in) :: q

    integer :: i, n_atoms

    X_sc%val = X%val*(alpha_0_sc%val/alpha_0%val)**q
    n_atoms = X%siz()
    if (allocated(X%dr) .or. allocated(alpha_0_sc%dr) .or. &
            allocated(alpha_0%dr)) then
        allocate (X_sc%dr(n_atoms, n_atoms, 3), source=0.d0)
    end if
    if (allocated(X%dr)) then
        do i = 1, n_atoms
            X_sc%dr(i, :, :) = X_sc%val(i)*X%dr(i, :, :)/X%val(i)
        end do
    end if
    if (allocated(alpha_0_sc%dr)) then
        do i = 1, n_atoms
            X_sc%dr(i, :, :) = X_sc%dr(i, :, :) + &
                X_sc%val(i)*q*alpha_0_sc%dr(i, :, :)/alpha_0_sc%val(i)
        end do
    end if
    if (allocated(alpha_0%dr)) then
        do i = 1, n_atoms
            X_sc%dr(i, :, :) = X_sc%dr(i, :, :) - &
                X_sc%val(i)*q*alpha_0%dr(i, :, :)/alpha_0%val(i)
        end do
    end if
end function


elemental function combine_C6(C6_i, C6_j, alpha_0_i, alpha_0_j) result(C6_ij)
    real(dp), intent(in) :: C6_i, C6_j, alpha_0_i, alpha_0_j
    real(dp) :: C6_ij

    C6_ij = 2*C6_i*C6_j/(alpha_0_j/alpha_0_i*C6_i+alpha_0_i/alpha_0_j*C6_j)
end function


! equation 12
type(vecn) function omega_eff(C6, alpha) result(omega)
    type(vecn), intent(in) :: C6, alpha

    integer :: i, n_atoms

    omega%val = 4.d0/3*C6%val/alpha%val**2
    n_atoms = C6%siz()
    if (allocated(C6%dr) .or. allocated(alpha%dr)) then
        allocate (omega%dr(n_atoms, n_atoms, 3), source=0.d0)
    end if
    if (allocated(C6%dr)) then
        do i = 1, n_atoms
            omega%dr(i, :, :) = omega%val(i)*C6%dr(i, :, :)/C6%val(i)
        end do
    end if
    if (allocated(alpha%dr)) then
        do i = 1, n_atoms
            omega%dr(i, :, :) = omega%dr(i, :, :) - &
                2*omega%val(i)*alpha%dr(i, :, :)/alpha%val(i)
        end do
    end if
end function


type(vecn) function get_sigma_selfint(alpha) result(sigma)
    type(vecn), intent(in) :: alpha

    integer :: i, n_atoms

    sigma%val = (sqrt(2.d0/pi)*alpha%val/3.d0)**(1.d0/3)
    if (allocated(alpha%dr)) then
        n_atoms = alpha%siz()
        allocate (sigma%dr(n_atoms, n_atoms, 3))
        do i = 1, n_atoms
            sigma%dr(i, :, :) = sigma%val(i)*alpha%dr(i, :, :)/(3*alpha%val(i))
        end do
    end if
end function


type(vecn) function get_C6_from_alpha(calc, alpha) result(C6)
    type(mbd_calc), intent(in) :: calc
    type(vecn), intent(in) :: alpha(0:)

    integer :: i_atom, i_freq, n_atoms

    n_atoms = alpha(0)%siz()
    allocate (C6%val(n_atoms), source=0.d0)
    do i_freq = 0, ubound(alpha, 1)
        C6%val = C6%val + &
            3.d0/pi*alpha(i_freq)%val**2*calc%omega_grid_w(i_freq)
    end do
    if (allocated(alpha(0)%dr)) then
        allocate (C6%dr(n_atoms, n_atoms, 3), source=0.d0)
        do i_freq = 0, ubound(alpha, 1)
            do i_atom = 1, n_atoms
                C6%dr(i_atom, :, :) = C6%dr(i_atom, :, :) + &
                    6.d0/pi*alpha(i_freq)%val(i_atom) * &
                    alpha(i_freq)%dr(i_atom, :, :)*calc%omega_grid_w(i_freq)
            end do
        end do
    end if
end function


function supercell_circum(sys, uc, radius) result(sc)
    type(mbd_system), intent(in) :: sys
    real(dp), intent(in) :: uc(3, 3), radius
    integer :: sc(3)

    real(dp) :: ruc(3, 3), layer_sep(3)
    integer :: i

    ruc = 2*pi*inverse(transpose(uc))
    forall (i = 1:3) &
        layer_sep(i) = sum(uc(:, i)*ruc(:, i)/sqrt(sum(ruc(:, i)**2)))
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
    real(dp) :: g_grid(3, n1*n2*n3)

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
        g_grid(:, i_kpt) = g_kpt_shifted/kpt_range
    end do
end function make_g_grid


logical function mbd_system_has_exc(sys)
    class(mbd_system), intent(in) :: sys

    mbd_system_has_exc = sys%calc%exc%label /= ''
end function


function make_k_grid(g_grid, uc) result(k_grid)
    real(dp), intent(in) :: g_grid(:, :), uc(3, 3)
    real(dp) :: k_grid(3, size(g_grid, 2))

    integer :: i_kpt
    real(dp) :: ruc(3, 3)

    ruc = 2*pi*inverse(transpose(uc))
    do i_kpt = 1, size(g_grid, 2)
        k_grid(:, i_kpt) = matmul(ruc, g_grid(:, i_kpt))
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
            calc%tm%timestamps(calc%tm%ts_aid) = &
                calc%tm%timestamps(calc%tm%ts_aid)+calc%tm%ts_cnt
            calc%tm%ts_counts(calc%tm%ts_aid) = &
                calc%tm%ts_counts(calc%tm%ts_aid)+1
        end if
    end if
end subroutine ts


function clock_rate() result(rate)
    integer :: cnt, rate, cnt_max

    call system_clock(cnt, rate, cnt_max)
end function clock_rate


integer function system_siz(this)
    class(mbd_system), intent(in) :: this

    if (allocated(this%coords)) then
        system_siz = size(this%coords, 2)
    else
        system_siz = 0
    end if
end function


end module mbd
