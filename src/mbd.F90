! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd

use mbd_common, only: tostr, print_matrix, dp, pi, mbd_exc, findval, lower, &
    MBD_EXC_NEG_EIGVALS, MBD_EXC_NEG_POL, MBD_EXC_UNIMPL, printer, &
    shift_cell
use mbd_system_type, only: mbd_system, mbd_calc, ang
use mbd_linalg, only: invh, inverse, eigh, eigvals, eigvalsh, outer, mmul
use mbd_types, only: mat3n3n, mat33, scalar, contract_cross_33
use mbd_parallel, only: mbd_blacs_grid, mbd_blacs, all_reduce
use mbd_defaults

implicit none

#ifndef MODULE_UNIT_TESTS
private
public :: mbd_damping, mbd_result, mbd_energy, dipole_matrix, mbd_scs_energy, &
    sigma_selfint, scale_TS, set_damping_parameters, &
    mbd_gradients, damping_fermi, test_frequency_grid
#endif

type :: mbd_damping
    character(len=20) :: version
    real(dp) :: beta = 0d0
    real(dp) :: a = MBD_DAMPING_A
    real(dp) :: ts_d = TS_DAMPING_D
    real(dp) :: ts_sr = 0d0
    real(dp) :: mayer_scaling = 1d0
    real(dp), allocatable :: r_vdw(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: damping_custom(:, :)
    real(dp), allocatable :: potential_custom(:, :, :, :)
end type mbd_damping

type :: mbd_result
    real(dp) :: energy
    real(dp), allocatable :: k_pts(:, :)
    real(dp), allocatable :: mode_eigs(:)
    real(dp), allocatable :: modes(:, :)
    real(dp), allocatable :: rpa_orders(:)
    real(dp), allocatable :: mode_eigs_k(:, :)
    complex(dp), allocatable :: modes_k(:, :, :)
    complex(dp), allocatable :: modes_k_single(:, :)
    real(dp), allocatable :: rpa_orders_k(:, :)
end type

type :: mbd_gradients
    real(dp), allocatable :: dcoords(:, :)  ! n_atoms by 3
    real(dp), allocatable :: dalpha(:)
    real(dp), allocatable :: dalpha_dyn(:, :)  ! n_atoms by 0:n_freq
    real(dp), allocatable :: dC6(:)
    real(dp), allocatable :: dr_vdw(:)
    real(dp), allocatable :: domega(:)
    real(dp), allocatable :: dV(:)
    real(dp), allocatable :: dV_free(:)
    real(dp), allocatable :: dX_free(:)
    contains
    procedure :: copy_alloc => gradients_copy_alloc
    procedure :: has_grad => gradients_has_grad
end type

contains

type(mbd_result) function mbd_scs_energy( &
        sys, variant, alpha_0, C6, damp, dene) result(res)
    type(mbd_system), intent(inout) :: sys
    character(len=*), intent(in) :: variant
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp
    type(mbd_gradients), intent(inout) :: dene

    real(dp), allocatable :: alpha_dyn(:, :), alpha_dyn_scs(:, :), &
        C6_scs(:), dC6_scs_dalpha_dyn_scs(:, :), &
        dene_dalpha_scs_dyn(:, :), freq_w(:)
    type(mbd_gradients), allocatable :: dalpha_dyn(:), dalpha_dyn_scs(:, :)
    type(mbd_gradients) :: dene_mbd, dr_vdw_scs
    type(mbd_damping) :: damp_scs, damp_mbd
    integer :: n_freq, i_freq, n_atoms, i_atom, my_i_atom
    character(len=15) :: damping_types(2)

    select case (variant)
    case ('scs')
        damping_types = [character(len=15) :: 'dip,gg', 'dip,1mexp']
    case ('rsscs')
        damping_types = [character(len=15) :: 'fermi,dip,gg', 'fermi,dip']
    end select
    n_freq = ubound(sys%calc%omega_grid, 1)
    n_atoms = sys%siz()
    call allocate_derivs()
    allocate (alpha_dyn(n_atoms, 0:n_freq))
    allocate (alpha_dyn_scs(n_atoms, 0:n_freq))
    alpha_dyn = alpha_dynamic_ts(sys%calc, alpha_0, C6, dalpha_dyn)
    damp_scs = damp
    damp_scs%version = damping_types(1)
    do i_freq = 0, n_freq
        alpha_dyn_scs(:, i_freq) = run_scs( &
            sys, alpha_dyn(:, i_freq), damp_scs, dalpha_dyn_scs(:, i_freq) &
        )
        if (sys%has_exc()) return
    end do
    C6_scs = get_C6_from_alpha(sys%calc, alpha_dyn_scs, dC6_scs_dalpha_dyn_scs)
    damp_mbd = damp
    damp_mbd%r_vdw = scale_TS( &
        damp%r_vdw, alpha_dyn_scs(:, 0), alpha_dyn(:, 0), 1d0/3, dr_vdw_scs &
    )
    damp_mbd%version = damping_types(2)
    res = mbd_energy(sys, alpha_dyn_scs(:, 0), C6_scs, damp_mbd, dene_mbd)
    if (sys%has_exc()) return
    if (.not. dene%has_grad()) return
    freq_w = sys%calc%omega_grid_w
    freq_w(0) = 1d0
    dene_dalpha_scs_dyn(:, 0) = dene_mbd%dalpha + dene_mbd%dr_vdw*dr_vdw_scs%dV
    do i_freq = 1, n_freq
        dene_dalpha_scs_dyn(:, i_freq) = &
            dene_mbd%dC6*dC6_scs_dalpha_dyn_scs(:, i_freq)
    end do
    if (allocated(dene%dcoords)) then
        dene%dcoords = 0d0
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = sys%blacs%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dcoords(sys%blacs%j_atom, :) = &
                    dene%dcoords(sys%blacs%j_atom, :) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dcoords
            end do
        end do
        call all_reduce(dene%dcoords, sys%blacs)
        dene%dcoords = dene%dcoords + dene_mbd%dcoords
    end if
    if (allocated(dene%dalpha)) then
        dene%dalpha = 0d0
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = sys%blacs%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dalpha(sys%blacs%j_atom) = dene%dalpha(sys%blacs%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dalpha * &
                    dalpha_dyn(i_freq)%dalpha(sys%blacs%j_atom)
            end do
        end do
        call all_reduce(dene%dalpha, sys%blacs)
        dene%dalpha = dene%dalpha + dene_mbd%dr_vdw*dr_vdw_scs%dV_free
    end if
    if (allocated(dene%dC6)) then
        dene%dC6 = 0d0
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = sys%blacs%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dC6(sys%blacs%j_atom) = dene%dC6(sys%blacs%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dalpha * &
                    dalpha_dyn(i_freq)%dC6(sys%blacs%j_atom)
            end do
        end do
        call all_reduce(dene%dC6, sys%blacs)
    end if
    if (allocated(dene%dr_vdw)) then
        dene%dr_vdw = 0d0
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = sys%blacs%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                dene%dr_vdw(sys%blacs%j_atom) = dene%dr_vdw(sys%blacs%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dr_vdw
            end do
        end do
        call all_reduce(dene%dr_vdw, sys%blacs)
        dene%dr_vdw = dene%dr_vdw + dene_mbd%dr_vdw*dr_vdw_scs%dX_free
    end if

    contains

    subroutine allocate_derivs()
        integer :: my_ncatoms

        if (allocated(dene%dcoords)) allocate (dene_mbd%dcoords(n_atoms, 3))
        if (dene%has_grad()) then
            allocate (dC6_scs_dalpha_dyn_scs(n_atoms, 0:n_freq))
            allocate (dr_vdw_scs%dV(n_atoms))
            allocate (dene_dalpha_scs_dyn(n_atoms, 0:n_freq))
            allocate (dene_mbd%dalpha(n_atoms))
            allocate (dene_mbd%dC6(n_atoms))
            allocate (dene_mbd%dr_vdw(n_atoms))
        end if
        if (allocated(dene%dalpha)) allocate (dr_vdw_scs%dV_free(n_atoms))
        if (allocated(dene%dr_vdw)) allocate (dr_vdw_scs%dX_free(n_atoms))
        allocate (dalpha_dyn(0:n_freq))
        allocate (dalpha_dyn_scs(size(sys%blacs%i_atom), 0:n_freq))
        my_ncatoms = size(sys%blacs%j_atom)
        do i_freq = 0, n_freq
            if (allocated(dene%dalpha)) &
                allocate (dalpha_dyn(i_freq)%dalpha(n_atoms))
            if (allocated(dene%dC6)) &
                allocate (dalpha_dyn(i_freq)%dC6(n_atoms))
            do my_i_atom = 1, size(dalpha_dyn_scs, 1)
                associate (da => dalpha_dyn_scs(my_i_atom, i_freq))
                    if (allocated(dene%dcoords)) &
                        allocate (da%dcoords(my_ncatoms, 3))
                    if (allocated(dene%dalpha) .or. allocated(dene%dC6)) &
                        allocate (da%dalpha(my_ncatoms))
                    if (allocated(dene%dr_vdw)) &
                        allocate (da%dr_vdw(my_ncatoms))
                end associate
            end do
        end do
    end subroutine
end function mbd_scs_energy

type(mat3n3n) function dipole_matrix(sys, damp, grad, k_point) result(dipmat)
    type(mbd_system), intent(inout) :: sys
    type(mbd_damping), intent(in) :: damp
    logical, intent(in) :: grad
    real(dp), intent(in), optional :: k_point(3)

    real(dp) :: R_cell(3), r(3), r_norm, R_vdw_ij, &
        sigma_ij, volume, ewald_alpha, real_space_cutoff, f_ij
    type(mat33) :: Tpp
    complex(dp) :: Tpp_c(3, 3)
    integer :: i_atom, j_atom, i_cell, idx_cell(3), range_cell(3), i, j, &
        n_atoms, my_i_atom, my_j_atom, my_nratoms, my_ncatoms
    logical :: do_ewald, is_periodic

    do_ewald = .false.
    is_periodic = allocated(sys%lattice)
    n_atoms = sys%siz()
    call dipmat%init(sys%blacs)
    my_nratoms = size(dipmat%blacs%i_atom)
    my_ncatoms = size(dipmat%blacs%j_atom)
    if (present(k_point)) then
        allocate (dipmat%cplx(3*my_nratoms, 3*my_ncatoms), source=(0d0, 0d0))
    else
        allocate (dipmat%re(3*my_nratoms, 3*my_ncatoms), source=0d0)
        if (grad) then
            allocate (dipmat%re_dr(3*my_nratoms, 3*my_ncatoms, 3), source=0d0)
            allocate (dipmat%re_dvdw(3*my_nratoms, 3*my_ncatoms), source=0d0)
            allocate (dipmat%re_dsigma(3*my_nratoms, 3*my_ncatoms), source=0d0)
        end if
    end if
    ! MPI code end
    if (is_periodic) then
        if (any(sys%vacuum_axis)) then
            real_space_cutoff = sys%calc%param%dipole_low_dim_cutoff
        else if (sys%calc%param%ewald_on) then
            if (grad) then
                sys%calc%exc%code = MBD_EXC_UNIMPL
                sys%calc%exc%msg = 'Forces not implemented for periodic systems'
                return
            end if
            do_ewald = .true.
            volume = max(abs(dble(product(eigvals(sys%lattice)))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1d0/3)
            real_space_cutoff = &
                6d0/ewald_alpha*sys%calc%param%ewald_real_cutoff_scaling
            sys%calc%info%ewald_alpha = &
                'Ewald: using alpha = ' // trim(tostr(ewald_alpha)) // &
                ', real cutoff = ' // trim(tostr(real_space_cutoff))
        else
            real_space_cutoff = sys%calc%param%dipole_cutoff
        end if
        range_cell = sys%supercell_circum(sys%lattice, real_space_cutoff)
    else
        range_cell(:) = 0
    end if
    if (is_periodic) then
        sys%calc%info%ewald_rsum = &
            'Ewald: summing real part in cell vector range of ' // &
            trim(tostr(1+2*range_cell(1))) // 'x' // &
            trim(tostr(1+2*range_cell(2))) // 'x' // &
            trim(tostr(1+2*range_cell(3)))
    end if
    call sys%clock(11)
    idx_cell = [0, 0, -1]
    do i_cell = 1, product(1+2*range_cell)
        call shift_cell(idx_cell, -range_cell, range_cell)
        if (is_periodic) then
            R_cell = matmul(sys%lattice, idx_cell)
        else
            R_cell(:) = 0d0
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
                if (is_periodic .and. r_norm > real_space_cutoff) cycle
                if (allocated(damp%R_vdw)) then
                    R_vdw_ij = sum(damp%R_vdw([i_atom, j_atom]))
                end if
                if (allocated(damp%sigma)) then
                    sigma_ij = damp%mayer_scaling * &
                        sqrt(sum(damp%sigma([i_atom, j_atom])**2))
                end if
                select case (damp%version)
                    case ("bare")
                        Tpp = T_bare_v2(r, grad)
                    case ("dip,1mexp")
                        Tpp%val = T_1mexp_coulomb(r, damp%beta*R_vdw_ij, damp%a)
                    case ("fermi,dip")
                        Tpp = T_damped(damping_fermi( &
                            r, damp%beta*R_vdw_ij, damp%a, grad &
                        ), T_bare_v2(r, grad))
                    case ("sqrtfermi,dip")
                        Tpp = T_damped(damping_sqrtfermi( &
                            r, damp%beta*R_vdw_ij, damp%a, grad &
                        ), T_bare_v2(r, grad))
                    case ("custom,dip")
                        Tpp%val = damp%damping_custom(i_atom, j_atom)*T_bare(r)
                    case ("dip,custom")
                        Tpp%val = damp%potential_custom(:, :, i_atom, j_atom)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma_ij, grad)
                    case ("fermi,dip,gg")
                        Tpp = T_damped(op1minus(damping_fermi( &
                            r, damp%beta*R_vdw_ij, damp%a, grad &
                        )), T_erf_coulomb(r, sigma_ij, grad))
                        do_ewald = .false.
                    case ("sqrtfermi,dip,gg")
                        Tpp = T_damped(op1minus(damping_sqrtfermi( &
                            r, damp%beta*R_vdw_ij, damp%a, grad &
                        )), T_erf_coulomb(r, sigma_ij, grad))
                        do_ewald = .false.
                    case ("custom,dip,gg")
                        f_ij = 1d0-damp%damping_custom(i_atom, j_atom)
                        Tpp = T_erf_coulomb(r, sigma_ij, grad)
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
                    Tpp_c = Tpp%val*exp(-cmplx(0d0, 1d0, 8)*( &
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
                    if (allocated(Tpp%dr)) then
                        associate (T => dipmat%re_dr(i+1:i+3, j+1:j+3, :))
                            T = T + Tpp%dr
                        end associate
                    end if
                    if (allocated(Tpp%dvdw)) then
                        associate (dTdRvdw => dipmat%re_dvdw(i+1:i+3, j+1:j+3))
                            dTdRvdw = dTdRvdw + Tpp%dvdw
                        end associate
                    end if
                    if (allocated(Tpp%dsigma)) then
                        associate (dTdsigma => dipmat%re_dsigma(i+1:i+3, j+1:j+3))
                            dTdsigma = dTdsigma + Tpp%dsigma
                        end associate
                    end if
                end if
            end do ! j_atom
        end do ! i_atom
    end do ! i_cell
    call sys%clock(-11)
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
    rec_space_cutoff = 10d0*alpha*sys%calc%param%ewald_rec_cutoff_scaling
    range_G_vector = sys%supercell_circum(rec_unit_cell, rec_space_cutoff)
    sys%calc%info%ewald_cutoff = 'Ewald: using reciprocal cutoff = ' // &
        trim(tostr(rec_space_cutoff))
    sys%calc%info%ewald_recsum = &
        'Ewald: summing reciprocal part in G vector range of ' // &
        trim(tostr(1+2*range_G_vector(1))) // 'x' // &
        trim(tostr(1+2*range_G_vector(2))) // 'x' // &
        trim(tostr(1+2*range_G_vector(3)))
    call sys%clock(12)
    idx_G_vector = [0, 0, -1]
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
        k_prefactor = 4*pi/volume*exp(-k_sq/(4*alpha**2))
        forall (i_xyz = 1:3, j_xyz = 1:3) &
                k_prefactor(i_xyz, j_xyz) = k_prefactor(i_xyz, j_xyz) &
                *k_total(i_xyz)*k_total(j_xyz)/k_sq
        do my_i_atom = 1, size(dipmat%blacs%i_atom)
            i_atom = dipmat%blacs%i_atom(my_i_atom)
            do my_j_atom = 1, size(dipmat%blacs%j_atom)
                j_atom = dipmat%blacs%j_atom(my_j_atom)
                r = sys%coords(:, i_atom)-sys%coords(:, j_atom)
                if (present(k_point)) then
                    Tpp_c = k_prefactor*exp(cmplx(0d0, 1d0, 8) &
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
    call sys%clock(-12)
end subroutine

subroutine test_frequency_grid(calc)
    type(mbd_calc), intent(inout) :: calc

    real(dp) :: alpha(1, 0:ubound(calc%omega_grid, 1)), C6(1), error
    type(mbd_gradients) :: dalpha(0:ubound(calc%omega_grid, 1))
    real(dp), allocatable :: dC6_dalpha(:, :)

    alpha = alpha_dynamic_ts(calc, [21d0], [99.5d0], dalpha)
    C6 = get_C6_from_alpha(calc, alpha, dC6_dalpha)
    error = abs(C6(1)/99.5d0-1d0)
    calc%info%freq_error = &
        "Relative quadrature error in C6 of carbon atom: " // &
        trim(tostr(error))
end subroutine

function run_scs(sys, alpha, damp, dalpha_scs) result(alpha_scs)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha(:)
    type(mbd_damping), intent(in) :: damp
    type(mbd_gradients), intent(inout) :: dalpha_scs(:)
    real(dp) :: alpha_scs(size(alpha))

    type(mat3n3n) :: alpha_full, dQ, T
    integer :: n_atoms, i_xyz, i_atom, my_i_atom
    type(mbd_damping) :: damp_local
    real(dp), allocatable :: dsij_dsi(:), dsigma_dalpha(:), &
        alpha_prime(:, :), B_prime(:, :), grads_i(:)

    n_atoms = sys%siz()
    damp_local = damp
    if (allocated(dalpha_scs(1)%dalpha)) allocate (dsigma_dalpha(n_atoms))
    damp_local%sigma = sigma_selfint(alpha, dsigma_dalpha)
    T = dipole_matrix(sys, damp_local, dalpha_scs(1)%has_grad())
    if (sys%has_exc()) return
    if (dalpha_scs(1)%has_grad()) then
        call alpha_full%copy_from(T)
    else
        call alpha_full%move_from(T)
    end if
    call alpha_full%add_diag(1d0/alpha)
    call sys%clock(32)
    call invh(alpha_full, sys%calc%exc)
    if (sys%has_exc()) return
    call sys%clock(-32)
    alpha_scs = alpha_full%contract_n33diag_cols()
    if (any(alpha_scs < 0)) then
        sys%calc%exc%code = MBD_EXC_NEG_POL
        sys%calc%exc%msg = 'Screening leads to negative polarizability'
        return
    end if
    if (.not. dalpha_scs(1)%has_grad()) return
    allocate (alpha_prime(3, 3*n_atoms), B_prime(3*n_atoms, 3), source=0d0)
    allocate (grads_i(n_atoms))
    call alpha_full%contract_n_transp('R', alpha_prime)
    call dQ%init_from(T)
    if (allocated(dalpha_scs(1)%dcoords)) then
        do i_xyz = 1, 3
            dQ%re = -T%re_dr(:, :, i_xyz)
            dQ = mmul(alpha_full, dQ)
            call dQ%contract_n_transp('C', B_prime)
            do i_atom = 1, n_atoms
                grads_i = contract_cross_33( &
                    i_atom, dQ, alpha_prime, alpha_full, B_prime &
                )
                my_i_atom = findval(sys%blacs%i_atom, i_atom)
                if (my_i_atom > 0) then
                    dalpha_scs(my_i_atom)%dcoords(:, i_xyz) = &
                        grads_i(sys%blacs%j_atom)
                end if
            end do
        end do
    end if
    if (allocated(dalpha_scs(1)%dalpha)) then
        dQ%re = T%re_dsigma
        do i_atom = 1, n_atoms
            dsij_dsi = damp_local%sigma(i_atom)*dsigma_dalpha(i_atom) / &
                sqrt(damp_local%sigma(i_atom)**2+damp_local%sigma**2)
            call dQ%mult_col(i_atom, dsij_dsi)
        end do
        call dQ%add_diag(-0.5d0/alpha**2)
        dQ = mmul(alpha_full, dQ)
        call dQ%contract_n_transp('C', B_prime)
        do i_atom = 1, n_atoms
            grads_i = contract_cross_33( &
                i_atom, dQ, alpha_prime, alpha_full, B_prime &
            )
            my_i_atom = findval(sys%blacs%i_atom, i_atom)
            if (my_i_atom > 0) then
                dalpha_scs(my_i_atom)%dalpha = grads_i(sys%blacs%j_atom)
            end if
        end do
    end if
    if (allocated(dalpha_scs(1)%dr_vdw)) then
        dQ%re = T%re_dvdw
        dQ = mmul(alpha_full, dQ)
        call dQ%contract_n_transp('C', B_prime)
        do i_atom = 1, n_atoms
            grads_i = contract_cross_33( &
                i_atom, dQ, alpha_prime, alpha_full, B_prime &
            )
            my_i_atom = findval(sys%blacs%i_atom, i_atom)
            if (my_i_atom > 0) then
                dalpha_scs(my_i_atom)%dr_vdw = grads_i(sys%blacs%j_atom)
            end if
        end do
    end if
end function run_scs

type(mbd_result) function mbd_energy(sys, alpha_0, C6, damp, dene) result(res)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp
    type(mbd_gradients), intent(inout) :: dene

    real(dp), allocatable :: alpha(:, :)
    type(mbd_gradients), allocatable :: dalpha(:)
    integer :: n_freq, n_atoms, i_freq

    if (.not. allocated(sys%lattice)) then
        if (.not. sys%do_rpa) then
            res = get_single_mbd_energy(sys, alpha_0, C6, damp, dene)
        else
            n_freq = ubound(sys%calc%omega_grid, 1)
            n_atoms = sys%siz()
            allocate (alpha(n_atoms, 0:n_freq), dalpha(0:n_freq))
            do i_freq = 0, n_freq
                call dene%copy_alloc(dalpha(i_freq))
            end do
            alpha = alpha_dynamic_ts(sys%calc, alpha_0, C6, dalpha)
            res = get_single_rpa_energy(sys, alpha, damp)
            ! TODO gradients
        end if
    else
        res = get_reciprocal_mbd_energy(sys, alpha_0, C6, damp)
    end if
end function mbd_energy

type(mbd_result) function get_single_mbd_energy( &
        sys, alpha_0, C6, damp, dene, k_point) result(res)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp
    type(mbd_gradients), intent(inout) :: dene
    real(dp), intent(in), optional :: k_point(3)

    type(mat3n3n) :: relay, dQ, T, modes, c_lambda12i_c
    real(dp), allocatable :: eigs(:), omega(:)
    type(mbd_gradients) :: domega
    integer :: i_xyz, n_negative_eigs, n_atoms
    logical :: grad
    character(120) :: msg

    n_atoms = sys%siz()
    grad = dene%has_grad()
    T = dipole_matrix(sys, damp, grad, k_point)
    if (sys%has_exc()) return
    if (grad) then
        call relay%copy_from(T)
    else
        call relay%move_from(T)
    end if
    call dene%copy_alloc(domega)
    omega = omega_eff(C6, alpha_0, domega)
    call relay%mult_cross(omega*sqrt(alpha_0))
    call relay%add_diag(omega**2)
    call sys%clock(21)
    if (sys%get_modes .or. grad) then
        call modes%alloc_from(relay)
        allocate (eigs(3*n_atoms))
        call eigh(modes, eigs, sys%calc%exc, src=relay)
        if (sys%get_modes) then
            if (allocated(modes%re)) then
                call move_alloc(modes%re, res%modes)
            else
                call move_alloc(modes%cplx, res%modes_k_single)
            end if
        end if
    else
        eigs = eigvalsh(relay, sys%calc%exc, destroy=.true.)
    end if
    if (sys%has_exc()) return
    call sys%clock(-21)
    if (sys%get_eigs) res%mode_eigs = eigs
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        msg = "CDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues"
        if (sys%calc%param%zero_negative_eigs) then
            where (eigs < 0) eigs = 0d0
            sys%calc%info%neg_eigvals = msg
        else
            sys%calc%exc%code = MBD_EXC_NEG_EIGVALS
            sys%calc%exc%msg = msg
        end if
    end if
    res%energy = 1d0/2*sum(sqrt(eigs))-3d0/2*sum(omega)
    if (.not. grad) return
    call c_lambda12i_c%copy_from(modes)
    call c_lambda12i_c%mult_cols_3n(eigs**(-1d0/4))
    c_lambda12i_c = mmul(c_lambda12i_c, c_lambda12i_c, transB=.true.)
    call dQ%init_from(T)
    if (allocated(dene%dcoords)) then
        do i_xyz = 1, 3
            dQ%re = T%re_dr(:, :, i_xyz)
            call dQ%mult_cross(omega*sqrt(alpha_0))
            dQ%re = c_lambda12i_c%re*dQ%re
            dene%dcoords(:, i_xyz) = 1d0/2*dQ%contract_n33_rows()
        end do
    end if
    if (allocated(dene%dalpha)) then
        dQ%re = T%re
        call dQ%mult_cross(omega*sqrt(alpha_0))
        call dQ%mult_rows(1d0/(2*alpha_0)+domega%dalpha/omega)
        call dQ%add_diag(omega*domega%dalpha)
        dQ%re = c_lambda12i_c%re*dQ%re
        dene%dalpha = 1d0/2*dQ%contract_n33_rows()-3d0/2*domega%dalpha
    end if
    if (allocated(dene%dC6)) then
        dQ%re = T%re
        call dQ%mult_cross(omega*sqrt(alpha_0))
        call dQ%mult_rows(domega%dC6/omega)
        call dQ%add_diag(omega*domega%dC6)
        dQ%re = c_lambda12i_c%re*dQ%re
        dene%dC6 = 1d0/2*dQ%contract_n33_rows()-3d0/2*domega%dC6
    end if
    if (allocated(dene%dr_vdw)) then
        dQ%re = T%re_dvdw
        call dQ%mult_cross(omega*sqrt(alpha_0))
        dQ%re = c_lambda12i_c%re*dQ%re
        dene%dr_vdw = 1d0/2*dQ%contract_n33_rows()
    end if
end function get_single_mbd_energy

type(mbd_result) function get_reciprocal_mbd_energy(sys, alpha_0, C6, damp) result(res)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_damping), intent(in) :: damp

    logical :: do_rpa
    integer :: i_kpt, n_kpts, n_atoms, n_freq
    real(dp) :: k_point(3)
    real(dp), allocatable :: alpha_ts(:, :)
    type(mbd_gradients), allocatable :: dalpha_ts(:)
    type(mbd_result) :: res_k
    type(mbd_gradients) :: dene_k

    n_atoms = sys%siz()
    res%k_pts = make_k_grid(make_g_grid( &
        sys%calc, sys%k_grid(1), sys%k_grid(2), sys%k_grid(3) &
    ), sys%lattice)
    n_kpts = size(res%k_pts, 2)
    n_freq = ubound(sys%calc%omega_grid, 1)
    do_rpa = sys%do_rpa

    allocate (alpha_ts(n_atoms, 0:n_freq), dalpha_ts(0:n_freq))
    alpha_ts = alpha_dynamic_ts(sys%calc, alpha_0, C6, dalpha_ts)
    res%energy = 0d0
    if (sys%get_eigs) &
        allocate (res%mode_eigs_k(3*n_atoms, n_kpts), source=0d0)
    if (sys%get_modes) &
        allocate (res%modes_k(3*n_atoms, 3*n_atoms, n_kpts), source=(0d0, 0d0))
    if (sys%get_rpa_orders) allocate ( &
        res%rpa_orders_k(sys%calc%param%rpa_order_max, n_kpts), source=0d0 &
    )
    do i_kpt = 1, n_kpts
        k_point = res%k_pts(:, i_kpt)
        if (do_rpa) then
            res_k = get_single_reciprocal_rpa_ene(sys, alpha_ts, k_point, damp)
            if (sys%get_rpa_orders) then
                res%rpa_orders_k(:, i_kpt) = res_k%rpa_orders
            end if
        else
            res_k = get_single_mbd_energy(sys, alpha_0, C6, damp, dene_k, k_point)
            if (sys%get_eigs) res%mode_eigs_k(:, i_kpt) = res_k%mode_eigs
            if (sys%get_modes) res%modes_k(:, :, i_kpt) = res_k%modes_k_single
        end if
        if (sys%has_exc()) return
        res%energy = res%energy + res_k%energy
    end do ! k_point loop
    res%energy = res%energy/size(res%k_pts, 2)
    if (sys%get_rpa_orders) res%rpa_orders = res%rpa_orders/n_kpts
end function get_reciprocal_mbd_energy

type(mbd_result) function get_single_rpa_energy(sys, alpha, damp) result(res)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha(:, 0:)
    type(mbd_damping), intent(in) :: damp

    type(mat3n3n) :: relay, AT
    complex(dp), allocatable :: eigs(:)
    integer :: i_freq, i, my_i_atom, n_order, n_negative_eigs
    type(mbd_damping) :: damp_alpha
    real(dp), allocatable :: dsigma_dalpha(:)

    res%energy = 0d0
    damp_alpha = damp
    allocate (eigs(3*sys%siz()))
    do i_freq = 0, ubound(sys%calc%omega_grid, 1)
        damp_alpha%sigma = sigma_selfint(alpha(:, i_freq), dsigma_dalpha)
        ! relay = T
        relay = dipole_matrix(sys, damp_alpha, .false.)
        do my_i_atom = 1, size(relay%blacs%i_atom)
            associate ( &
                    i_atom => relay%blacs%i_atom(my_i_atom), &
                    relay_sub => relay%re(3*(my_i_atom-1)+1:, :) &
            )
                relay_sub(:3, :) = relay_sub(:3, :)*alpha(i_atom, i_freq)
            end associate
        end do
        ! relay = alpha*T
        if (sys%get_rpa_orders) AT = relay
        ! relay = 1+alpha*T
        call relay%add_diag_scalar(1d0)
        call sys%clock(23)
        eigs = eigvals(relay, sys%calc%exc, destroy=.true.)
        call sys%clock(-23)
        if (sys%has_exc()) return
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            sys%calc%exc%code = MBD_EXC_NEG_EIGVALS
            sys%calc%exc%msg = "1+AT matrix has " // &
                trim(tostr(n_negative_eigs)) // " negative eigenvalues"
            return
        end if
        res%energy = res%energy + &
            1d0/(2*pi)*sum(log(dble(eigs)))*sys%calc%omega_grid_w(i_freq)
        if (sys%get_rpa_orders) then
            call sys%clock(24)
            eigs = eigvals(AT, sys%calc%exc, destroy=.true.)
            call sys%clock(-24)
            if (sys%has_exc()) return
            allocate (res%rpa_orders(sys%calc%param%rpa_order_max))
            do n_order = 2, sys%calc%param%rpa_order_max
                res%rpa_orders(n_order) = res%rpa_orders(n_order) &
                    +(-1d0/(2*pi)*(-1)**n_order &
                    *sum(dble(eigs)**n_order)/n_order) &
                    *sys%calc%omega_grid_w(i_freq)
            end do
        end if
    end do
end function get_single_rpa_energy

type(mbd_result) function get_single_reciprocal_rpa_ene(sys, alpha, k_point, damp) &
        result(res)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: alpha(0:, :)
    real(dp), intent(in) :: k_point(3)
    type(mbd_damping), intent(in) :: damp

    type(mat3n3n) :: relay, AT
    complex(dp), allocatable :: eigs(:)
    integer :: i_atom, i_freq, i, n_order, n_negative_eigs
    type(mbd_damping) :: damp_alpha
    real(dp), allocatable :: dsigma_dalpha(:)

    res%energy = 0d0
    damp_alpha = damp
    allocate (eigs(3*sys%siz()))
    do i_freq = 0, ubound(sys%calc%omega_grid, 1)
        damp_alpha%sigma = sigma_selfint(alpha(:, i_freq), dsigma_dalpha)
        ! relay = T
        relay = dipole_matrix(sys, damp_alpha, .false., k_point)
        do i_atom = 1, sys%siz()
            i = 3*(i_atom-1)
            relay%cplx(i+1:i+3, :i) = alpha(i_freq, i_atom) * &
                conjg(transpose(relay%cplx(:i, i+1:i+3)))
        end do
        do i_atom = 1, sys%siz()
            i = 3*(i_atom-1)
            relay%cplx(i+1:i+3, i+1:) = &
                alpha(i_freq, i_atom)*relay%cplx(i+1:i+3, i+1:)
        end do
        ! relay = alpha*T
        if (sys%get_rpa_orders) AT = relay
        do i = 1, 3*sys%siz()
            relay%cplx(i, i) = 1d0+relay%cplx(i, i) ! relay = 1+alpha*T
        end do
        call sys%clock(25)
        eigs = eigvals(relay%cplx, sys%calc%exc, destroy=.true.)
        if (sys%has_exc()) return
        call sys%clock(-25)
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            sys%calc%exc%code = MBD_EXC_NEG_EIGVALS
            sys%calc%exc%msg = "1+AT matrix has " // &
                trim(tostr(n_negative_eigs)) // " negative eigenvalues"
            return
        end if
        res%energy = res%energy + &
            1d0/(2*pi)*dble(sum(log(eigs)))*sys%calc%omega_grid_w(i_freq)
        if (sys%get_rpa_orders) then
            call sys%clock(26)
            eigs = eigvals(AT%cplx, sys%calc%exc, destroy=.true.)
            if (sys%has_exc()) return
            call sys%clock(-26)
            do n_order = 2, sys%calc%param%rpa_order_max
                res%rpa_orders(n_order) = res%rpa_orders(n_order) + &
                    (-1d0)/(2*pi)*(-1)**n_order * &
                    dble(sum(eigs**n_order))/n_order * &
                    sys%calc%omega_grid_w(i_freq)
            end do
        end if
    end do
end function get_single_reciprocal_rpa_ene

function T_bare(rxyz) result(T)
    real(dp), intent(in) :: rxyz(3)
    real(dp) :: T(3, 3)

    integer :: i, j
    real(dp) :: r_sq, r_5

    r_sq = sum(rxyz(:)**2)
    r_5 = sqrt(r_sq)**5
    do i = 1, 3
        T(i, i) = (3d0*rxyz(i)**2-r_sq)/r_5
        do j = i+1, 3
            T(i, j) = 3d0*rxyz(i)*rxyz(j)/r_5
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

    C = (3*erfc(a*r)+(2*a*r/sqrt(pi))*(3d0+2*(a*r)**2)*exp(-(a*r)**2))/r**5
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
    f%val = 1d0/(1+exp(-d*(eta-1)))
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

type(mat33) function T_damped(f, T) result(fT)
    type(scalar), intent(in) :: f
    type(mat33), intent(in) :: T

    integer :: c

    fT%val = f%val*T%val
    if (allocated(f%dr) .or. allocated(T%dr)) &
        allocate (fT%dr(3, 3, 3), source=0d0)
    if (allocated(f%dvdw) .or. allocated(T%dvdw)) &
        allocate (fT%dvdw(3, 3), source=0d0)
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
    rr_r5 = outer(r, r)/r_5
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
    zeta_1 = 1d0-exp(-r_sigma)-a*r_sigma*exp(-r_sigma)
    zeta_2 = -r_sigma*a*exp(-r_sigma)*(1+a*(-1+r_sigma))
    T = zeta_1*T_bare(rxyz)-zeta_2*outer(rxyz, rxyz)/sqrt(sum(rxyz**2))**5
end function

subroutine set_damping_parameters(xc, ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, &
        mbd_ts_erf_beta, mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta)
    character(len=*), intent(in) :: xc
    real(dp), intent(out) :: &
        ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, mbd_ts_erf_beta, &
        mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta

    ts_d = 20d0
    ts_s_r = 1d0
    mbd_scs_a = 2d0
    mbd_ts_a = 6d0
    mbd_ts_erf_beta = 1d0
    mbd_ts_fermi_beta = 1d0
    mbd_rsscs_a = 6d0
    mbd_rsscs_beta = 1d0
    select case (lower(xc))
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
end subroutine set_damping_parameters

elemental function terf(r, r0, a)
    real(dp), intent(in) :: r, r0, a
    real(dp) :: terf

    terf = 0.5d0*(erf(a*(r+r0))+erf(a*(r-r0)))
end function

function alpha_dynamic_ts(calc, alpha_0, C6, dalpha) result(alpha)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(mbd_gradients), intent(inout) :: dalpha(0:)
    real(dp) :: alpha(size(alpha_0), 0:ubound(calc%omega_grid, 1))

    integer :: i_freq, n_atoms
    real(dp), allocatable :: omega(:)
    type(mbd_gradients) :: domega

    n_atoms = size(alpha_0)
    call dalpha(0)%copy_alloc(domega)
    omega = omega_eff(C6, alpha_0, domega)
    do i_freq = 0, ubound(alpha, 2)
        if (allocated(dalpha(i_freq)%dalpha) .or. &
                allocated(dalpha(i_freq)%dC6)) then
            allocate (dalpha(i_freq)%domega(n_atoms))
        end if
        alpha(:, i_freq) = alpha_osc(&
            alpha_0, omega, calc%omega_grid(i_freq), dalpha(i_freq) &
        )
        if (allocated(dalpha(i_freq)%dalpha)) then
            dalpha(i_freq)%dalpha = dalpha(i_freq)%dalpha + &
                dalpha(i_freq)%domega*domega%dalpha
        end if
        if (allocated(dalpha(i_freq)%dC6)) then
            dalpha(i_freq)%dC6 = dalpha(i_freq)%domega*domega%dC6
        end if
        if (allocated(dalpha(i_freq)%domega)) deallocate(dalpha(i_freq)%domega)
    end do
end function

! equation 14
function alpha_osc(alpha_0, omega, u, dalpha) result(alpha)
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: omega(:)
    real(dp), intent(in) :: u
    type(mbd_gradients), intent(inout) :: dalpha
    real(dp) :: alpha(size(alpha_0))

    alpha = alpha_0/(1+(u/omega)**2)
    if (allocated(dalpha%dalpha)) dalpha%dalpha = alpha/alpha_0
    if (allocated(dalpha%domega)) dalpha%domega = &
        alpha*2d0/omega/(1d0+(omega/u)**2)
end function

! equation 13
function scale_TS(X_free, V, V_free, q, dX) result(X)
    real(dp), intent(in) :: X_free(:), V(:), V_free(:)
    real(dp), intent(in) :: q
    type(mbd_gradients), intent(inout) :: dX
    real(dp) :: X(size(X_free))

    X = X_free*(V/V_free)**q
    if (allocated(dX%dX_free)) dX%dX_free = X/X_free
    if (allocated(dX%dV)) dX%dV = X*q/V
    if (allocated(dX%dV_free)) dX%dV_free = -X*q/V_free
end function



! equation 12
function omega_eff(C6, alpha, domega) result(omega)
    real(dp), intent(in) :: C6(:)
    real(dp), intent(in) :: alpha(:)
    type(mbd_gradients), intent(inout) :: domega
    real(dp) :: omega(size(C6))

    omega = 4d0/3*C6/alpha**2
    if (allocated(domega%dC6)) domega%dC6 = omega/C6
    if (allocated(domega%dalpha)) domega%dalpha = -2*omega/alpha
end function

function sigma_selfint(alpha, dsigma_dalpha) result(sigma)
    real(dp), intent(in) :: alpha(:)
    real(dp), intent(inout), allocatable :: dsigma_dalpha(:)
    real(dp) :: sigma(size(alpha))

    sigma = (sqrt(2d0/pi)*alpha/3d0)**(1d0/3)
    if (allocated(dsigma_dalpha)) dsigma_dalpha = sigma/(3*alpha)
end function

function get_C6_from_alpha(calc, alpha, dC6_dalpha) result(C6)
    type(mbd_calc), intent(in) :: calc
    real(dp), intent(in) :: alpha(:, 0:)
    real(dp), intent(inout), allocatable :: dC6_dalpha(:, :)
    real(dp) :: C6(size(alpha, 1))

    integer :: i_freq, n_atoms

    n_atoms = size(alpha, 1)
    C6 = 0d0
    do i_freq = 0, ubound(alpha, 2)
        C6 = C6 + 3d0/pi*alpha(:, i_freq)**2*calc%omega_grid_w(i_freq)
    end do
    if (.not. allocated(dC6_dalpha)) return
    dC6_dalpha = 0d0
    do i_freq = 0, ubound(alpha, 2)
        dC6_dalpha(:, i_freq) = dC6_dalpha(:, i_freq) + 6d0/pi*alpha(:, i_freq)
    end do
end function

function make_g_grid(calc, n1, n2, n3) result(g_grid)
    type(mbd_calc), intent(in) :: calc
    integer, intent(in) :: n1, n2, n3
    real(dp) :: g_grid(3, n1*n2*n3)

    integer :: g_kpt(3), i_kpt, kpt_range(3)
    real(dp) :: g_kpt_shifted(3)

    g_kpt = [0, 0, -1]
    kpt_range = [n1, n2, n3]
    do i_kpt = 1, n1*n2*n3
        call shift_cell (g_kpt, [0, 0, 0], kpt_range-1)
        g_kpt_shifted = dble(g_kpt)+calc%param%k_grid_shift
        where (2*g_kpt_shifted > kpt_range)
            g_kpt_shifted = g_kpt_shifted-dble(kpt_range)
        end where
        g_grid(:, i_kpt) = g_kpt_shifted/kpt_range
    end do
end function make_g_grid

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

subroutine gradients_copy_alloc(this, other)
    class(mbd_gradients), intent(in) :: this
    type(mbd_gradients), intent(out) :: other

    if (allocated(this%dcoords)) &
        allocate (other%dcoords(size(this%dcoords, 2), 3))
    if (allocated(this%dalpha)) &
        allocate (other%dalpha(size(this%dalpha)))
    if (allocated(this%dC6)) allocate (other%dC6(size(this%dC6)))
    if (allocated(this%dr_vdw)) allocate (other%dr_vdw(size(this%dr_vdw)))
    if (allocated(this%domega)) allocate (other%domega(size(this%domega)))
end subroutine

logical function gradients_has_grad(this)
    class(mbd_gradients), intent(in) :: this

    gradients_has_grad = allocated(this%dcoords) .or. &
        allocated(this%dalpha) .or. allocated(this%dC6) .or. &
        allocated(this%dr_vdw) .or. allocated(this%domega)
end function

end module mbd
