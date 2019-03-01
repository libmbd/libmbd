! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_methods
!! Obtaining MBD energies.

use mbd_constants
use mbd_damping, only: damping_t
use mbd_formulas, only: omega_qho, alpha_dyn_qho, scale_with_ratio, C6_from_alpha
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_t, grad_request_t
use mbd_hamiltonian, only: get_mbd_hamiltonian_energy
use mbd_lapack, only: eigvals, inverse
use mbd_rpa, only: get_mbd_rpa_energy
use mbd_scs, only: run_scs
use mbd_utils, only: result_t, tostr, quad_pt_t, shift_idx
#ifdef WITH_SCALAPACK
use mbd_blacs, only: blacs_all_reduce
#endif
#ifdef WITH_MPI
use mbd_mpi, only: mpi_all_reduce
#endif

implicit none

private
public :: get_mbd_energy, get_mbd_scs_energy

contains

type(result_t) function get_mbd_energy(geom, alpha_0, C6, damp, grad) result(res)
    !! Get MBD energy.
    !!
    !! For a nonperiodic system, the method just transforms \(C_6\) coefficients
    !! to frequencies, and performs a single call to
    !! [[get_mbd_hamiltonian_energy]]. For a periodic system, the method
    !! integrates the energy over the frist Brillouin zone.
    !!
    !! $$
    !! E=\int_\text{FBZ}\mathrm d\mathbf q\,E(\mathbf
    !! q)\approx\frac1{N_k}\sum_i^{N_k}E(\mathbf q_i)
    !! \\ \mathbf q_i=\boldsymbol{\mathcal B}\mathbf n_i,\qquad\partial\mathbf
    !! q_i=-\big((\partial\boldsymbol{\mathcal
    !! A})\boldsymbol{\mathcal A}^{-1}\big)^\mathrm T\mathbf q_i
    !! $$
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    type(grad_request_t), intent(in) :: grad

    real(dp), allocatable :: alpha(:, :), omega(:), k_pts(:, :), dkdlattice(:, :, :, :)
    type(grad_t), allocatable :: dalpha(:)
    integer :: n_kpts, i_kpt, a
    type(result_t) :: res_k
    type(grad_t) :: domega
    type(grad_request_t) :: grad_ham

    omega = omega_qho(C6, alpha_0, domega, grad)
    if (geom%do_rpa) then
        alpha = alpha_dyn_qho(alpha_0, omega, geom%freq, dalpha, grad_request_t())
    end if
    grad_ham = grad
    if (grad%dC6 .or. grad%dalpha) grad_ham%domega = .true.
    if (grad%dlattice) grad_ham%dq = .true.
    if (.not. allocated(geom%lattice)) then
        if (.not. geom%do_rpa) then
            res = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, grad_ham)
            if (grad%dC6) res%dE%dC6 = res%dE%domega*domega%dC6
            if (grad%dalpha) res%dE%dalpha = res%dE%dalpha + res%dE%domega*domega%dalpha
            if (allocated(res%dE%domega)) deallocate (res%dE%domega)
        else
            res = get_mbd_rpa_energy(geom, alpha, damp)
            ! TODO gradients
        end if
    else
        call make_k_pts( &
            k_pts, geom%k_grid, geom%lattice, geom%param%k_grid_shift, &
            dkdlattice, grad%dlattice &
        )
        n_kpts = size(k_pts, 2)
        res%energy = 0d0
        if (geom%get_eigs) &
            allocate (res%mode_eigs_k(3*geom%siz(), n_kpts), source=0d0)
        if (geom%get_modes) &
            allocate (res%modes_k(3*geom%siz(), 3*geom%siz(), n_kpts), source=(0d0, 0d0))
        if (geom%get_rpa_orders) allocate ( &
            res%rpa_orders_k(geom%param%rpa_order_max, n_kpts), source=0d0 &
        )
        if (grad%dcoords) allocate (res%dE%dcoords(geom%siz(), 3), source=0d0)
        if (grad%dlattice) allocate (res%dE%dlattice(3, 3), source=0d0)
        if (grad%dalpha) allocate (res%dE%dalpha(geom%siz()), source=0d0)
        if (grad%dC6) allocate (res%dE%dC6(geom%siz()), source=0d0)
        if (grad%dR_vdw) allocate (res%dE%dR_vdw(geom%siz()), source=0d0)
        do i_kpt = 1, n_kpts
#ifdef WITH_MPI
            if (geom%parallel_mode == 'k_points') then
                if (modulo(i_kpt, geom%mpi_size) /= geom%mpi_rank) cycle
            end if
#endif
            call geom%clock(51)
            associate (k_pt => k_pts(:, i_kpt))
                if (.not. geom%do_rpa) then
                    res_k = get_mbd_hamiltonian_energy( &
                        geom, alpha_0, omega, damp, grad_ham, k_pt &
                    )
                else
                    res_k = get_mbd_rpa_energy(geom, alpha, damp, k_pt)
                end if
            end associate
            call geom%clock(-51)
            if (geom%has_exc()) return
            if (geom%get_eigs) then
                res%mode_eigs_k(:, i_kpt) = res_k%mode_eigs
            end if
            if (geom%get_modes) then
                res%modes_k(:, :, i_kpt) = res_k%modes_k_single
            end if
            if (geom%get_rpa_orders) then
                res%rpa_orders_k(:, i_kpt) = res_k%rpa_orders
            end if
            res%energy = res%energy + res_k%energy/n_kpts
            if (grad%dcoords) res%dE%dcoords = res%dE%dcoords + res_k%dE%dcoords/n_kpts
            if (grad%dlattice) then
                res%dE%dlattice = res%dE%dlattice + res_k%dE%dlattice/n_kpts
                do a = 1, 3
                    res%dE%dlattice = res%dE%dlattice &
                        + res_k%dE%dq(a)*dkdlattice(a, i_kpt, :, :)/n_kpts
                end do
            end if
            if (grad%dalpha) then
                res%dE%dalpha = res%dE%dalpha &
                    + (res_k%dE%dalpha + res_k%dE%domega*domega%dalpha)/n_kpts
            end if
            if (grad%dC6) res%dE%dC6 = res%dE%dC6 + res_k%dE%domega*domega%dC6/n_kpts
            if (grad%dR_vdw) res%dE%dR_vdw = res%dE%dR_vdw + res_k%dE%dR_vdw/n_kpts
        end do
#ifdef WITH_MPI
        if (geom%parallel_mode == 'k_points') then
            call mpi_all_reduce(res%energy, geom%mpi_comm)
            if (grad%dcoords) call mpi_all_reduce(res%dE%dcoords, geom%mpi_comm)
            if (grad%dlattice) call mpi_all_reduce(res%dE%dlattice, geom%mpi_comm)
            if (grad%dalpha) call mpi_all_reduce(res%dE%dalpha, geom%mpi_comm)
            if (grad%dC6) call mpi_all_reduce(res%dE%dC6, geom%mpi_comm)
            if (grad%dR_vdw) call mpi_all_reduce(res%dE%dR_vdw, geom%mpi_comm)
        end if
#endif
    end if
end function

type(result_t) function get_mbd_scs_energy(geom, variant, alpha_0, C6, damp, grad) result(res)
    !! Get screened MBD energy.
    type(geom_t), intent(inout) :: geom
    character(len=*), intent(in) :: variant
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    type(grad_request_t), intent(in) :: grad

    real(dp), allocatable :: alpha_dyn(:, :), alpha_dyn_scs(:, :), &
        C6_scs(:), dC6_scs_dalpha_dyn_scs(:, :), &
        dene_dalpha_scs_dyn(:, :), freq_w(:), omega(:)
    type(grad_t), allocatable :: dalpha_dyn(:), dalpha_dyn_scs(:, :)
    type(grad_t) :: dene_mbd, dr_vdw_scs, domega
    type(grad_request_t) :: grad_scs
    type(damping_t) :: damp_scs, damp_mbd
    integer :: n_freq, i_freq, n_atoms, i_atom, my_i_atom
    character(len=15) :: damping_types(2)

    call geom%clock(90)
    select case (variant)
    case ('scs')
        damping_types = [character(len=15) :: 'dip,gg', 'dip,1mexp']
    case ('rsscs')
        damping_types = [character(len=15) :: 'fermi,dip,gg', 'fermi,dip']
    end select
    n_freq = ubound(geom%freq, 1)
    n_atoms = geom%siz()
    allocate (alpha_dyn(n_atoms, 0:n_freq))
    allocate (alpha_dyn_scs(n_atoms, 0:n_freq))
    allocate (dalpha_dyn_scs(size(geom%idx%i_atom), 0:n_freq))
    if (grad%any()) allocate (dene_dalpha_scs_dyn(n_atoms, 0:n_freq))
    omega = omega_qho(C6, alpha_0, domega, grad)
    alpha_dyn = alpha_dyn_qho( &
        alpha_0, omega, geom%freq, dalpha_dyn, &
        grad_request_t(dalpha=grad%dalpha, domega=grad%dalpha .or. grad%dC6) &
    )
    grad_scs = grad_request_t( &
        dcoords=grad%dcoords, &
        dlattice=grad%dlattice, &
        dalpha=grad%dalpha .or. grad%dC6, &
        dr_vdw=grad%dr_vdw &
    )
    damp_scs = damp
    damp_scs%version = damping_types(1)
    call geom%clock(50)
    do i_freq = 0, n_freq
        alpha_dyn_scs(:, i_freq) = run_scs( &
            geom, alpha_dyn(:, i_freq), damp_scs, dalpha_dyn_scs(:, i_freq), grad_scs &
        )
        if (geom%has_exc()) return
    end do
    call geom%clock(-50)
    C6_scs = C6_from_alpha(alpha_dyn_scs, geom%freq, dC6_scs_dalpha_dyn_scs, grad%any())
    damp_mbd = damp
    damp_mbd%r_vdw = scale_with_ratio( &
        damp%r_vdw, alpha_dyn_scs(:, 0), alpha_dyn(:, 0), 1d0/3, dr_vdw_scs, &
        grad_request_t(dV=grad%any(), dV_free=grad%dalpha, dX_free=grad%dr_vdw) &
    )
    damp_mbd%version = damping_types(2)
    res = get_mbd_energy(geom, alpha_dyn_scs(:, 0), C6_scs, damp_mbd, &
        grad_request_t( &
            dcoords=grad%dcoords, dlattice=grad%dlattice, &
            dalpha=grad%any(), dC6=grad%any(), dr_vdw=grad%any() &
        ) &
    )
    dene_mbd = res%dE
    res%dE = grad_t()
    call geom%clock(-90)
    if (geom%has_exc()) return
    if (.not. grad%any()) return
    call geom%clock(91)
    allocate (freq_w(0:ubound(geom%freq, 1)))
    freq_w = geom%freq%weight
    freq_w(0) = 1d0
    dene_dalpha_scs_dyn(:, 0) = dene_mbd%dalpha + dene_mbd%dr_vdw*dr_vdw_scs%dV
    do i_freq = 1, n_freq
        dene_dalpha_scs_dyn(:, i_freq) = &
            dene_mbd%dC6*dC6_scs_dalpha_dyn_scs(:, i_freq)
    end do
    if (grad%dcoords) then
        allocate (res%dE%dcoords(n_atoms, 3), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                res%dE%dcoords(geom%idx%j_atom, :) &
                    = res%dE%dcoords(geom%idx%j_atom, :) &
                    + freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) &
                    * dalpha_dyn_scs(my_i_atom, i_freq)%dcoords
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call blacs_all_reduce(res%dE%dcoords, geom%blacs)
#endif
        res%dE%dcoords = res%dE%dcoords + dene_mbd%dcoords
    end if
    if (grad%dlattice) then
        allocate (res%dE%dlattice(3, 3), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            if (.not. any(i_atom == geom%idx%j_atom)) cycle
            do i_freq = 0, n_freq
                res%dE%dlattice = res%dE%dlattice &
                    + freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) &
                    * dalpha_dyn_scs(my_i_atom, i_freq)%dlattice
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call blacs_all_reduce(res%dE%dlattice, geom%blacs)
#endif
        res%dE%dlattice = res%dE%dlattice + dene_mbd%dlattice
    end if
    if (grad%dalpha) then
        allocate (res%dE%dalpha(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                res%dE%dalpha(geom%idx%j_atom) = res%dE%dalpha(geom%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dalpha * ( &
                        dalpha_dyn(i_freq)%dalpha(geom%idx%j_atom) &
                        + dalpha_dyn(i_freq)%domega(geom%idx%j_atom) &
                        * domega%dalpha(geom%idx%j_atom) &
                    )
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call blacs_all_reduce(res%dE%dalpha, geom%blacs)
#endif
        res%dE%dalpha = res%dE%dalpha + dene_mbd%dr_vdw*dr_vdw_scs%dV_free
    end if
    if (grad%dC6) then
        allocate (res%dE%dC6(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                res%dE%dC6(geom%idx%j_atom) = res%dE%dC6(geom%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dalpha * &
                    dalpha_dyn(i_freq)%domega(geom%idx%j_atom) &
                    * domega%dC6(geom%idx%j_atom)
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call blacs_all_reduce(res%dE%dC6, geom%blacs)
#endif
    end if
    if (grad%dr_vdw) then
        allocate (res%dE%dr_vdw(n_atoms), source=0d0)
        do my_i_atom = 1, size(dalpha_dyn_scs, 1)
            i_atom = geom%idx%i_atom(my_i_atom)
            do i_freq = 0, n_freq
                res%dE%dr_vdw(geom%idx%j_atom) = res%dE%dr_vdw(geom%idx%j_atom) + &
                    freq_w(i_freq)*dene_dalpha_scs_dyn(i_atom, i_freq) * &
                    dalpha_dyn_scs(my_i_atom, i_freq)%dr_vdw
            end do
        end do
#ifdef WITH_SCALAPACK
        if (geom%idx%parallel) call blacs_all_reduce(res%dE%dr_vdw, geom%blacs)
#endif
        res%dE%dr_vdw = res%dE%dr_vdw + dene_mbd%dr_vdw*dr_vdw_scs%dX_free
    end if
    call geom%clock(-91)
end function

real(dp) function test_frequency_grid(freq) result(error)
    !! Calculate relative quadrature error in C6 of a carbon atom
    type(quad_pt_t), intent(in) :: freq(:)

    real(dp) :: alpha(1, 0:ubound(freq, 1)), C6(1)
    type(grad_t), allocatable :: dalpha(:)
    type(grad_request_t) :: grad

    alpha = alpha_dyn_qho([21d0], [99.5d0], freq, dalpha, grad)
    C6 = C6_from_alpha(alpha, freq)
    error = abs(C6(1)/99.5d0-1d0)
end function

! This used to be a function returning the k_pts array, but that was causing
! segfaults with some compilers. I suspect some combination of the product()
! in the dimension specification and assignemnt to allocatable array.
subroutine make_k_pts(k_pts, k_grid, lattice, shift, dkdlattice, grad)
    real(dp), allocatable, intent(out) :: k_pts(:, :)
    integer, intent(in) :: k_grid(3)
    real(dp), intent(in) :: lattice(3, 3)
    real(dp), intent(in) :: shift
    real(dp), allocatable, intent(out) :: dkdlattice(:, :, :, :)
    logical, intent(in) :: grad

    integer :: n_kpt(3), i_kpt, i_latt, a, n_kpts
    real(dp) :: n_kpt_shifted(3), latt_inv(3, 3)

    n_kpts = product(k_grid)
    allocate (k_pts(3, n_kpts))
    n_kpt = [0, 0, -1]
    do i_kpt = 1, n_kpts
        call shift_idx(n_kpt, [0, 0, 0], k_grid-1)
        n_kpt_shifted = dble(n_kpt)+shift
        where (2*n_kpt_shifted > k_grid) n_kpt_shifted = n_kpt_shifted-dble(k_grid)
        k_pts(:, i_kpt) = n_kpt_shifted/k_grid
    end do
    latt_inv = inverse(lattice)
    k_pts = matmul(2*pi*transpose(latt_inv), k_pts)
    if (grad) then
        allocate (dkdlattice(3, n_kpts, 3, 3))
        forall (i_kpt = 1:n_kpts, i_latt = 1:3, a = 1:3)
            dkdlattice(:, i_kpt, i_latt, a) = &
                -latt_inv(i_latt, :)*k_pts(a, i_kpt)
        end forall
    end if
end subroutine

end module
