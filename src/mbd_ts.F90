! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Workaround for bad OpenMP ifort: github.com/libmbd/libmbd/issues/39
#ifndef IS_BAD_OPENMP_IFORT
#define IS_BAD_OPENMP_IFORT defined(__INTEL_COMPILER) && (__INTEL_COMPILER > 1700) && (__INTEL_COMPILER <= 2021)
#endif

module mbd_ts
!! Obtaining TS energies.

use mbd_constants
use mbd_utils, only: shift_idx, tostr, result_t, diff3
use mbd_damping, only: damping_t, damping_fermi
use mbd_geom, only: geom_t, supercell_circum
use mbd_gradients, only: grad_request_t, grad_scalar_t
use mbd_lapack, only: eigvals, inverse
#ifdef WITH_MPI
use mbd_mpi, only: mpi_all_reduce
#endif

implicit none

private
public :: get_ts_energy

contains

type(result_t) function get_ts_energy(geom, alpha_0, C6, damp, grad) result(res)
    !! Get TS energy.
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    type(grad_request_t), intent(in) :: grad

    real(dp) :: C6_ij, Rnij(3), Rnij_norm, R_vdw_ij, ene_ij, Rn(3), f_damp
    integer :: i_cell, i_atom, j_atom, range_n(3), n(3), n_atoms, i_latt
    logical :: is_periodic, do_ewald
    type(grad_request_t) :: grad_ij
    type(grad_scalar_t) :: df, dC6, dphi, dene_ij

    do_ewald = .false.
    is_periodic = allocated(geom%lattice)
    n_atoms = geom%siz()
    grad_ij = grad
    grad_ij%dcoords = grad%dcoords .or. grad%dlattice
    if (is_periodic) then
        do_ewald = geom%gamm > 0d0
        range_n = supercell_circum(geom%lattice, geom%real_space_cutoff)
    else
        range_n(:) = 0
    end if
    if (grad%dcoords) allocate (res%dE%dcoords(n_atoms, 3), source=0d0)
    if (grad%dlattice) allocate (res%dE%dlattice(3, 3), source=0d0)
    if (grad%dC6) allocate (res%dE%dC6(n_atoms), source=0d0)
    if (grad%dalpha) allocate (res%dE%dalpha(n_atoms), source=0d0)
    if (grad%dr_vdw) allocate (res%dE%dr_vdw(n_atoms), source=0d0)
    res%energy = 0d0
    n = [0, 0, -1]
    each_cell: do i_cell = 1, product(1 + 2 * range_n)
        call shift_idx(n, -range_n, range_n)
        if (is_periodic) then
            Rn = matmul(geom%lattice, n)
        else
            Rn(:) = 0d0
        end if
        each_atom: do i_atom = 1, geom%siz()
#ifdef WITH_MPI
            if (modulo(i_atom, geom%mpi_size) /= geom%mpi_rank) cycle
#endif
            each_atom_pair: do j_atom = 1, i_atom
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                Rnij = geom%coords(:, i_atom) - geom%coords(:, j_atom) - Rn
                Rnij_norm = sqrt(sum(Rnij**2))
                if (is_periodic .and. Rnij_norm > geom%real_space_cutoff) cycle
                C6_ij = combine_C6( &
                    C6(i_atom), C6(j_atom), alpha_0(i_atom), alpha_0(j_atom), dC6, grad &
                )
                if (allocated(damp%r_vdw)) then
                    R_vdw_ij = damp%r_vdw(i_atom) + damp%r_vdw(j_atom)
                end if
                select case (damp%version)
                    case ("fermi")
                        f_damp = damping_fermi( &
                            Rnij, damp%ts_sr * R_vdw_ij, damp%ts_d, df, grad_ij &
                        )
                    case ("fermi2")
                        f_damp = damping_fermi(Rnij, damp%ts_sr * R_vdw_ij, damp%ts_d)**2
                    case ("custom")
                        f_damp = damp%damping_custom(i_atom, j_atom)
                end select
                ene_ij = -C6_ij * f_damp / Rnij_norm**6
                if (grad_ij%dcoords) &
                    dene_ij%dr = ene_ij * (df%dr / f_damp - 6 * Rnij / Rnij_norm**2)
                if (grad_ij%dr_vdw) dene_ij%dvdw = ene_ij / f_damp * df%dvdw * damp%ts_sr
                if (do_ewald) then
                    ene_ij = ene_ij - C6_ij * ( &
                        disp_real(Rnij_norm, geom%gamm, dphi, grad_ij) &
                        - 1d0 / Rnij_norm**6 &
                    )
                    if (grad_ij%dcoords) dene_ij%dr = dene_ij%dr &
                        - C6_ij * (dphi%dr_1 + 6 / Rnij_norm**7) * Rnij / Rnij_norm
                end if
                if (i_atom == j_atom) then
                    ene_ij = ene_ij / 2
                    if (grad_ij%dcoords) dene_ij%dr = dene_ij%dr / 2
                    if (grad_ij%dr_vdw) dene_ij%dvdw = dene_ij%dvdw / 2
                end if
                res%energy = res%energy + ene_ij
                if (.not. grad%any()) cycle
                if (grad%dcoords) then
                    res%dE%dcoords(i_atom, :) = res%dE%dcoords(i_atom, :) + dene_ij%dr
                    res%dE%dcoords(j_atom, :) = res%dE%dcoords(j_atom, :) - dene_ij%dr
                end if
                if (grad%dlattice) then
                    do concurrent(i_latt=1:3)
                        res%dE%dlattice(i_latt, :) = res%dE%dlattice(i_latt, :) &
                            - dene_ij%dr * n(i_latt)
                    end do
                end if
                if (grad%dC6) then
                    res%dE%dC6(i_atom) = res%dE%dC6(i_atom) + ene_ij / C6_ij * dC6%dC6i
                    res%dE%dC6(j_atom) = res%dE%dC6(j_atom) + ene_ij / C6_ij * dC6%dC6j
                end if
                if (grad%dalpha) then
                    res%dE%dalpha(i_atom) = res%dE%dalpha(i_atom) + ene_ij / C6_ij * dC6%da0i
                    res%dE%dalpha(j_atom) = res%dE%dalpha(j_atom) + ene_ij / C6_ij * dC6%da0j
                end if
                if (grad%dr_vdw) then
                    res%dE%dr_vdw(i_atom) = res%dE%dr_vdw(i_atom) + dene_ij%dvdw
                    res%dE%dr_vdw(j_atom) = res%dE%dr_vdw(j_atom) + dene_ij%dvdw
                end if
            end do each_atom_pair
        end do each_atom
    end do each_cell
    if (do_ewald) call add_ewald_ts_parts(geom, alpha_0, C6, res, grad)
#ifdef WITH_MPI
    call mpi_all_reduce(res%energy, geom%mpi_comm)
    if (grad%dcoords) call mpi_all_reduce(res%dE%dcoords, geom%mpi_comm)
    if (grad%dlattice) call mpi_all_reduce(res%dE%dlattice, geom%mpi_comm)
    if (grad%dalpha) call mpi_all_reduce(res%dE%dalpha, geom%mpi_comm)
    if (grad%dC6) call mpi_all_reduce(res%dE%dC6, geom%mpi_comm)
    if (grad%dR_vdw) call mpi_all_reduce(res%dE%dR_vdw, geom%mpi_comm)
#endif
end function

subroutine add_ewald_ts_parts(geom, alpha_0, C6, res, grad)
    type(geom_t), intent(in) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(result_t), intent(inout) :: res
    type(grad_request_t), intent(in) :: grad

    real(dp) :: rec_latt(3, 3), volume, Rij(3), k(3), phi, dkdAk_proj, &
        k_norm, k_Rij, latt_inv(3, 3), C6_ij, exp_kR, ene_ij, dkdA(3)
    integer :: i_atom, j_atom, m(3), i_m, range_m(3), i_latt, i_xyz
    type(grad_scalar_t) :: dC6, dphi, dene_ij

    latt_inv = inverse(geom%lattice)
    rec_latt = 2 * pi * transpose(latt_inv)
    volume = abs(dble(product(eigvals(geom%lattice))))
    range_m = supercell_circum(rec_latt, geom%rec_space_cutoff)
    m = [0, 0, -1]
    each_recip_vec: do i_m = 1, product(1 + 2 * range_m)
        call shift_idx(m, -range_m, range_m)
        k = matmul(rec_latt, m)
        k_norm = sqrt(sum(k**2))
        if (k_norm > geom%rec_space_cutoff) cycle
        each_atom: do i_atom = 1, geom%siz()
#ifdef WITH_MPI
            if (modulo(i_atom, geom%mpi_size) /= geom%mpi_rank) cycle
#endif
            each_atom_pair: do j_atom = 1, i_atom
                C6_ij = combine_C6( &
                    C6(i_atom), C6(j_atom), alpha_0(i_atom), alpha_0(j_atom), dC6, grad &
                )
                Rij = geom%coords(:, i_atom) - geom%coords(:, j_atom)
                k_Rij = dot_product(k, Rij)
                exp_kR = cos(k_Rij)
                phi = disp_rec(k_norm, geom%gamm, dphi, grad)
                ene_ij = -C6_ij * phi / volume * exp_kR
                if (i_atom == j_atom) ene_ij = ene_ij / 2
                res%energy = res%energy + ene_ij
                if (grad%dcoords .and. i_atom /= j_atom) then
                    dene_ij%dr = ene_ij / exp_kR * sin(k_Rij) * k
                    res%dE%dcoords(i_atom, :) = res%dE%dcoords(i_atom, :) - dene_ij%dr
                    res%dE%dcoords(j_atom, :) = res%dE%dcoords(j_atom, :) + dene_ij%dr
                end if
                if (grad%dlattice) then
#ifdef IS_BAD_OPENMP_IFORT
                    do i_latt=1, 3
                        do i_xyz=1, 3
#else
                    do concurrent(i_latt=1:3, i_xyz=1:3)
#endif
                            dkdA = -latt_inv(i_latt, :) * k(i_xyz)
                            if (k_norm > 0d0) then
                                dkdAk_proj = dot_product(dkdA, k) / k_norm
                            else
                                dkdAk_proj = 0d0
                            end if
                            res%dE%dlattice(i_latt, i_xyz) = res%dE%dlattice(i_latt, i_xyz) &
                                - ene_ij * latt_inv(i_latt, i_xyz) &
                                - ene_ij / exp_kR * sin(k_Rij) * dot_product(dkdA, Rij) &
                                + ene_ij / phi * dphi%dk_1 * dkdAk_proj
#ifdef IS_BAD_OPENMP_IFORT
                        end do ! i_xyz loop
#endif
                    end do
                end if
                if (grad%dC6) then
                    res%dE%dC6(i_atom) = res%dE%dC6(i_atom) + ene_ij / C6_ij * dC6%dC6i
                    res%dE%dC6(j_atom) = res%dE%dC6(j_atom) + ene_ij / C6_ij * dC6%dC6j
                end if
                if (grad%dalpha) then
                    res%dE%dalpha(i_atom) = res%dE%dalpha(i_atom) + ene_ij / C6_ij * dC6%da0i
                    res%dE%dalpha(j_atom) = res%dE%dalpha(j_atom) + ene_ij / C6_ij * dC6%da0j
                end if
            end do each_atom_pair
        end do each_atom
    end do each_recip_vec
    do i_atom = 1, geom%siz()
#ifdef WITH_MPI
        if (modulo(i_atom, geom%mpi_size) /= geom%mpi_rank) cycle
#endif
        res%energy = res%energy + geom%gamm**6 / 12 * C6(i_atom)  ! self energy
        if (grad%dC6) then
            res%dE%dC6(i_atom) = res%dE%dC6(i_atom) + geom%gamm**6 / 12
        end if
    end do
end subroutine

real(dp) function combine_C6( &
        C6_i, C6_j, alpha_0_i, alpha_0_j, dC6, grad) result(C6_ij)
    real(dp), intent(in) :: C6_i, C6_j, alpha_0_i, alpha_0_j
    type(grad_scalar_t), intent(out) :: dC6
    type(grad_request_t), intent(in) :: grad

    C6_ij = 2 * C6_i * C6_j &
        / (alpha_0_j / alpha_0_i * C6_i + alpha_0_i / alpha_0_j * C6_j)
    if (grad%dC6) then
        dC6%dC6i = C6_ij**2 * alpha_0_i / (2 * C6_i**2 * alpha_0_j)
        dC6%dC6j = C6_ij**2 * alpha_0_j / (2 * C6_j**2 * alpha_0_i)
    end if
    if (grad%dalpha) then
        dC6%da0i = C6_ij * (1 / alpha_0_i - C6_ij / (C6_i * alpha_0_j))
        dC6%da0j = C6_ij * (1 / alpha_0_j - C6_ij / (C6_j * alpha_0_i))
    end if
end function

real(dp) function disp_real(r, gamm, dphi, grad) result(phi)
    real(dp), intent(in) :: r, gamm
    type(grad_scalar_t), intent(out) :: dphi
    type(grad_request_t), intent(in) :: grad

    real(dp) :: gamm_r

    gamm_r = gamm * r
    phi = (2 + 2 * gamm_r**2 + gamm_r**4) * gamm**6 &
        / (2 * exp(gamm_r**2) * gamm_r**6)
    if (grad%dcoords) then
        dphi%dr_1 = -gamm**7 * (6 + 6 * gamm_r**2 + 3 * gamm_r**4 + gamm_r**6) &
            / (exp(gamm_r**2) * gamm_r**7)
    end if
end function

real(dp) function disp_rec(k, gamm, dphi, grad) result(phi)
    real(dp), intent(in) :: k, gamm
    type(grad_scalar_t), intent(out) :: dphi
    type(grad_request_t), intent(in) :: grad

    real(dp) :: k_gamm

    k_gamm = k / gamm
    phi = pi**1.5d0 * gamm**3 / 12 * ( &
        (-2 * (-2 + k_gamm**2)) / exp(k_gamm**2 / 4) &
        + k_gamm**3 * sqrt(pi) * erfc(k_gamm / 2) &
    )
    if (grad%dlattice) then
        dphi%dk_1 = k_gamm * pi**1.5d0 * gamm**2 / 4 * ( &
            -2 / exp(k_gamm**2 / 4) + k_gamm * sqrt(pi) * erfc(k_gamm / 2) &
        )
    end if
end function

end module
