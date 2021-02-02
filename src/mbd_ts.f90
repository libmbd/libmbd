! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_ts
!! Obtaining TS energies.

use mbd_constants
use mbd_utils, only: shift_idx, tostr, result_t, diff3
use mbd_damping, only: damping_t, damping_fermi
use mbd_geom, only: geom_t, supercell_circum
use mbd_gradients, only: grad_request_t
use mbd_lapack, only: eigvals, inverse

implicit none

private
public :: get_ts_energy, get_ts_energy_num_grad

contains

type(result_t) function get_ts_energy_num_grad(geom, alpha_0, C6, damp, grad) result(res)
    !! Get TS energy and numerical gradients.
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    type(grad_request_t), intent(in) :: grad

    integer :: i_atom, i_xyz, i_step, i_latt
    real(dp) :: delta
    real(dp), allocatable :: ene_diffed(:)
    real(dp), allocatable :: coords_orig(:, :), lattice_orig(:, :)

    res%energy = get_ts_energy(geom, alpha_0, C6, damp)
    if (.not. grad%any()) return
    delta = geom%param%ts_num_grad_delta
    allocate (ene_diffed(-1:1))
    if (grad%dcoords) then
        allocate (res%dE%dcoords(geom%siz(), 3))
        do i_atom = 1, geom%siz()
            do i_xyz = 1, 3
                do i_step = -1, 1
                    if (i_step == 0) cycle
                    coords_orig = geom%coords
                    geom%coords(i_xyz, i_atom) = geom%coords(i_xyz, i_atom) + i_step * delta
                    ene_diffed(i_step) = get_ts_energy(geom, alpha_0, C6, damp)
                    geom%coords = coords_orig
                end do
                res%dE%dcoords(i_atom, i_xyz) = diff3(ene_diffed, delta)
            end do
        end do
    end if
    if (grad%dlattice) then
        allocate (res%dE%dlattice(3, 3))
        do i_latt = 1, 3
            do i_xyz = 1, 3
                do i_step = -1, 1
                    if (i_step == 0) cycle
                    lattice_orig = geom%lattice
                    geom%lattice(i_xyz, i_latt) = geom%lattice(i_xyz, i_latt) + i_step * delta
                    ene_diffed(i_step) = get_ts_energy(geom, alpha_0, C6, damp)
                    geom%lattice = lattice_orig
                end do
                res%dE%dlattice(i_latt, i_xyz) = diff3(ene_diffed, delta)
            end do
        end do
    end if
end function

function get_ts_energy(geom, alpha_0, C6, damp) result(ene)
    !! Get TS energy.
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    real(dp) :: ene

    real(dp) :: C6_ij, Rnij(3), Rnij_norm, R_vdw_ij, ene_ij, Rn(3), &
        f_damp
    integer :: i_cell, i_atom, j_atom, range_n(3), n(3)
    logical :: is_periodic, do_ewald

    do_ewald = .false.
    is_periodic = allocated(geom%lattice)
    if (is_periodic) then
        do_ewald = geom%gamm > 0d0
        range_n = supercell_circum(geom%lattice, geom%real_space_cutoff)
    else
        range_n(:) = 0
    end if
    ene = 0d0
    n = [0, 0, -1]
    each_cell: do i_cell = 1, product(1 + 2 * range_n)
        call shift_idx(n, -range_n, range_n)
        if (is_periodic) then
            Rn = matmul(geom%lattice, n)
        else
            Rn(:) = 0d0
        end if
        each_atom: do i_atom = 1, geom%siz()
            each_atom_pair: do j_atom = 1, i_atom
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                Rnij = geom%coords(:, i_atom) - geom%coords(:, j_atom) - Rn
                Rnij_norm = sqrt(sum(Rnij**2))
                if (is_periodic .and. Rnij_norm > geom%real_space_cutoff) cycle
                C6_ij = combine_C6( &
                    C6(i_atom), C6(j_atom), alpha_0(i_atom), alpha_0(j_atom) &
                )
                if (allocated(damp%r_vdw)) then
                    R_vdw_ij = damp%r_vdw(i_atom) + damp%r_vdw(j_atom)
                end if
                select case (damp%version)
                    case ("fermi")
                        f_damp = damping_fermi(Rnij, damp%ts_sr * R_vdw_ij, damp%ts_d)
                    case ("fermi2")
                        f_damp = damping_fermi(Rnij, damp%ts_sr * R_vdw_ij, damp%ts_d)**2
                    case ("custom")
                        f_damp = damp%damping_custom(i_atom, j_atom)
                end select
                ene_ij = -C6_ij * f_damp / Rnij_norm**6
                if (do_ewald) then
                    ene_ij = ene_ij &
                        + C6_ij * (disp_real(Rnij_norm, geom%gamm) + 1d0 / Rnij_norm**6)
                end if
                if (i_atom == j_atom) ene_ij = ene_ij / 2
                ene = ene + ene_ij
            end do each_atom_pair
        end do each_atom
    end do each_cell
    if (do_ewald) call add_ewald_ts_parts(geom, alpha_0, C6, ene)
end function

subroutine add_ewald_ts_parts(geom, alpha_0, C6, ene)
    type(geom_t), intent(in) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    real(dp), intent(inout) :: ene

    real(dp) :: rec_latt(3, 3), volume, Rij(3), k(3), &
        k_norm, k_Rij, latt_inv(3, 3), C6_ij, exp_kR, ene_ij
    integer :: i_atom, j_atom, m(3), i_m, range_m(3)

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
            each_atom_pair: do j_atom = 1, i_atom
                C6_ij = combine_C6( &
                    C6(i_atom), C6(j_atom), alpha_0(i_atom), alpha_0(j_atom) &
                )
                Rij = geom%coords(:, i_atom) - geom%coords(:, j_atom)
                k_Rij = dot_product(k, Rij)
                exp_kR = cos(k_Rij)
                ene_ij = C6_ij * disp_rec(k_norm, geom%gamm) / volume * exp_kR
                if (i_atom == j_atom) ene_ij = ene_ij / 2
                ene = ene + ene_ij
            end do each_atom_pair
        end do each_atom
    end do each_recip_vec
    ene = ene + geom%gamm**6 / 12 * sum(C6)  ! self energy
end subroutine

elemental function combine_C6(C6_i, C6_j, alpha_0_i, alpha_0_j) result(C6_ij)
    real(dp), intent(in) :: C6_i, C6_j, alpha_0_i, alpha_0_j
    real(dp) :: C6_ij

    C6_ij = 2 * C6_i * C6_j / (alpha_0_j / alpha_0_i * C6_i + alpha_0_i / alpha_0_j * C6_j)
end function

elemental real(dp) function disp_real(r, gamm)
    real(dp), intent(in) :: r, gamm

    disp_real = -( &
        1 / r**6 + gamm**2 / r**4 + gamm**4 / (2 * r**2) &
    ) * exp(-((gamm * r)**2))
end function

elemental real(dp) function disp_rec(k, gamm)
    real(dp), intent(in) :: k, gamm

    disp_rec = -(pi**(3d0 / 2) / 12) * ( &
        sqrt(pi) * k**3 * erfc(k / (2 * gamm)) &
        + (4 * gamm**3 - 2 * gamm * k**2) &
        * exp(-(k**2) / (4 * gamm**2)) &
    )
end function

end module
