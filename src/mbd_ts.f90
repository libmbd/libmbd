! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_ts
!! Obtaining TS energies.

use mbd_constants
use mbd_utils, only: shift_idx, tostr
use mbd_damping, only: damping_t, damping_fermi
use mbd_geom, only: geom_t

implicit none

private
public :: ts_energy

contains

function ts_energy(geom, alpha_0, C6, damp) result(ene)
    !! Get TS energy.
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha_0(:)
    real(dp), intent(in) :: C6(:)
    type(damping_t), intent(in) :: damp
    real(dp) :: ene

    real(dp) :: C6_ij, r(3), r_norm, R_vdw_ij, ene_shell, ene_pair, R_cell(3), &
        f_damp
    integer :: i_shell, i_cell, i_atom, j_atom, range_cell(3), idx_cell(3)
    real(dp), parameter :: shell_thickness = 10d0
    logical :: is_periodic

    is_periodic = allocated(geom%lattice)
    ene = 0d0
    i_shell = 0
    do
        i_shell = i_shell+1
        ene_shell = 0d0
        if (is_periodic) then
            range_cell = geom%supercell_circum(geom%lattice, i_shell*shell_thickness)
        else
            range_cell = [0, 0, 0]
        end if
        idx_cell = [0, 0, -1]
        do i_cell = 1, product(1+2*range_cell)
            call shift_idx(idx_cell, -range_cell, range_cell)
            if (is_periodic) then
                R_cell = matmul(geom%lattice, idx_cell)
            else
                R_cell = [0d0, 0d0, 0d0]
            end if
            do i_atom = 1, geom%siz()
                do j_atom = 1, i_atom
                    if (i_cell == 1) then
                        if (i_atom == j_atom) cycle
                    end if
                    r = geom%coords(:, i_atom)-geom%coords(:, j_atom)-R_cell
                    r_norm = sqrt(sum(r**2))
                    if (r_norm > geom%calc%param%ts_cutoff_radius) cycle
                    if (is_periodic) then
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
                            f_damp = damping_fermi(r, damp%ts_sr*R_vdw_ij, damp%ts_d)
                        case ("fermi2")
                            f_damp = damping_fermi(r, damp%ts_sr*R_vdw_ij, damp%ts_d)**2
                        case ("custom")
                            f_damp = damp%damping_custom(i_atom, j_atom)
                    end select
                    ene_pair = -C6_ij*f_damp/r_norm**6
                    if (i_atom == j_atom) then
                        ene_shell = ene_shell+ene_pair/2
                    else
                        ene_shell = ene_shell+ene_pair
                    endif
                end do ! j_atom
            end do ! i_atom
        end do ! i_cell
        ene = ene+ene_shell
        if (.not. is_periodic) exit
        if (i_shell > 1 .and. &
                abs(ene_shell) < geom%calc%param%ts_energy_accuracy) then
            call geom%calc%print( &
                "Periodic TS converged in " // trim(tostr(i_shell)) &
                // " shells, " // trim(tostr(i_shell*shell_thickness/ang)) &
                // " angstroms" &
            )
            exit
        endif
    end do ! i_shell
end function

elemental function combine_C6(C6_i, C6_j, alpha_0_i, alpha_0_j) result(C6_ij)
    real(dp), intent(in) :: C6_i, C6_j, alpha_0_i, alpha_0_j
    real(dp) :: C6_ij

    C6_ij = 2*C6_i*C6_j/(alpha_0_j/alpha_0_i*C6_i+alpha_0_i/alpha_0_j*C6_j)
end function

end module
