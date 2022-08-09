! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_density

use mbd_constants
use mbd_geom, only: geom_t
use mbd_linalg, only: diag
use mbd_lapack, only: eigvalsh, inverse

implicit none

contains

function eval_mbd_nonint_density(geom, pts, charges, masses, omegas) result(rho)
    type(geom_t), intent(in) :: geom
    real(dp), intent(in) :: pts(:, :), charges(:), masses(:), omegas(:)
    real(dp) :: rho(size(pts, 2))

    integer :: i_pt, i_atom, n_atoms
    real(dp), dimension(:), allocatable :: pre, kernel, rsq

    pre = charges * (masses * omegas / pi)**(3.d0 / 2)
    kernel = masses * omegas
    n_atoms = geom%siz()
    rho(:) = 0.d0
    allocate (rsq(n_atoms))
    do i_pt = 1, size(pts, 2)
        do concurrent(i_atom=1:n_atoms)
            rsq(i_atom) = sum((pts(:, i_pt) - geom%coords(:, i_atom))**2)
        end do
        rho(i_pt) = sum(pre * exp(-kernel * rsq))
    end do
end function

function eval_mbd_int_density(geom, pts, charges, masses, omegas, modes) result(rho)
    type(geom_t), intent(in) :: geom
    real(dp), intent(in) :: pts(:, :), charges(:), masses(:), omegas(:), modes(:, :)
    real(dp) :: rho(size(pts, 2))

    integer :: i_pt, i_atom, n_atoms, i, i_xyz, j_xyz, self(3)
    real(dp) :: rdiffsq(3, 3), rdiff(3)
    integer, allocatable :: other(:)
    real(dp), allocatable :: pre(:), factor(:), omegas_p(:, :), kernel(:, :, :)

    omegas_p = matmul(matmul(modes, diag(omegas)), transpose(modes))
    n_atoms = geom%siz()
    allocate (kernel(3, 3, n_atoms), source=0.d0)
    allocate (pre(n_atoms), source=0.d0)
    do i_atom = 1, n_atoms
        self = (/(3 * (i_atom - 1) + i, i=1, 3)/)
        other = (/(i, i=1, 3 * (i_atom - 1)), (i, i=3 * i_atom + 1, 3 * n_atoms)/)
        kernel(:, :, i_atom) = masses(i_atom) &
            * (omegas_p(self, self) &
                - matmul(matmul(omegas_p(self, other), inverse(omegas_p(other, other))), &
                    omegas_p(other, self)))
        pre(i_atom) = charges(i_atom) * (masses(i_atom) / pi)**(3.d0 / 2) &
            * sqrt(product(omegas) / product(eigvalsh(omegas_p(other, other))))
    end do
    rho(:) = 0.d0
    allocate (factor(n_atoms))
    do i_pt = 1, size(pts, 2)
        do i_atom = 1, n_atoms
            rdiff = pts(:, i_pt) - geom%coords(:, i_atom)
            do concurrent(i_xyz=1:3, j_xyz=1:3)
                rdiffsq(i_xyz, j_xyz) = rdiff(i_xyz) * rdiff(j_xyz)
            end do
            factor(i_atom) = sum(kernel(:, :, i_atom) * rdiffsq(:, :))
        end do
        rho(i_pt) = sum(pre * exp(-factor))
    end do
end function

end module
