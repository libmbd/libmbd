! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_parallel_impl

use mbd_parallel

implicit none

private
public :: init_blacs_grid, init_blacs

contains

subroutine init_blacs_grid(this)
    type(mbd_blacs_grid), intent(out) :: this

    this%nprows = 1
    this%npcols = 1
    this%my_prow = 1
    this%my_pcol = 1
end subroutine

subroutine init_blacs(this, n_atoms)
    type(mbd_blacs), intent(out) :: this
    integer, intent(in) :: n_atoms

    integer :: i

    this%i_atom = [(i, i = 1, n_atoms)]
    this%j_atom = this%i_atom
    this%n_atoms = n_atoms
end subroutine

end module

