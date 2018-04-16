! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_parallel

implicit none

private
public :: mbd_blacs, mbd_blacs_grid

type :: mbd_blacs_grid
    integer :: ctx
    integer :: nprows
    integer :: npcols
    integer :: my_prow
    integer :: my_pcol
contains
    procedure :: init => mbd_blacs_grid_init
    procedure :: destroy => mbd_blacs_grid_destroy
end type

type :: mbd_blacs
    integer, allocatable :: i_atom(:)
    integer, allocatable :: j_atom(:)
    integer :: n_atoms
    integer :: desc(9)
    type(mbd_blacs_grid) :: grid
contains
    procedure :: init => mbd_blacs_init
end type

contains

subroutine mbd_blacs_grid_init(this)
    class(mbd_blacs_grid), intent(out) :: this

    this%nprows = 1
    this%npcols = 1
    this%my_prow = 1
    this%my_pcol = 1
end subroutine

subroutine mbd_blacs_grid_destroy(this)
    class(mbd_blacs_grid), intent(inout) :: this
end subroutine

subroutine mbd_blacs_init(this, n_atoms, grid)
    class(mbd_blacs), intent(out) :: this
    integer, intent(in) :: n_atoms
    type(mbd_blacs_grid), intent(in) :: grid

    integer :: i

    this%grid = grid
    this%i_atom = [(i, i = 1, n_atoms)]
    this%j_atom = this%i_atom
    this%n_atoms = n_atoms
end subroutine

end module
