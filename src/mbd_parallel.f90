! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_parallel

implicit none

private
public :: mbd_blacs, mbd_blacs_grid

type :: mbd_blacs
  integer, allocatable :: i_atom(:)
  integer, allocatable :: j_atom(:)
  integer :: n_atoms
end type

type :: mbd_blacs_grid
    integer :: ctx
    integer :: nprows
    integer :: npcols
    integer :: my_prow
    integer :: my_pcol
end type

end module
