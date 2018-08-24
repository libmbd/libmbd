! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_parallel

use mbd_common, only: dp

implicit none

private
public :: mbd_blacs, mbd_blacs_grid, all_reduce

type :: mbd_blacs_grid
    integer :: ctx = -1
    integer :: nprows
    integer :: npcols
    integer :: my_prow
    integer :: my_pcol
contains
    procedure :: init => mbd_blacs_grid_init
    final :: mbd_blacs_grid_destroy
end type

type :: mbd_blacs
    integer, allocatable :: i_atom(:)
    integer, allocatable :: j_atom(:)
    integer :: n_atoms
    integer :: desc(9)
    integer :: ctx = -1
contains
    procedure :: init => mbd_blacs_init
    procedure :: parallel => mbd_blacs_parallel
end type

interface all_reduce
    module procedure all_reduce_re_1d_
    module procedure all_reduce_re_2d_
end interface

#ifdef WITH_SCALAPACK
external :: BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDINFO, &
    BLACS_GRIDEXIT, NUMROC, DESCINIT, DGSUM2D
integer :: NUMROC
#endif

contains


subroutine mbd_blacs_grid_init(this)
    class(mbd_blacs_grid), intent(inout) :: this

#ifdef WITH_SCALAPACK
    integer :: my_task, n_tasks, nprows

    call BLACS_PINFO(my_task, n_tasks)
    do nprows = int(sqrt(dble(n_tasks))), 1, -1
        if (mod(n_tasks, nprows) == 0) exit
    enddo
    this%nprows = nprows
    this%npcols = n_tasks/this%nprows
    call BLACS_GET(-1, 0, this%ctx)
    call BLACS_GRIDINIT(this%ctx, 'R', this%nprows, this%npcols)
    call BLACS_GRIDINFO( &
        this%ctx, this%nprows, this%npcols, this%my_prow, this%my_pcol &
    )
#else
    this%nprows = 1
    this%npcols = 1
    this%my_prow = 1
    this%my_pcol = 1
#endif
end subroutine

subroutine mbd_blacs_grid_destroy(this)
    type(mbd_blacs_grid) :: this
#ifdef WITH_SCALAPACK

    if (this%ctx /= -1) call BLACS_GRIDEXIT(this%ctx)
#endif
end subroutine

subroutine mbd_blacs_init(this, n_atoms, grid, parallel)
    class(mbd_blacs), intent(out) :: this
    type(mbd_blacs_grid), intent(in) :: grid
    integer, intent(in) :: n_atoms
    logical, intent(in) :: parallel

    integer :: blocksize, my_nratoms, my_ncatoms, ierr, i

    this%ctx = grid%ctx
    this%n_atoms = n_atoms
    if (parallel) then
#ifdef WITH_SCALAPACK
        blocksize = 3
        my_nratoms = NUMROC(n_atoms, blocksize/3, grid%my_prow, 0, grid%nprows)
        my_ncatoms = NUMROC(n_atoms, blocksize/3, grid%my_pcol, 0, grid%npcols)
        call DESCINIT( &
            this%desc, 3*n_atoms, 3*n_atoms, blocksize, blocksize, 0, 0, &
            grid%ctx, 3*my_nratoms, ierr &
        )
        this%i_atom = get_idx_map( &
            grid%my_prow, grid%nprows, n_atoms, blocksize/3, my_nratoms &
        )
        this%j_atom = get_idx_map( &
            grid%my_pcol, grid%npcols, n_atoms, blocksize/3, my_ncatoms &
        )
#else
        this%desc(1) = -1
        this%i_atom = [(i, i = 1, n_atoms)]
        this%j_atom = this%i_atom
#endif
    else
        this%desc(1) = -1
        this%i_atom = [(i, i = 1, n_atoms)]
        this%j_atom = this%i_atom
    end if
end subroutine

logical function mbd_blacs_parallel(this)
    class(mbd_blacs), intent(in) :: this

    mbd_blacs_parallel = this%n_atoms /= size(this%i_atom) .or. &
        this%n_atoms /= size(this%j_atom)
end function

function get_idx_map(my_task, n_tasks, n, blocksize, nidx) result(idx_map)
    integer, intent(in) :: my_task, n_tasks, n, blocksize, nidx
    integer :: idx_map(nidx)

    integer :: i, i_block, n_in_block, my_i

    i_block = 0
    n_in_block = 0
    my_i = 1
    do i = 1, n
        if (mod(i_block, n_tasks) == my_task) then
            idx_map(my_i) = i
            my_i = my_i + 1
        end if
        n_in_block = n_in_block + 1
        if (n_in_block == blocksize) then
            n_in_block = 0
            i_block = i_block + 1
        end if
    end do
end function

subroutine all_reduce_re_1d_(A, blacs)
    real(dp), intent(inout) :: A(:)
    type(mbd_blacs), intent(in) :: blacs
#ifdef WITH_SCALAPACK

    call DGSUM2D(blacs%ctx, 'A', ' ', size(A), 1, A, size(A), -1, -1)
#endif
end subroutine

subroutine all_reduce_re_2d_(A, blacs)
    real(dp), intent(inout) :: A(:, :)
    type(mbd_blacs), intent(in) :: blacs
#ifdef WITH_SCALAPACK

    call DGSUM2D( &
        blacs%ctx, 'A', ' ', size(A, 1), size(A, 2), A, size(A, 1), -1, -1 &
    )
#endif
end subroutine

end module
