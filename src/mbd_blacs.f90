! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_blacs

use mbd_constants

implicit none

private
public :: blacs_all_reduce

integer, parameter :: DLEN = 9

type, public :: blacs_grid_t
    integer :: ctx
    integer :: comm
    integer :: nprows
    integer :: npcols
    integer :: my_prow
    integer :: my_pcol
contains
    procedure :: init => blacs_grid_init
    procedure :: destroy => blacs_grid_destroy
end type

type, public :: blacs_desc_t
    integer, allocatable :: i_atom(:)
    integer, allocatable :: j_atom(:)
    integer :: n_atoms
    integer :: desc(DLEN)
    integer :: ctx
    integer :: blocksize
    integer :: comm = -1
contains
    procedure :: init => blacs_desc_init
end type

interface blacs_all_reduce
    module procedure all_reduce_real_scalar
    module procedure all_reduce_complex_scalar
    module procedure all_reduce_real_1d
    module procedure all_reduce_real_2d
    module procedure all_reduce_complex_1d
    module procedure all_reduce_complex_2d
end interface


interface

    subroutine blacs_pinfo(id, nproc)
        integer, intent(out) :: id, nproc
    end subroutine blacs_pinfo

    subroutine blacs_gridinit(ictxt, order, nprow, npcol)
        integer, intent(inout) :: ictxt
        character, intent(in) :: order
        integer, intent(in) :: nprow, npcol
    end subroutine blacs_gridinit

    subroutine blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)
        integer, intent(in) :: ictxt,nprow, npcol
        integer, intent(out) :: myrow, mycol
    end subroutine blacs_gridinfo

    subroutine blacs_gridexit(ictxt)
        integer, intent(in) :: ictxt
    end subroutine blacs_gridexit

    integer function numroc(nn, nb, iproc, isrcproc, nprocs)
        integer, intent(in) :: nn, nb, iproc, isrcproc, nprocs
    end function numroc

    subroutine descinit(desc, mm, nn, mb, nb, irsrc, icsrc, ictxt, lld, info)
        import :: DLEN
        integer, intent(out) :: desc(DLEN)
        integer, intent(in) :: mm, nn, mb, nb, irsrc, icsrc, ictxt, lld
        integer, intent(out) :: info
    end subroutine descinit

    subroutine blacs_get(ictxt, what, val)
        integer, intent(in) :: ictxt, what
        integer, intent(out) :: val
    end subroutine blacs_get

    subroutine dgsum2d(ictxt, scope, top, mm, nn, aa, lda, rdest, cdest)
        import :: dp
        integer, intent(in) :: ictxt
        character, intent(in) :: scope, top
        integer, intent(in) :: mm, nn
        integer, intent(in) :: lda
        real(dp), intent(inout) :: aa(lda,*)
        integer, intent(in) :: rdest, cdest
    end subroutine dgsum2d

    subroutine zgsum2d(ictxt, scope, top, mm, nn, aa, lda, rdest, cdest)
        import :: dp
        integer, intent(in) :: ictxt
        character, intent(in) :: scope, top
        integer, intent(in) :: mm, nn
        integer, intent(in) :: lda
        complex(dp), intent(inout) :: aa(lda,*)
        integer, intent(in) :: rdest, cdest
    end subroutine zgsum2d

end interface


contains

subroutine blacs_grid_init(this, comm)
    class(blacs_grid_t), intent(inout) :: this
    integer, intent(in), optional :: comm

    integer :: my_task, n_tasks, nprows

    call BLACS_PINFO(my_task, n_tasks)
    do nprows = int(sqrt(dble(n_tasks))), 1, -1
        if (mod(n_tasks, nprows) == 0) exit
    enddo
    this%nprows = nprows
    this%npcols = n_tasks/this%nprows
    if (present(comm)) then
        this%ctx = comm
        this%comm = comm
    else
        call BLACS_GET(0, 0, this%ctx)
    end if
    call BLACS_GRIDINIT(this%ctx, 'R', this%nprows, this%npcols)
    call BLACS_GRIDINFO( &
        this%ctx, this%nprows, this%npcols, this%my_prow, this%my_pcol &
    )
end subroutine

! TODO this should be made a destructor once support for gfortran 4.9 is dropped
subroutine blacs_grid_destroy(this)
    class(blacs_grid_t), intent(inout) :: this

    call BLACS_GRIDEXIT(this%ctx)
end subroutine

subroutine blacs_desc_init(this, n_atoms, grid, max_atoms_per_block)
    class(blacs_desc_t), intent(out) :: this
    type(blacs_grid_t), intent(in) :: grid
    integer, intent(in) :: n_atoms
    integer, intent(in) :: max_atoms_per_block

    integer :: my_nratoms, my_ncatoms, ierr, atoms_per_block, n_proc

    this%comm = grid%comm
    this%ctx = grid%ctx
    this%n_atoms = n_atoms
    n_proc = max(grid%nprows, grid%npcols)
    if (n_proc == 1) return
    atoms_per_block = (n_atoms - 1) / (n_proc - 1)
    if (atoms_per_block == 0) return
    atoms_per_block = min(atoms_per_block, max_atoms_per_block)
    my_nratoms = NUMROC(n_atoms, atoms_per_block, grid%my_prow, 0, grid%nprows)
    my_ncatoms = NUMROC(n_atoms, atoms_per_block, grid%my_pcol, 0, grid%npcols)
    this%blocksize = 3*atoms_per_block
    call DESCINIT( &
        this%desc, 3*n_atoms, 3*n_atoms, this%blocksize, this%blocksize, 0, 0, &
        grid%ctx, 3*my_nratoms, ierr &
    )
    this%i_atom = idx_map( &
        grid%my_prow, grid%nprows, n_atoms, atoms_per_block, my_nratoms &
    )
    this%j_atom = idx_map( &
        grid%my_pcol, grid%npcols, n_atoms, atoms_per_block, my_ncatoms &
    )
end subroutine

function idx_map(my_task, n_tasks, n, blocksize, nidx)
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

subroutine all_reduce_real_scalar(x, blacs)
    real(dp), intent(inout) :: x
    type(blacs_desc_t), intent(in) :: blacs

    real(dp) :: xtemp(1)

    xtemp(1) = x
    call DGSUM2D(blacs%ctx, 'A', ' ', 1, 1, xtemp, 1, -1, -1)
    x = xtemp(1)

end subroutine

subroutine all_reduce_complex_scalar(x, blacs)
    complex(dp), intent(inout) :: x
    type(blacs_desc_t), intent(in) :: blacs

    complex(dp) :: xtemp(1)

    xtemp(1) = x
    call ZGSUM2D(blacs%ctx, 'A', ' ', 1, 1, xtemp, 1, -1, -1)
    x = xtemp(1)

end subroutine

subroutine all_reduce_real_1d(A, blacs)
    real(dp), intent(inout) :: A(:)
    type(blacs_desc_t), intent(in) :: blacs

    call DGSUM2D(blacs%ctx, 'A', ' ', size(A), 1, A, size(A), -1, -1)
end subroutine

subroutine all_reduce_real_2d(A, blacs)
    real(dp), contiguous, target, intent(inout) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs

    real(dp), pointer :: ptrA(:)

    ptrA(1 : size(A)) => A
    call DGSUM2D( &
        blacs%ctx, 'A', ' ', size(A, 1), size(A, 2), ptrA, size(A, 1), -1, -1 &
    )
end subroutine

subroutine all_reduce_complex_1d(A, blacs)
    complex(dp), intent(inout) :: A(:)
    type(blacs_desc_t), intent(in) :: blacs

    call ZGSUM2D(blacs%ctx, 'A', ' ', size(A), 1, A, size(A), -1, -1)
end subroutine

subroutine all_reduce_complex_2d(A, blacs)
    complex(dp), contiguous, target, intent(inout) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs

    complex(dp), pointer :: ptrA(:)

    ptrA(1 : size(A)) => A
    call ZGSUM2D( &
        blacs%ctx, 'A', ' ', size(A, 1), size(A, 2), ptrA, size(A, 1), -1, -1 &
    )
end subroutine

end module
