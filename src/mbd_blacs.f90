! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_blacs

use mbd_constants

implicit none

private
public :: blacs_all_reduce, sys2blacs_handle, free_blacs_system_handle

type, public :: blacs_grid_t
    integer :: ctx = -1
    integer :: sys_ctx = -1
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
    integer :: desc(9)
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
    ! The following interfaces were written by hand based on
    ! https://www.netlib.org/blacs/BLACS/QRef.html and https://www.ibm.com/docs/
    subroutine BLACS_PINFO(MYPNUM, NPROCS)
        integer :: MYPNUM, NPROCS
    end
    subroutine BLACS_GRIDINIT(ICONTXT, ORDER, NPROW, NPCOL)
        integer :: ICONTXT, NPROW, NPCOL
        character :: ORDER
    end
    subroutine BLACS_GRIDINFO(ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL)
        integer :: ICONTXT, NPROW, NPCOL, MYPROW, MYPNUM
    end
    subroutine BLACS_GRIDEXIT(ICONTXT)
        integer :: ICONTXT
    end
    subroutine BLACS_GET(ICONTXT, WHAT, VAL)
        integer :: ICONTXT, WHAT, VAL
    end subroutine
    integer function SYS2BLACS_HANDLE(COMM)
        integer :: COMM
    end
    subroutine FREE_BLACS_SYSTEM_HANDLE(IHANDLE)
        integer :: IHANDLE
    end
    integer function NUMROC(n, nb, iproc, isrcproc, nprocs)
        integer :: n, nb, iproc, isrcproc, nprocs
    end
    subroutine DGSUM2D(ICONTXT, SCOPE, TOP, M, N, A, LDA, RDEST, CDEST)
        integer :: ICONTXT, RDEST, CDEST, M, N, LDA
        character :: SCOPE, TOP
        double precision :: A(*)
    end
    subroutine ZGSUM2D(ICONTXT, SCOPE, TOP, M, N, A, LDA, RDEST, CDEST)
        import :: dp
        integer :: ICONTXT, RDEST, CDEST, M, N, LDA
        character :: SCOPE, TOP
        complex(dp) :: A(*)
    end

    ! The following interfaces were taken straight from the ScaLAPACK codebase
    SUBROUTINE DESCINIT(DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO)
        INTEGER ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
        INTEGER DESC(*)
    END
end interface

contains

subroutine blacs_grid_init(this, n_tasks, ctx)
    class(blacs_grid_t), intent(inout) :: this
    integer, intent(in), optional :: n_tasks
    integer, intent(in), optional :: ctx

    integer :: my_task, n_tasks_, nprows

    if (present(n_tasks)) then
        n_tasks_ = n_tasks
    else
        call BLACS_PINFO(my_task, n_tasks_)
    end if
    do nprows = int(sqrt(dble(n_tasks_))), 1, -1
        if (mod(n_tasks_, nprows) == 0) exit
    end do
    this%nprows = nprows
    this%npcols = n_tasks_ / this%nprows
    if (present(ctx)) then
        this%ctx = ctx
        this%sys_ctx = ctx
    else
        call BLACS_GET(0, 0, this%ctx)
        this%sys_ctx = -1
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
    if (this%sys_ctx /= -1) then
        call FREE_BLACS_SYSTEM_HANDLE(this%sys_ctx)
        this%sys_ctx = -1
    end if
end subroutine

subroutine blacs_desc_init(this, n_atoms, grid, max_atoms_per_block)
    class(blacs_desc_t), intent(out) :: this
    type(blacs_grid_t), intent(in) :: grid
    integer, intent(in) :: n_atoms
    integer, intent(in) :: max_atoms_per_block

    integer :: my_nratoms, my_ncatoms, ierr, atoms_per_block, n_proc

    this%ctx = grid%ctx
    this%n_atoms = n_atoms
    n_proc = max(grid%nprows, grid%npcols)
    if (n_proc == 1) return
    atoms_per_block = (n_atoms - 1) / (n_proc - 1)
    if (atoms_per_block == 0) return
    atoms_per_block = min(atoms_per_block, max_atoms_per_block)
    my_nratoms = NUMROC(n_atoms, atoms_per_block, grid%my_prow, 0, grid%nprows)
    my_ncatoms = NUMROC(n_atoms, atoms_per_block, grid%my_pcol, 0, grid%npcols)
    this%blocksize = 3 * atoms_per_block
    call DESCINIT( &
        this%desc, 3 * n_atoms, 3 * n_atoms, this%blocksize, this%blocksize, 0, 0, &
        grid%ctx, 3 * my_nratoms, ierr &
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

    real(dp), pointer :: x_arr(:, :)

    x_arr(1, 1) = x
    call DGSUM2D(blacs%ctx, 'A', ' ', 1, 1, x_arr, 1, -1, -1)
    x = x_arr(1, 1)
end subroutine

subroutine all_reduce_complex_scalar(x, blacs)
    complex(dp), intent(inout) :: x
    type(blacs_desc_t), intent(in) :: blacs

    complex(dp) :: x_arr(1, 1)

    x_arr(1, 1) = x
    call ZGSUM2D(blacs%ctx, 'A', ' ', 1, 1, x_arr, 1, -1, -1)
    x = x_arr(1, 1)
end subroutine

subroutine all_reduce_real_1d(A, blacs)
    real(dp), target, intent(inout) :: A(:)
    type(blacs_desc_t), intent(in) :: blacs

    real(dp), pointer :: A_p(:, :)

    A_p(1:size(A), 1:1) => A
    call DGSUM2D(blacs%ctx, 'A', ' ', size(A), 1, A_p, size(A), -1, -1)
end subroutine

subroutine all_reduce_real_2d(A, blacs)
    real(dp), intent(inout) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs

    call DGSUM2D( &
        blacs%ctx, 'A', ' ', size(A, 1), size(A, 2), A, size(A, 1), -1, -1 &
    )
end subroutine

subroutine all_reduce_complex_1d(A, blacs)
    complex(dp), target, intent(inout) :: A(:)
    type(blacs_desc_t), intent(in) :: blacs

    complex(dp), pointer :: A_p(:, :)

    A_p(1:size(A), 1:1) => A
    call ZGSUM2D(blacs%ctx, 'A', ' ', size(A), 1, A_p, size(A), -1, -1)
end subroutine

subroutine all_reduce_complex_2d(A, blacs)
    complex(dp), intent(inout) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs

    call ZGSUM2D( &
        blacs%ctx, 'A', ' ', size(A, 1), size(A, 2), A, size(A, 1), -1, -1 &
    )
end subroutine

end module
