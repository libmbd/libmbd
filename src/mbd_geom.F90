! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_geom
!! Representing a molecule or a crystal unit cell.

use mbd_constants
use mbd_calc, only: calc_t
use mbd_lapack, only: inverse
use mbd_utils, only: shift_idx, atom_index_t
#ifdef WITH_SCALAPACK
use mbd_blacs, only: blacs_desc_t, blacs_grid_t
#endif
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

private

type, public :: geom_t
    !! Represents a molecule or a crystal unit cell.
    type(calc_t), pointer :: calc
    real(dp), allocatable :: coords(:, :)  ! 3 by n_atoms
    logical :: vacuum_axis(3) = [.false., .false., .false.]
    real(dp), allocatable :: lattice(:, :)  ! vectors in columns
    integer :: k_grid(3)  ! TODO make allocatable
    real(dp), allocatable :: k_pts(:, :)
    real(dp), allocatable :: gamma_ew
        ! TODO create ewald_t type

    integer :: supercell(3)
    character(len=10) :: parallel_mode = 'auto'
        !! Type of parallelization:
        !!
        !! - `atoms`: distribute matrices over all MPI tasks using ScaLAPACK,
        !! solve eigenproblems sequentialy.
        !! - `k_points`: parallelize over k-points (each MPI task solves entire
        !! eigenproblems for its k-points)
    type(atom_index_t) :: idx
#ifdef WITH_SCALAPACK
    ! TODO makes these two private (see use in mbd_methods, mbd_dipole)
    type(blacs_desc_t) :: blacs
    type(blacs_grid_t) :: blacs_grid
#endif
#ifdef WITH_MPI
    integer :: comm = MPI_COMM_WORLD
#endif
    contains
    procedure :: init => geom_init
    procedure :: destroy => geom_destroy
    procedure :: siz => geom_siz
    procedure :: has_exc => geom_has_exc
    procedure :: supercell_circum => geom_supercell_circum
    procedure :: ensure_k_pts => geom_ensure_k_pts
    procedure :: clock => geom_clock
end type

contains

subroutine geom_init(this, calc)
    class(geom_t), intent(inout) :: this
    type(calc_t), target, intent(in) :: calc

    integer :: i_atom

    this%calc => calc
    if (this%parallel_mode == 'auto') this%parallel_mode = 'atoms'
        ! TODO put some logic here
#ifdef WITH_SCALAPACK
    this%idx%parallel = this%parallel_mode == 'atoms'
    if (this%idx%parallel) then
#   ifdef WITH_MPI
        call this%blacs_grid%init(this%comm)
#   else
        call this%blacs_grid%init()
#   endif
        call this%blacs%init(this%siz(), this%blacs_grid)
        this%idx%i_atom = this%blacs%i_atom
        this%idx%j_atom = this%blacs%j_atom
    else
        this%idx%i_atom = [(i_atom, i_atom = 1, this%siz())]
        this%idx%j_atom = this%idx%i_atom
    end if
#else
    this%idx%i_atom = [(i_atom, i_atom = 1, this%siz())]
    this%idx%j_atom = this%idx%i_atom
#endif
    this%idx%n_atoms = this%siz()
end subroutine

subroutine geom_destroy(this)
    class(geom_t), intent(inout) :: this

#ifdef WITH_SCALAPACK
    if (this%idx%parallel) call this%blacs_grid%destroy()
#endif
end subroutine

integer function geom_siz(this) result(siz)
    class(geom_t), intent(in) :: this

    if (allocated(this%coords)) then
        siz = size(this%coords, 2)
    else
        siz = 0
    end if
end function

logical function geom_has_exc(this) result(has_exc)
    class(geom_t), intent(in) :: this

    has_exc = this%calc%exc%code /= 0
end function

function geom_supercell_circum(this, uc, radius) result(sc)
    class(geom_t), intent(in) :: this
    real(dp), intent(in) :: uc(3, 3), radius
    integer :: sc(3)

    real(dp) :: ruc(3, 3), layer_sep(3)
    integer :: i

    ruc = 2*pi*inverse(transpose(uc))
    forall (i = 1:3) &
        layer_sep(i) = sum(uc(:, i)*ruc(:, i)/sqrt(sum(ruc(:, i)**2)))
    sc = ceiling(radius/layer_sep+0.5d0)
    where (this%vacuum_axis) sc = 0
end function

subroutine geom_ensure_k_pts(this)
    class(geom_t), intent(inout) :: this

    if (allocated(this%k_pts)) return
    if (.not. allocated(this%lattice)) return
    this%k_pts = make_k_pts(this%k_grid, this%lattice, this%calc%param%k_grid_shift)
end subroutine

function make_k_pts(k_grid, lattice, shift) result(k_pts)
    integer, intent(in) :: k_grid(3)
    real(dp), intent(in) :: lattice(3, 3)
    real(dp), intent(in), optional :: shift
    real(dp) :: k_pts(3, product(k_grid))

    integer :: n_kpt(3), i_kpt
    real(dp) :: n_kpt_shifted(3)

    n_kpt = [0, 0, -1]
    do i_kpt = 1, product(k_grid)
        call shift_idx(n_kpt, [0, 0, 0], k_grid-1)
        n_kpt_shifted = dble(n_kpt)
        if (present(shift)) n_kpt_shifted = n_kpt_shifted+shift
        where (2*n_kpt_shifted > k_grid) n_kpt_shifted = n_kpt_shifted-dble(k_grid)
        k_pts(:, i_kpt) = n_kpt_shifted/k_grid
    end do
    k_pts = matmul(2*pi*inverse(transpose(lattice)), k_pts)
end function

subroutine geom_clock(this, id)
    class(geom_t), intent(inout) :: this
    integer, intent(in) :: id

    call this%calc%clock%clock(id)
end subroutine

end module
