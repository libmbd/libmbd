! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LEGENDRE_PREC
#define LEGENDRE_PREC 8
#endif
#include "defaults.h"

module mbd_geom
!! Representing a molecule or a crystal unit cell.

use mbd_constants
use mbd_lapack, only: eigvals, inverse
use mbd_utils, only: shift_idx, atom_index_t, quad_pt_t, exception_t, tostr, clock_t
#ifdef WITH_SCALAPACK
use mbd_blacs, only: blacs_desc_t, blacs_grid_t
#endif
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

private
public :: supercell_circum, get_freq_grid

type, public :: param_t
    !! Calculation-wide paramters.
    real(dp) :: ts_energy_accuracy = TS_ENERGY_ACCURACY
    real(dp) :: ts_cutoff_radius = 50d0*ang
    real(dp) :: dipole_cutoff = 400d0*ang  ! used only when Ewald is off
    real(dp) :: ewald_real_cutoff_scaling = 1d0
    real(dp) :: ewald_rec_cutoff_scaling = 1d0
    real(dp) :: k_grid_shift = K_GRID_SHIFT
    logical :: ewald_on = .true.
    logical :: zero_negative_eigvals = .false.
    logical :: rpa_rescale_eigs = .false.
    integer :: rpa_order_max = 10
    integer :: n_freq = N_FREQUENCY_GRID
end type

type, public :: geom_t
    !! Represents a molecule or a crystal unit cell.
    !!
    !! The documented variables should be set before calling the initializer.
    real(dp), allocatable :: coords(:, :)
        !! (\(3\times N\), a.u.) Atomic coordinates.
    real(dp), allocatable :: lattice(:, :)
        !! (\(3\times 3\), a.u.) Lattice vectors in columns, unallocated if not
        !! periodic.
    integer, allocatable :: k_grid(:)
        !! Number of \(k\)-points along reciprocal axes.
    real(dp), allocatable :: custom_k_pts(:, :)
        !! Custom \(k\)-point grid.
    character(len=10) :: parallel_mode = 'auto'
        !! Type of parallelization:
        !!
        !! - `atoms`: distribute matrices over all MPI tasks using ScaLAPACK,
        !! solve eigenproblems sequentialy.
        !! - `k_points`: parallelize over k-points (each MPI task solves entire
        !! eigenproblems for its k-points)
    logical :: get_eigs = .false.
        !! Whether to keep MBD eigenvalues
    logical :: get_modes = .false.
        !! Whether to calculate MBD eigenvectors
    logical :: do_rpa = .false.
        !! Whether to calculate MBD energy by frequency integration
    logical :: get_rpa_orders = .false.
        !! Whether to calculate RPA orders
#ifdef WITH_MPI
    integer :: mpi_comm = MPI_COMM_WORLD
        !! MPI communicator
#endif
    ! The following components are set by the initializer and should be
    ! considered read-only
    type(clock_t) :: timer
    type(exception_t) :: exc
    logical :: muted = .false.
    type(quad_pt_t), allocatable :: freq(:)
    real(dp) :: gamm = 0d0
    real(dp) :: real_space_cutoff
    real(dp) :: rec_space_cutoff
    type(param_t) :: param
    type(atom_index_t) :: idx
#ifdef WITH_SCALAPACK
    ! TODO makes these two private (see use in mbd_methods, mbd_dipole)
    type(blacs_desc_t) :: blacs
    type(blacs_grid_t) :: blacs_grid
#endif
#ifdef WITH_MPI
    integer :: mpi_size = -1
    integer :: mpi_rank = -1
#endif
    contains
    procedure :: init => geom_init
    procedure :: destroy => geom_destroy
    procedure :: siz => geom_siz
    procedure :: has_exc => geom_has_exc
    procedure :: clock => geom_clock
end type

contains

subroutine geom_init(this)
    class(geom_t), intent(inout) :: this

    integer :: i_atom, n_kpts
    real(dp) :: volume
#ifdef WITH_MPI
    integer :: ierr
#endif

    associate (n => this%param%n_freq)
        allocate (this%freq(0:n))
        call get_freq_grid(n, this%freq(1:n)%val, this%freq(1:n)%weight)
    end associate
    this%freq(0)%val = 0d0
    this%freq(0)%weight = 0d0
    call this%timer%init(100)
    if (allocated(this%lattice)) then
        volume = abs(dble(product(eigvals(this%lattice))))
        if (this%param%ewald_on) then
            this%gamm = 2.5d0/volume**(1d0/3)
            this%real_space_cutoff = 6d0/this%gamm*this%param%ewald_real_cutoff_scaling
            this%rec_space_cutoff = 10d0*this%gamm*this%param%ewald_rec_cutoff_scaling
        else
            this%real_space_cutoff = this%param%dipole_cutoff
        end if
    end if
    if (this%parallel_mode == 'auto') then
        this%parallel_mode = 'atoms'
        n_kpts = -1
        if (allocated(this%custom_k_pts)) then
            n_kpts = size(this%custom_k_pts, 2)
        else if (allocated(this%k_grid)) then
            n_kpts = product(this%k_grid)
        end if
        if (allocated(this%lattice) .and. n_kpts > 0) then
            if (this%siz()**2 < n_kpts) then
                this%parallel_mode = 'k_points'
            end if
        end if
    end if
#ifdef WITH_MPI
    call MPI_COMM_SIZE(this%mpi_comm, this%mpi_size, ierr)
    call MPI_COMM_RANK(this%mpi_comm, this%mpi_rank, ierr)
#endif
#ifdef WITH_SCALAPACK
    this%idx%parallel = this%parallel_mode == 'atoms' .and. this%siz() > 1
    if (this%idx%parallel) then
#   ifdef WITH_MPI
        call this%blacs_grid%init(this%mpi_comm)
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
    deallocate (this%freq)
    deallocate (this%timer%timestamps)
    deallocate (this%timer%counts)
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

    has_exc = this%exc%code /= 0
end function

function supercell_circum(lattice, radius) result(sc)
    real(dp), intent(in) :: lattice(3, 3)
    real(dp), intent(in) :: radius
    integer :: sc(3)

    real(dp) :: ruc(3, 3), layer_sep(3)
    integer :: i

    ruc = 2*pi*inverse(transpose(lattice))
    forall (i = 1:3) &
        layer_sep(i) = sum(lattice(:, i)*ruc(:, i)/sqrt(sum(ruc(:, i)**2)))
    sc = ceiling(radius/layer_sep+0.5d0)
end function

subroutine geom_clock(this, id)
    class(geom_t), intent(inout) :: this
    integer, intent(in) :: id

    call this%timer%clock(id)
end subroutine

subroutine get_freq_grid(n, x, w, L)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n), w(n)
    real(dp), intent(in), optional :: L

    real(dp) :: L_

    if (present(L)) then
        L_ = L
    else
        L_ = 0.6d0
    end if
    call gauss_legendre(n, x, w)
    w = 2*L_/(1-x)**2*w
    x = L_*(1+x)/(1-x)
    w = w(n:1:-1)
    x = x(n:1:-1)
end subroutine

subroutine gauss_legendre(n, r, w)
    integer, intent(in) :: n
    real(dp), intent(out) :: r(n), w(n)

    integer, parameter :: q = LEGENDRE_PREC
    integer, parameter :: n_iter = 1000
    real(q) :: x, f, df, dx
    integer :: k, iter, i
    real(q) :: Pk(0:n), Pk1(0:n-1), Pk2(0:n-2)

    if (n == 1) then
        r(1) = 0d0
        w(1) = 2d0
        return
    end if
    Pk2(0) = 1._q  ! k = 0
    Pk1(0:1) = [0._q, 1._q]  ! k = 1
    do k = 2, n
        Pk(0:k) = ((2*k-1) * &
            [0._q, Pk1(0:k-1)]-(k-1)*[Pk2(0:k-2), 0._q, 0._q])/k
        if (k < n) then
            Pk2(0:k-1) = Pk1(0:k-1)
            Pk1(0:k) = Pk(0:k)
        end if
    end do
    ! now Pk contains k-th Legendre polynomial
    do i = 1, n
        x = cos(pi*(i-0.25_q)/(n+0.5_q))
        do iter = 1, n_iter
            df = 0._q
            f = Pk(n)
            do k = n-1, 0, -1
                df = f + x*df
                f = Pk(k) + x*f
            end do
            dx = f/df
            x = x-dx
            if (abs(dx) < 10*epsilon(dx)) exit
        end do
        r(i) = dble(x)
        w(i) = dble(2/((1-x**2)*df**2))
    end do
end subroutine

end module
