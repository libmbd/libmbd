! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LEGENDRE_PREC
#define LEGENDRE_PREC 8
#endif
module mbd_system_type

use mbd_constants
use mbd_common, only: mbd_exc, printer, tostr
use mbd_lapack, only: inverse
use mbd_matrix_type, only: mbd_index
#ifdef WITH_SCALAPACK
use mbd_blacs, only: mbd_blacs_desc, mbd_blacs_grid
#endif
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

#ifndef MODULE_UNIT_TESTS
private
public :: mbd_system, mbd_calc, ang, clock_rate
#endif
real(dp), parameter :: ang = 1.8897259886d0
integer, parameter :: n_timestamps = 100

type :: mbd_param
    real(dp) :: ts_energy_accuracy = TS_ENERGY_ACCURACY
    real(dp) :: ts_cutoff_radius = 50d0*ang
    real(dp) :: dipole_low_dim_cutoff = 100d0*ang
    real(dp) :: dipole_cutoff = 400d0*ang  ! used only when Ewald is off
    real(dp) :: ewald_real_cutoff_scaling = 1d0
    real(dp) :: ewald_rec_cutoff_scaling = 1d0
    real(dp) :: k_grid_shift = K_GRID_SHIFT
    logical :: ewald_on = .true.
    logical :: zero_negative_eigs = .false.
    integer :: rpa_order_max = 10
    integer :: n_frequency_grid = N_FREQUENCY_GRID
end type

type :: mbd_timing
    logical :: measure_time = .true.
    integer :: timestamps(n_timestamps), ts_counts(n_timestamps)
    integer :: ts_cnt, ts_rate, ts_cnt_max, ts_aid
end type mbd_timing

type :: mbd_info
    character(len=120) :: ewald_alpha = '', ewald_rsum = '', ts_conv = '', &
        ewald_cutoff = '', ewald_recsum = '', freq_n = '', freq_error = '', &
        neg_eigvals = ''
    contains
    procedure :: print => mbd_info_print
end type mbd_info

type :: mbd_calc
    type(mbd_param) :: param
    type(mbd_timing) :: tm
    real(dp), allocatable :: omega_grid(:)
    real(dp), allocatable :: omega_grid_w(:)
    type(mbd_exc) :: exc
    type(mbd_info) :: info
    contains
    procedure :: init_grid => mbd_calc_init_grid
end type mbd_calc

type :: mbd_system
    type(mbd_calc), pointer :: calc
    real(dp), allocatable :: coords(:, :)  ! 3 by n_atoms
    logical :: vacuum_axis(3) = [.false., .false., .false.]
    real(dp), allocatable :: lattice(:, :)  ! vectors in columns
    integer :: k_grid(3)
    integer :: supercell(3)
    logical :: do_rpa = .false.
    logical :: get_eigs = .false.
    logical :: get_modes = .false.
    logical :: get_rpa_orders = .false.
    !> Type of parallelization: `"atoms"` or `"k_points"`.
    !>
    !> - `"atoms"`: distribute matrices over all MPI tasks using ScaLAPACK, solve
    !> eigenproblems sequentialy.
    !> - `"k_points"`: parallelize over k-points (each MPI task solves entire
    !> eigenproblems for its k-points)
    character(len=10) :: parallel_mode = 'auto'
    type(mbd_index) :: idx
#ifdef WITH_SCALAPACK
    type(mbd_blacs_desc) :: blacs
    type(mbd_blacs_grid) :: blacs_grid
#endif
#ifdef WITH_MPI
    integer :: comm = MPI_COMM_WORLD
#endif
    contains
    procedure :: init => mbd_system_init
    procedure :: destroy => mbd_system_destroy
    procedure :: siz => mbd_system_siz
    procedure :: has_exc => mbd_system_has_exc
    procedure :: supercell_circum => mbd_system_supercell_circum
    procedure :: clock => mbd_system_clock
end type mbd_system

contains

subroutine mbd_system_init(this, calc)
    class(mbd_system), intent(inout) :: this
    type(mbd_calc), target, intent(in) :: calc

    integer :: i_atom

    this%calc => calc
    ! TODO put some logic here
    if (this%parallel_mode == 'auto') this%parallel_mode = 'atoms'
#ifdef WITH_SCALAPACK
    this%idx%parallel = this%parallel_mode == 'atoms'
    if (this%idx%parallel) then
#ifdef WITH_MPI
        call this%blacs_grid%init(this%comm)
#else
        call this%blacs_grid%init()
#endif
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

subroutine mbd_system_destroy(this)
    class(mbd_system), intent(inout) :: this
#ifdef WITH_SCALAPACK

    if (this%idx%parallel) call this%blacs_grid%destroy()
#endif
end subroutine

integer function mbd_system_siz(this) result(siz)
    class(mbd_system), intent(in) :: this

    if (allocated(this%coords)) then
        siz = size(this%coords, 2)
    else
        siz = 0
    end if
end function

logical function mbd_system_has_exc(this) result(has_exc)
    class(mbd_system), intent(in) :: this

    has_exc = this%calc%exc%code /= 0
end function

function mbd_system_supercell_circum(this, uc, radius) result(sc)
    class(mbd_system), intent(in) :: this
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

subroutine mbd_system_clock(this, id, always)
    class(mbd_system), intent(inout) :: this
    integer, intent(in) :: id
    logical, intent(in), optional :: always

    associate (tm => this%calc%tm)
        if (tm%measure_time .or. present(always)) then
            call system_clock(tm%ts_cnt, tm%ts_rate, tm%ts_cnt_max)
            if (id > 0) then
                tm%timestamps(id) = tm%timestamps(id)-tm%ts_cnt
            else
                tm%ts_aid = abs(id)
                tm%timestamps(tm%ts_aid) = &
                    tm%timestamps(tm%ts_aid)+tm%ts_cnt
                tm%ts_counts(tm%ts_aid) = &
                    tm%ts_counts(tm%ts_aid)+1
            end if
        end if
    end associate
end subroutine

subroutine mbd_info_print(this, info)
    class(mbd_info), intent(in) :: this
    procedure(printer) :: info

    if (this%freq_n /= '') call info(this%freq_n)
    if (this%freq_error /= '') call info(this%freq_error)
    if (this%ewald_alpha /= '') call info(this%ewald_alpha)
    if (this%ewald_rsum /= '') call info(this%ewald_rsum)
    if (this%ewald_cutoff /= '') call info(this%ewald_cutoff)
    if (this%ewald_recsum /= '') call info(this%ewald_recsum)
    if (this%ts_conv /= '') call info(this%ts_conv)
    if (this%neg_eigvals /= '') call info(this%neg_eigvals)
end subroutine

subroutine mbd_calc_init_grid(this)
    class(mbd_calc), intent(inout) :: this

    integer :: n

    n = this%param%n_frequency_grid
    allocate (this%omega_grid(0:n))
    allocate (this%omega_grid_w(0:n))
    this%omega_grid(0) = 0d0
    this%omega_grid_w(0) = 0d0
    call get_omega_grid(n, 0.6d0, this%omega_grid(1:n), this%omega_grid_w(1:n))
    this%info%freq_n = &
        "Initialized a radial integration grid of " // trim(tostr(n)) // &
        " points."
end subroutine

subroutine get_omega_grid(n, L, x, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: L
    real(dp), intent(out) :: x(n), w(n)

    call gauss_legendre(n, x, w)
    w = 2*L/(1-x)**2*w
    x = L*(1+x)/(1-x)
    w = w(n:1:-1)
    x = x(n:1:-1)
end subroutine get_omega_grid

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

function clock_rate() result(rate)
    integer :: cnt, rate, cnt_max

    call system_clock(cnt, rate, cnt_max)
end function clock_rate

end module
