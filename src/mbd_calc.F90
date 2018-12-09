! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LEGENDRE_PREC
#define LEGENDRE_PREC 8
#endif
module mbd_calc

use mbd_constants
use mbd_utils, only: tostr, exception_t, clock_t

implicit none

private
public :: calc_t, get_freq_grid

type :: param_t
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

type :: info_t
    character(len=120) :: ewald_alpha = '', ewald_rsum = '', ts_conv = '', &
        ewald_cutoff = '', ewald_recsum = '', freq_n = '', freq_error = '', &
        neg_eigvals = ''
    contains
    procedure :: print => info_print
end type

type :: calc_t
    type(param_t) :: param
    type(clock_t) :: clock
    real(dp), allocatable :: omega_grid(:)
    real(dp), allocatable :: omega_grid_w(:)
    type(exception_t) :: exc
    type(info_t) :: info
    logical :: do_rpa = .false.
    logical :: get_eigs = .false.
    logical :: get_modes = .false.
    logical :: get_rpa_orders = .false.
    contains
    procedure :: init => calc_init
end type

abstract interface
    subroutine printer(msg)
        character(len=*), intent(in) :: msg
    end subroutine
end interface

contains

subroutine info_print(this, info)
    class(info_t), intent(in) :: this
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

subroutine calc_init(this)
    class(calc_t), intent(inout) :: this

    integer :: n

    call this%clock%init(100)
    n = this%param%n_frequency_grid
    allocate (this%omega_grid(0:n))
    allocate (this%omega_grid_w(0:n))
    this%omega_grid(0) = 0d0
    this%omega_grid_w(0) = 0d0
    call get_freq_grid(n, this%omega_grid(1:n), this%omega_grid_w(1:n))
    this%info%freq_n = &
        "Initialized a radial integration grid of " // trim(tostr(n)) // &
        " points."
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
end subroutine get_freq_grid

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
