! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_gradients

use mbd_constants

implicit none

private
public :: grad_t, grad_matrix_re_t, grad_matrix_cplx_t, grad_request_t, &
    grad_scalar_t

type :: grad_t
    real(dp), allocatable :: dcoords(:, :)  ! n_atoms by 3
    real(dp), allocatable :: dalpha(:)
    real(dp), allocatable :: dalpha_dyn(:, :)  ! n_atoms by 0:n_freq
    real(dp), allocatable :: dC6(:)
    real(dp), allocatable :: dr_vdw(:)
    real(dp), allocatable :: domega(:)
    real(dp), allocatable :: dV(:)
    real(dp), allocatable :: dV_free(:)
    real(dp), allocatable :: dX_free(:)
end type

type :: grad_request_t
    logical :: dcoords = .false.
    logical :: dalpha = .false.
    logical :: dalpha_dyn = .false.
    logical :: dC6 = .false.
    logical :: dr_vdw = .false.
    logical :: domega = .false.
    logical :: dsigma = .false.
    logical :: dV = .false.
    logical :: dV_free = .false.
    logical :: dX_free = .false.
    contains
    procedure :: any => grad_request_any
end type

type :: grad_matrix_re_t
    real(dp), allocatable :: dr(:, :, :)
    real(dp), allocatable :: dvdw(:, :)
    real(dp), allocatable :: dsigma(:, :)
end type

type :: grad_matrix_cplx_t
    complex(dp), allocatable :: dr(:, :, :)
    complex(dp), allocatable :: dvdw(:, :)
    complex(dp), allocatable :: dsigma(:, :)
end type

type :: grad_scalar_t
    real(dp), allocatable :: dr(:)  ! explicit derivative
    real(dp), allocatable :: dvdw
end type

contains

logical function grad_request_any(this) result(any)
    class(grad_request_t), intent(in) :: this

    any = this%dcoords .or. this%dalpha .or. this%dC6 .or. &
        this%dr_vdw .or. this%domega
end function

end module
