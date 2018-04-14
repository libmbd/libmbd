! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_types

use mbd_common, only: dp

implicit none

private
public :: mat3n3n, mat33, vecn, scalar
public :: eye, diag, operator(.cprod.), add_diag, repeatn, &
    symmetrize, mult_small, multed_small, operator(.cadd.), cross_self_add, &
    cross_self_prod

type :: mat3n3n
    real(dp), allocatable :: re(:, :)
    complex(dp), allocatable :: cplx(:, :)
    real(dp), allocatable :: re_dr(:, :, :)
    real(dp), allocatable :: re_dvdw(:, :)
    real(dp), allocatable :: re_dsigma(:, :)
    contains
    procedure  :: siz => mat3n3n_siz
end type

type :: mat33
    real(dp) :: val(3, 3)
    ! explicit derivative, [abc] ~ dval_{ab}/dR_c
    real(dp), allocatable :: dr(:, :, :)
    real(dp), allocatable :: dvdw(:, :)
    real(dp), allocatable :: dsigma(:, :)
end type

type :: vecn
    real(dp), allocatable :: val(:)
    real(dp), allocatable :: dr(:, :, :)
    contains
    procedure  :: siz => vecn_siz
end type

interface vecn
    module procedure vecn_constructor_no_dr
end interface

type :: scalar
    real(dp) :: val
    real(dp), allocatable :: dr(:)  ! explicit derivative
    real(dp), allocatable :: dvdw
end type

interface operator(.cprod.)
    module procedure cart_prod_
end interface

interface operator(.cadd.)
    module procedure cart_add_
end interface

interface diag
    module procedure get_diag_
    module procedure get_diag_cmplx_
    module procedure make_diag_
end interface

interface add_diag
    module procedure add_diag_scalar_
    module procedure add_diag_vec_
end interface

contains

type(vecn) function vecn_constructor_no_dr(x) result(vec)
    real(dp), intent(in) :: x(:)

    vec%val = x
end function

integer function mat3n3n_siz(this, ndim)
    class(mat3n3n), intent(in) :: this
    integer, intent(in) :: ndim

    if (allocated(this%re)) then
        mat3n3n_siz = size(this%re, ndim)
    elseif (allocated(this%cplx)) then
        mat3n3n_siz = size(this%cplx, ndim)
    else
        mat3n3n_siz = 0
    end if
end function

integer function vecn_siz(this)
    class(vecn), intent(in) :: this

    if (allocated(this%val)) then
        vecn_siz = size(this%val)
    else
        vecn_siz = 0
    end if
end function

function eye(n) result(A)
    integer, intent(in) :: n
    real(dp) :: A(n, n)

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:n) A(i, i) = 1.d0
end function

function cart_prod_(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))

    integer :: i, j

    do i = 1, size(a)
        do j = 1, size(b)
            c(i, j) = a(i)*b(j)
        end do
    end do
end function

function cart_add_(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))

    integer :: i, j

    do i = 1, size(a)
        do j = 1, size(b)
            c(i, j) = a(i)+b(j)
        end do
    end do
end function

function cross_self_prod(a) result(c)
    real(dp), intent(in) :: a(:)
    real(dp) :: c(size(a), size(a))

    c = cart_prod_(a, a)
end function

function cross_self_add(a) result(c)
    real(dp), intent(in) :: a(:)
    real(dp) :: c(size(a), size(a))

    c = cart_add_(a, a)
end function


subroutine add_diag_scalar_(A, d)
    type(mat3n3n), intent(inout) :: A
    real(dp), intent(in) :: d

    integer :: i

    call add_diag_vec_(A, [(d, i = 1, A%siz(1))])
end subroutine

subroutine add_diag_vec_(A, d)
    type(mat3n3n), intent(inout) :: A
    real(dp), intent(in) :: d(:)

    integer :: i

    if (allocated(A%re)) then
        do i = 1, size(d)
            A%re(i, i) = A%re(i, i) + d(i)
        end do
    end if
    if (allocated(A%cplx)) then
        do i = 1, size(d)
            A%cplx(i, i) = A%cplx(i, i) + d(i)
        end do
    end if
end subroutine

subroutine mult_small(A, B)
    real(dp), intent(inout) :: A(:, :)
    real(dp), intent(in) :: B(:, :)

    integer :: i, i3, j, j3

    forall (i = 1:size(B, 1), i3 = 1:3, j = 1:size(B, 1), j3 = 1:3)
        A((i-1)*3+i3, (j-1)*3+j3) = B(i, j)*A((i-1)*3+i3, (j-1)*3+j3)
    end forall
end subroutine

function multed_small(A, B)
    real(dp), intent(in) :: A(:, :)
    real(dp), intent(in) :: B(:, :)
    real(dp) :: multed_small(size(A, 1), size(A, 1))

    multed_small = A
    call mult_small(multed_small, B)
end function


function repeatn(x, n)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n
    real(dp) :: repeatn(n*size(x))

    integer :: i, j

    repeatn = [([(x(i), j = 1, n)], i = 1, size(x))]
end function


function get_diag_(A) result(d)
    real(dp), intent(in) :: A(:, :)
    real(dp) :: d(size(A, 1))

    integer :: i

    forall (i = 1:size(A, 1)) d(i) = A(i, i)
end function


function get_diag_cmplx_(A) result(d)
    complex(dp), intent(in) :: A(:, :)
    complex(dp) :: d(size(A, 1))

    integer :: i

    forall (i = 1:size(A, 1)) d(i) = A(i, i)
end function


function make_diag_(d) result(A)
    real(dp), intent(in) :: d(:)
    real(dp) :: A(size(d), size(d))

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:size(d)) A(i, i) = d(i)
end function

function symmetrize(A)
    real(dp), intent(in) :: A(:, :)
    real(dp) :: symmetrize(size(A, 1), size(A, 1))

    symmetrize = A + transpose(A)
end function

end module
