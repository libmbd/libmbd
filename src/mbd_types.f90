! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_types

use mbd_common, only: dp
use mbd_parallel, only: mbd_blacs

implicit none

private
public :: mat3n3n, mat33, vecn, scalar
public :: operator(.cprod.), add_diag, symmetrize, mult_small, multed_small, &
    cross_self_add, cross_self_prod

type :: mat3n3n
    real(dp), allocatable :: re(:, :)
    complex(dp), allocatable :: cplx(:, :)
    real(dp), allocatable :: re_dr(:, :, :)
    real(dp), allocatable :: re_dvdw(:, :)
    real(dp), allocatable :: re_dsigma(:, :)
    type(mbd_blacs) :: blacs
    contains
    procedure  :: siz => mat3n3n_siz
    procedure :: init => mat3n3n_init
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

subroutine mat3n3n_init(this, n_atoms)
    class(mat3n3n), intent(inout) :: this
    integer, intent(in) :: n_atoms

    integer :: i

    this%blacs%i_atom = [(i, i = 1, n_atoms)]
    this%blacs%j_atom = this%blacs%i_atom
    this%blacs%n_atoms = n_atoms
end subroutine

integer function vecn_siz(this)
    class(vecn), intent(in) :: this

    if (allocated(this%val)) then
        vecn_siz = size(this%val)
    else
        vecn_siz = 0
    end if
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

    call add_diag_vec_(A, [(d, i = 1, A%blacs%n_atoms)])
end subroutine

subroutine add_diag_vec_(A, d)
    type(mat3n3n), intent(inout) :: A
    real(dp), intent(in) :: d(:)

    integer :: my_i_atom, my_j_atom, i

    if (allocated(A%re)) then
        do my_i_atom = 1, size(A%blacs%i_atom)
            do my_j_atom = 1, size(A%blacs%j_atom)
                associate ( &
                        i_atom => A%blacs%i_atom(my_i_atom), &
                        j_atom => A%blacs%j_atom(my_j_atom), &
                        A_diag => A%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    if (i_atom /= j_atom) cycle
                    do i = 1, 3
                        A_diag(i, i) = A_diag(i, i) + d(i_atom)
                    end do
                end associate
            end do
        end do
    end if
    if (allocated(A%cplx)) then
        do my_i_atom = 1, size(A%blacs%i_atom)
            do my_j_atom = 1, size(A%blacs%j_atom)
                associate ( &
                        i_atom => A%blacs%i_atom(my_i_atom), &
                        j_atom => A%blacs%j_atom(my_j_atom), &
                        A_diag => A%cplx(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    if (i_atom /= j_atom) cycle
                    do i = 1, 3
                        A_diag(i, i) = A_diag(i, i) + d(i_atom)
                    end do
                end associate
            end do
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

function symmetrize(A)
    real(dp), intent(in) :: A(:, :)
    real(dp) :: symmetrize(size(A, 1), size(A, 1))

    symmetrize = A + transpose(A)
end function

end module
