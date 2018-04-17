! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_types

use mbd_common, only: dp
use mbd_parallel, only: mbd_blacs, mbd_blacs_grid

implicit none

private
public :: mat3n3n, mat33, vecn, scalar
public :: operator(.cprod.)

type :: mat3n3n
    real(dp), allocatable :: re(:, :)
    complex(dp), allocatable :: cplx(:, :)
    real(dp), allocatable :: re_dr(:, :, :)
    real(dp), allocatable :: re_dvdw(:, :)
    real(dp), allocatable :: re_dsigma(:, :)
    type(mbd_blacs) :: blacs
    contains
    procedure :: siz => mat3n3n_siz
    procedure :: init => mat3n3n_init
    procedure :: add => mat3n3n_add
    procedure :: add_diag => mat3n3n_add_diag
    procedure :: add_diag_scalar => mat3n3n_add_diag_scalar
    procedure :: mult_cross => mat3n3n_mult_cross
    procedure :: multed_cross => mat3n3n_multed_cross
    procedure :: mult_cross_add => mat3n3n_mult_cross_add
    procedure :: mult_dsigma => mat3n3n_mult_dsigma
    procedure :: copy_from => mat3n3n_copy_from
    procedure :: move_from => mat3n3n_move_from
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
    module procedure vecn_constructor_no_dr_
end interface

type :: scalar
    real(dp) :: val
    real(dp), allocatable :: dr(:)  ! explicit derivative
    real(dp), allocatable :: dvdw
end type

interface operator(.cprod.)
    module procedure cart_prod_
end interface

contains

type(vecn) function vecn_constructor_no_dr_(x) result(vec)
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

subroutine mat3n3n_init(this, n_atoms, blacs_grid)
    class(mat3n3n), intent(out) :: this
    integer, intent(in) :: n_atoms
    type(mbd_blacs_grid), intent(in) :: blacs_grid

    call this%blacs%init(n_atoms, blacs_grid)
end subroutine

subroutine mat3n3n_copy_from(this, other)
    class(mat3n3n), intent(out) :: this
    type(mat3n3n), intent(in) :: other

    if (allocated(other%re)) then
        this%re = other%re
    else
        this%cplx = other%cplx
    end if
    this%blacs = other%blacs
end subroutine

subroutine mat3n3n_move_from(this, other)
    class(mat3n3n), intent(out) :: this
    type(mat3n3n), intent(inout) :: other

    if (allocated(other%re)) then
        call move_alloc(other%re, this%re)
    else
        call move_alloc(other%cplx, this%cplx)
    end if
    this%blacs = other%blacs
end subroutine

integer function vecn_siz(this)
    class(vecn), intent(in) :: this

    if (allocated(this%val)) then
        vecn_siz = size(this%val)
    else
        vecn_siz = 0
    end if
end function

subroutine mat3n3n_add(this, other)
    class(mat3n3n), intent(inout) :: this
    class(mat3n3n), intent(in) :: other

    if (allocated(this%re) .and. allocated(other%re)) then
        this%re = this%re + other%re
    else
        stop 1
    end if
end subroutine

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

subroutine mat3n3n_add_diag_scalar(this, d)
    class(mat3n3n), intent(inout) :: this
    real(dp), intent(in) :: d

    integer :: i

    call mat3n3n_add_diag(this, [(d, i = 1, this%blacs%n_atoms)])
end subroutine

subroutine mat3n3n_add_diag(this, d)
    class(mat3n3n), intent(inout) :: this
    real(dp), intent(in) :: d(:)

    integer :: my_i_atom, my_j_atom, i

    if (allocated(this%re)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            do my_j_atom = 1, size(this%blacs%j_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        j_atom => this%blacs%j_atom(my_j_atom), &
                        this_diag => this%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    if (i_atom /= j_atom) cycle
                    do i = 1, 3
                        this_diag(i, i) = this_diag(i, i) + d(i_atom)
                    end do
                end associate
            end do
        end do
    end if
    if (allocated(this%cplx)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            do my_j_atom = 1, size(this%blacs%j_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        j_atom => this%blacs%j_atom(my_j_atom), &
                        this_diag => this%cplx(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    if (i_atom /= j_atom) cycle
                    do i = 1, 3
                        this_diag(i, i) = this_diag(i, i) + d(i_atom)
                    end do
                end associate
            end do
        end do
    end if
end subroutine

subroutine mat3n3n_mult_cross(this, b, c)
    class(mat3n3n), intent(inout) :: this
    real(dp), intent(in) :: b(:)
    real(dp), intent(in), optional :: c(:)

    integer :: my_i_atom, my_j_atom

    if (allocated(this%re)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            do my_j_atom = 1, size(this%blacs%j_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        j_atom => this%blacs%j_atom(my_j_atom), &
                        this_sub => this%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    if (present(c)) then
                        this_sub(:3, :3) = this_sub(:3, :3) * &
                            (b(i_atom)*c(j_atom)+c(i_atom)*b(j_atom))
                    else
                        this_sub(:3, :3) = this_sub(:3, :3)*b(i_atom)*b(j_atom)
                    end if
                end associate
            end do
        end do
    end if
    if (allocated(this%cplx)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            do my_j_atom = 1, size(this%blacs%j_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        j_atom => this%blacs%j_atom(my_j_atom), &
                        this_sub => this%cplx(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    if (present(c)) then
                        this_sub(:3, :3) = this_sub(:3, :3) * &
                            (b(i_atom)*c(j_atom)+c(i_atom)*b(j_atom))
                    else
                        this_sub(:3, :3) = this_sub(:3, :3)*b(i_atom)*b(j_atom)
                    end if
                end associate
            end do
        end do
    end if
end subroutine

subroutine mat3n3n_mult_cross_add(this, b)
    class(mat3n3n), intent(inout) :: this
    real(dp), intent(in) :: b(:)

    integer :: my_i_atom, my_j_atom

    if (allocated(this%re)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            do my_j_atom = 1, size(this%blacs%j_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        j_atom => this%blacs%j_atom(my_j_atom), &
                        this_sub => this%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    this_sub(:3, :3) = this_sub(:3, :3)*(b(i_atom)+b(j_atom))
                end associate
            end do
        end do
    end if
    if (allocated(this%cplx)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            do my_j_atom = 1, size(this%blacs%j_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        j_atom => this%blacs%j_atom(my_j_atom), &
                        this_sub => this%cplx(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    this_sub(:3, :3) = this_sub(:3, :3)*(b(i_atom)+b(j_atom))
                end associate
            end do
        end do
    end if
end subroutine

subroutine mat3n3n_mult_dsigma(this, sigma, dsigma)
    class(mat3n3n), intent(inout) :: this
    real(dp), intent(in) :: sigma(:), dsigma(:)

    integer :: my_i_atom, my_j_atom

    if (allocated(this%re)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            do my_j_atom = 1, size(this%blacs%j_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        j_atom => this%blacs%j_atom(my_j_atom), &
                        this_sub => this%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    this_sub(:3, :3) = this_sub(:3, :3) * &
                        (sigma(i_atom)*dsigma(i_atom)+sigma(j_atom)*dsigma(j_atom)) / &
                        sqrt(sigma(i_atom)**2+sigma(j_atom)**2)
                end associate
            end do
        end do
    end if
    if (allocated(this%cplx)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            do my_j_atom = 1, size(this%blacs%j_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        j_atom => this%blacs%j_atom(my_j_atom), &
                        this_sub => this%cplx(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    this_sub(:3, :3) = this_sub(:3, :3) * &
                        (sigma(i_atom)*dsigma(i_atom)+sigma(j_atom)*dsigma(j_atom)) / &
                        sqrt(sigma(i_atom)**2+sigma(j_atom)**2)
                end associate
            end do
        end do
    end if
end subroutine

type(mat3n3n) function mat3n3n_multed_cross(this, b, c) result(res)
    class(mat3n3n), intent(in) :: this
    real(dp), intent(in) :: b(:), c(:)

    ! expanded res = this assignment because of compiler bug in gfortran 4.9
    if (allocated(this%re)) res%re = this%re
    if (allocated(this%cplx)) res%cplx = this%cplx
    res%blacs = this%blacs
    call res%mult_cross(b, c)
end function

end module
