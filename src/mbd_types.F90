! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_types

use mbd_common, only: dp, findval
use mbd_parallel, only: mbd_blacs, mbd_blacs_grid, all_reduce

implicit none

private
public :: mat3n3n, mat33, scalar, contract_cross_33

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
    procedure :: mult_rows => mat3n3n_mult_rows
    procedure :: mult_cols_3n => mat3n3n_mult_cols_3n
    procedure :: mult_col => mat3n3n_mult_col
    procedure :: contract_n_transp => mat3n3n_contract_n_transp
    procedure :: contract_n33diag_cols => mat3n3n_contract_n33diag_cols
    procedure :: contract_n33_rows => mat3n3n_contract_n33_rows
    procedure :: copy_from => mat3n3n_copy_from
    procedure :: move_from => mat3n3n_move_from
    procedure :: init_from => mat3n3n_init_from
    procedure :: alloc_from => mat3n3n_alloc_from
end type

type :: mat33
    real(dp) :: val(3, 3)
    ! explicit derivative, [abc] ~ dval_{ab}/dR_c
    real(dp), allocatable :: dr(:, :, :)
    real(dp), allocatable :: dvdw(:, :)
    real(dp), allocatable :: dsigma(:, :)
end type

type :: scalar
    real(dp) :: val
    real(dp), allocatable :: dr(:)  ! explicit derivative
    real(dp), allocatable :: dvdw
end type

contains

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

subroutine mat3n3n_init(this, blacs)
    class(mat3n3n), intent(out) :: this
    type(mbd_blacs), intent(in) :: blacs

    this%blacs = blacs
end subroutine

subroutine mat3n3n_init_from(this, other)
    class(mat3n3n), intent(out) :: this
    type(mat3n3n), intent(in) :: other

    this%blacs = other%blacs
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

subroutine mat3n3n_alloc_from(this, other)
    class(mat3n3n), intent(out) :: this
    type(mat3n3n), intent(in) :: other

    integer :: n1, n2

    n1 = other%siz(1)
    n2 = other%siz(2)
    if (allocated(other%re)) then
        allocate (this%re(n1, n2))
    else
        allocate (this%cplx(n1, n2))
    end if
    this%blacs = other%blacs
end subroutine

subroutine mat3n3n_add(this, other)
    class(mat3n3n), intent(inout) :: this
    class(mat3n3n), intent(in) :: other

    if (allocated(this%re) .and. allocated(other%re)) then
        this%re = this%re + other%re
    else
        stop 1
    end if
end subroutine

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

subroutine mat3n3n_mult_rows(this, b)
    class(mat3n3n), intent(inout) :: this
    real(dp), intent(in) :: b(:)

    integer :: my_i_atom

    if (allocated(this%re)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            associate ( &
                    i_atom => this%blacs%i_atom(my_i_atom), &
                    this_sub => this%re(3*(my_i_atom-1)+1:, :) &
            )
                this_sub(:3, :) = this_sub(:3, :)*b(i_atom)
            end associate
        end do
    end if
    if (allocated(this%cplx)) then
        do my_i_atom = 1, size(this%blacs%i_atom)
            associate ( &
                    i_atom => this%blacs%i_atom(my_i_atom), &
                    this_sub => this%cplx(3*(my_i_atom-1)+1:, :) &
            )
                this_sub(:3, :) = this_sub(:3, :)*b(i_atom)
            end associate
        end do
    end if
end subroutine

subroutine mat3n3n_mult_cols_3n(this, b)
    class(mat3n3n), intent(inout) :: this
    real(dp), intent(in) :: b(:)

    integer :: my_j_atom, i

    if (allocated(this%re)) then
        do my_j_atom = 1, size(this%blacs%j_atom)
            associate ( &
                    b_sub => b(3*(this%blacs%j_atom(my_j_atom)-1)+1:), &
                    this_sub => this%re(:, 3*(my_j_atom-1)+1:) &
            )
                forall (i = 1:3) this_sub(:, i) = this_sub(:, i)*b_sub(i)
            end associate
        end do
    end if
    if (allocated(this%cplx)) then
        do my_j_atom = 1, size(this%blacs%j_atom)
            associate ( &
                    b_sub => b(3*(this%blacs%j_atom(my_j_atom)-1)+1:), &
                    this_sub => this%cplx(:, 3*(my_j_atom-1)+1:) &
            )
                forall (i = 1:3) this_sub(:, i) = this_sub(:, i)*b_sub(i)
            end associate
        end do
    end if
end subroutine

subroutine mat3n3n_mult_col(this, idx, a)
    class(mat3n3n), intent(inout) :: this
    integer, intent(in) :: idx
    real(dp), intent(in) :: a(:)

    integer :: my_i_atom, my_j_atom

    if (allocated(this%re)) then
        do my_j_atom = 1, size(this%blacs%j_atom)
            if (this%blacs%j_atom(my_j_atom) /= idx) cycle
            do my_i_atom = 1, size(this%blacs%i_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        this_sub => this%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    this_sub(:3, :3) = this_sub(:3, :3)*a(i_atom)
                end associate
            end do
        end do
    end if
    if (allocated(this%cplx)) then
        do my_j_atom = 1, size(this%blacs%j_atom)
            if (this%blacs%j_atom(my_j_atom) /= idx) cycle
            do my_i_atom = 1, size(this%blacs%i_atom)
                associate ( &
                        i_atom => this%blacs%i_atom(my_i_atom), &
                        this_sub => this%cplx(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
                )
                    this_sub(:3, :3) = this_sub(:3, :3)*a(i_atom)
                end associate
            end do
        end do
    end if
end subroutine

subroutine mat3n3n_contract_n_transp(this, dir, res)
    class(mat3n3n), intent(in) :: this
    character(len=*), intent(in) :: dir
    real(dp), intent(out), target :: res(:, :)

    integer :: my_i_atom, my_j_atom
    real(dp), pointer :: res_sub(:, :)

    res(:, :) = 0d0
    do my_i_atom = 1, size(this%blacs%i_atom)
        do my_j_atom = 1, size(this%blacs%j_atom)
            select case (dir(1:1))
            case ('R')
                res_sub => res(:, 3*(this%blacs%i_atom(my_i_atom)-1)+1:)
            case ('C')
                res_sub => res(3*(this%blacs%j_atom(my_j_atom)-1)+1:, :)
            end select
            associate ( &
                    this_sub => this%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                res_sub(:3, :3) = res_sub(:3, :3) + transpose(this_sub(:3, :3))
            end associate
        end do
    end do
    call all_reduce(res, this%blacs)
end subroutine

function contract_cross_33(k_atom, A, A_prime, B, B_prime) result(res)
    integer, intent(in) :: k_atom
    type(mat3n3n), intent(in) :: A, B
    real(dp), intent(in) :: A_prime(:, :), B_prime(:, :)
    real(dp) :: res(A%blacs%n_atoms)

    integer :: my_i_atom, my_j_atom, i_atom, j_atom

    res(:) = 0d0
    my_i_atom = findval(A%blacs%i_atom, k_atom)
    if (my_i_atom > 0) then
        do my_j_atom = 1, size(A%blacs%j_atom)
            j_atom = A%blacs%j_atom(my_j_atom)
            associate ( &
                    A_sub => A%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:), &
                    A_prime_sub => A_prime(:, 3*(j_atom-1)+1:) &
            )
                res(j_atom) = -1d0/3*sum(A_sub(:3, :3)*A_prime_sub(:, :3))
            end associate
        end do
    end if
    my_j_atom = findval(A%blacs%j_atom, k_atom)
    if (my_j_atom > 0) then
        do my_i_atom = 1, size(A%blacs%i_atom)
            i_atom = A%blacs%i_atom(my_i_atom)
            associate ( &
                    B_sub => B%re(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:), &
                    B_prime_sub => B_prime(3*(i_atom-1)+1:, :) &
            )
                res(i_atom) = res(i_atom) + &
                    (-1d0/3)*sum(B_prime_sub(:3, :)*B_sub(:3, :3))
            end associate
        end do
    end if
    call all_reduce(res, A%blacs)
end function

function mat3n3n_contract_n33diag_cols(A) result(res)
    class(mat3n3n), intent(in) :: A
    real(dp) :: res(A%blacs%n_atoms)

    integer :: i_xyz, my_j_atom, j_atom

    res(:) = 0d0
    do my_j_atom = 1, size(A%blacs%j_atom)
        j_atom = A%blacs%j_atom(my_j_atom)
        do i_xyz = 1, 3
            res(j_atom) = res(j_atom) + &
                sum(A%re(i_xyz::3, 3*(my_j_atom-1)+i_xyz))
        end do
    end do
    res = res/3
    call all_reduce(res, A%blacs)
end function

function mat3n3n_contract_n33_rows(A) result(res)
    class(mat3n3n), intent(in) :: A
    real(dp) :: res(A%blacs%n_atoms)

    integer :: my_i_atom, i_atom

    res(:) = 0d0
    do my_i_atom = 1, size(A%blacs%i_atom)
        i_atom = A%blacs%i_atom(my_i_atom)
        associate (A_sub => A%re(3*(my_i_atom-1)+1:, :))
            res(i_atom) = res(i_atom) + sum(A_sub(:3, :))
        end associate
    end do
    call all_reduce(res, A%blacs)
end function

end module
