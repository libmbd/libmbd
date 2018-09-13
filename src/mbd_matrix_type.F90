! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MBD_INCLUDED
module mbd_matrix_type

use mbd_constants
use mbd_common, only: findval, exception => mbd_exc
use mbd_lapack, only: mmul, invh, invh, eigh, eigvals, eigvalsh
#ifdef WITH_SCALAPACK
use mbd_blacs, only: mbd_blacs_desc, mbd_blacs_grid, all_reduce
use mbd_scalapack, only: pmmul, pinvh, pinvh, peigh, peigvalsh
#endif

implicit none

private
public :: mbd_matrix_real, mbd_matrix_complex, contract_cross_33, mbd_index

type :: mbd_index
    integer, allocatable :: i_atom(:)
    integer, allocatable :: j_atom(:)
    integer :: n_atoms
#ifdef WITH_SCALAPACK
    logical :: parallel
#endif
end type

type :: mbd_matrix_real
    real(dp), allocatable :: val(:, :)
    type(mbd_index) :: idx
#ifdef WITH_SCALAPACK
    type(mbd_blacs_desc) :: blacs
#endif
    contains
    procedure :: siz => mbd_matrix_real_siz
    procedure :: init => mbd_matrix_real_init
    procedure :: add => mbd_matrix_real_add
    procedure :: add_diag => mbd_matrix_real_add_diag
    procedure :: add_diag_scalar => mbd_matrix_real_add_diag_scalar
    procedure :: mult_cross => mbd_matrix_real_mult_cross
    procedure :: mult_rows => mbd_matrix_real_mult_rows
    procedure :: mult_cols_3n => mbd_matrix_real_mult_cols_3n
    procedure :: mult_col => mbd_matrix_real_mult_col
    procedure :: mmul => mbd_matrix_real_mmul
    procedure :: invh => mbd_matrix_real_invh
    procedure :: eigh => mbd_matrix_real_eigh
    procedure :: eigvals => mbd_matrix_real_eigvals
    procedure :: eigvalsh => mbd_matrix_real_eigvalsh
    procedure :: contract_n_transp => mbd_matrix_real_contract_n_transp
    procedure :: contract_n33diag_cols => mbd_matrix_real_contract_n33diag_cols
    procedure :: contract_n33_rows => mbd_matrix_real_contract_n33_rows
    procedure :: copy_from => mbd_matrix_real_copy_from
    procedure :: move_from => mbd_matrix_real_move_from
    procedure :: init_from => mbd_matrix_real_init_from
    procedure :: alloc_from => mbd_matrix_real_alloc_from
end type

type :: mbd_matrix_complex
    complex(dp), allocatable :: val(:, :)
    type(mbd_index) :: idx
#ifdef WITH_SCALAPACK
    type(mbd_blacs_desc) :: blacs
#endif
    contains
    procedure :: siz => mbd_matrix_complex_siz
    procedure :: init => mbd_matrix_complex_init
    procedure :: add => mbd_matrix_complex_add
    procedure :: add_diag => mbd_matrix_complex_add_diag
    procedure :: add_diag_scalar => mbd_matrix_complex_add_diag_scalar
    procedure :: mult_cross => mbd_matrix_complex_mult_cross
    procedure :: mult_rows => mbd_matrix_complex_mult_rows
    procedure :: mult_cols_3n => mbd_matrix_complex_mult_cols_3n
    procedure :: mult_col => mbd_matrix_complex_mult_col
    ! procedure :: mmul => mbd_matrix_complex_mmul
    ! procedure :: invh => mbd_matrix_complex_invh
    procedure :: eigh => mbd_matrix_complex_eigh
    procedure :: eigvals => mbd_matrix_complex_eigvals
    procedure :: eigvalsh => mbd_matrix_complex_eigvalsh
    procedure :: contract_n_transp => mbd_matrix_complex_contract_n_transp
    procedure :: contract_n33diag_cols => mbd_matrix_complex_contract_n33diag_cols
    procedure :: contract_n33_rows => mbd_matrix_complex_contract_n33_rows
    procedure :: copy_from => mbd_matrix_complex_copy_from
    procedure :: move_from => mbd_matrix_complex_move_from
    procedure :: init_from => mbd_matrix_complex_init_from
    procedure :: alloc_from => mbd_matrix_complex_alloc_from
end type

interface contract_cross_33
    module procedure contract_cross_33_real
    module procedure contract_cross_33_complex
end interface

contains

#endif

#ifndef MBD_TYPE
#define MBD_TYPE 0
#endif

#if MBD_TYPE == 0
integer function mbd_matrix_real_siz(this, ndim) result(siz)
    class(mbd_matrix_real), intent(in) :: this
#elif MBD_TYPE == 1
integer function mbd_matrix_complex_siz(this, ndim) result(siz)
    class(mbd_matrix_complex), intent(in) :: this
#endif
    integer, intent(in) :: ndim

    siz = size(this%val, ndim)
end function

#if MBD_TYPE == 0
#   ifdef WITH_SCALAPACK
subroutine mbd_matrix_real_init(this, idx, blacs)
#   else
subroutine mbd_matrix_real_init(this, idx)
#   endif
    class(mbd_matrix_real), intent(out) :: this
#elif MBD_TYPE == 1
#   ifdef WITH_SCALAPACK
subroutine mbd_matrix_complex_init(this, idx, blacs)
#   else
subroutine mbd_matrix_complex_init(this, idx)
#   endif
    class(mbd_matrix_complex), intent(out) :: this
#endif
    type(mbd_index), intent(in) :: idx
#ifdef WITH_SCALAPACK
    type(mbd_blacs_desc), intent(in) :: blacs
#endif

    this%idx = idx
#ifdef WITH_SCALAPACK
    this%blacs = blacs
#endif
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_init_from(this, other)
    class(mbd_matrix_real), intent(out) :: this
    type(mbd_matrix_real), intent(in) :: other
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_init_from(this, other)
    class(mbd_matrix_complex), intent(out) :: this
    type(mbd_matrix_complex), intent(in) :: other
#endif

    this%idx = other%idx
#ifdef WITH_SCALAPACK
    this%blacs = other%blacs
#endif
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_copy_from(this, other)
    class(mbd_matrix_real), intent(out) :: this
    type(mbd_matrix_real), intent(in) :: other
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_copy_from(this, other)
    class(mbd_matrix_complex), intent(out) :: this
    type(mbd_matrix_complex), intent(in) :: other
#endif
    this%val = other%val
    this%idx = other%idx
#ifdef WITH_SCALAPACK
    this%blacs = other%blacs
#endif
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_move_from(this, other)
    class(mbd_matrix_real), intent(out) :: this
    type(mbd_matrix_real), intent(inout) :: other
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_move_from(this, other)
    class(mbd_matrix_complex), intent(out) :: this
    type(mbd_matrix_complex), intent(inout) :: other
#endif

    call move_alloc(other%val, this%val)
    this%idx = other%idx
#ifdef WITH_SCALAPACK
    this%blacs = other%blacs
#endif
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_alloc_from(this, other)
    class(mbd_matrix_real), intent(out) :: this
    type(mbd_matrix_real), intent(in) :: other
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_alloc_from(this, other)
    class(mbd_matrix_complex), intent(out) :: this
    type(mbd_matrix_complex), intent(in) :: other
#endif

    integer :: n1, n2

    n1 = other%siz(1)
    n2 = other%siz(2)
    allocate (this%val(n1, n2))
    this%idx = other%idx
#ifdef WITH_SCALAPACK
    this%blacs = other%blacs
#endif
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_add(this, other)
    class(mbd_matrix_real), intent(inout) :: this
    class(mbd_matrix_real), intent(in) :: other
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_add(this, other)
    class(mbd_matrix_complex), intent(inout) :: this
    class(mbd_matrix_complex), intent(in) :: other
#endif

    this%val = this%val + other%val
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_add_diag_scalar(this, d)
    class(mbd_matrix_real), intent(inout) :: this
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_add_diag_scalar(this, d)
    class(mbd_matrix_complex), intent(inout) :: this
#endif
    real(dp), intent(in) :: d

    integer :: i

    call this%add_diag([(d, i = 1, this%idx%n_atoms)])
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_add_diag(this, d)
    class(mbd_matrix_real), intent(inout) :: this
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_add_diag(this, d)
    class(mbd_matrix_complex), intent(inout) :: this
#endif
    real(dp), intent(in) :: d(:)

    integer :: my_i_atom, my_j_atom, i

    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_diag => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                if (i_atom /= j_atom) cycle
                do i = 1, 3
                    this_diag(i, i) = this_diag(i, i) + d(i_atom)
                end do
            end associate
        end do
    end do
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_mult_cross(this, b, c)
    class(mbd_matrix_real), intent(inout) :: this
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_mult_cross(this, b, c)
    class(mbd_matrix_complex), intent(inout) :: this
#endif
    real(dp), intent(in) :: b(:)
    real(dp), intent(in), optional :: c(:)

    integer :: my_i_atom, my_j_atom

    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
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
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_mult_rows(this, b)
    class(mbd_matrix_real), intent(inout) :: this
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_mult_rows(this, b)
    class(mbd_matrix_complex), intent(inout) :: this
#endif
    real(dp), intent(in) :: b(:)

    integer :: my_i_atom

    do my_i_atom = 1, size(this%idx%i_atom)
        associate ( &
                i_atom => this%idx%i_atom(my_i_atom), &
                this_sub => this%val(3*(my_i_atom-1)+1:, :) &
        )
            this_sub(:3, :) = this_sub(:3, :)*b(i_atom)
        end associate
    end do
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_mult_cols_3n(this, b)
    class(mbd_matrix_real), intent(inout) :: this
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_mult_cols_3n(this, b)
    class(mbd_matrix_complex), intent(inout) :: this
#endif
    real(dp), intent(in) :: b(:)

    integer :: my_j_atom, i

    do my_j_atom = 1, size(this%idx%j_atom)
        associate ( &
                b_sub => b(3*(this%idx%j_atom(my_j_atom)-1)+1:), &
                this_sub => this%val(:, 3*(my_j_atom-1)+1:) &
        )
            forall (i = 1:3) this_sub(:, i) = this_sub(:, i)*b_sub(i)
        end associate
    end do
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_mult_col(this, idx, a)
    class(mbd_matrix_real), intent(inout) :: this
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_mult_col(this, idx, a)
    class(mbd_matrix_complex), intent(inout) :: this
#endif
    integer, intent(in) :: idx
    real(dp), intent(in) :: a(:)

    integer :: my_i_atom, my_j_atom

    do my_j_atom = 1, size(this%idx%j_atom)
        if (this%idx%j_atom(my_j_atom) /= idx) cycle
        do my_i_atom = 1, size(this%idx%i_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                this_sub(:3, :3) = this_sub(:3, :3)*a(i_atom)
            end associate
        end do
    end do
end subroutine

#if MBD_TYPE == 0
subroutine mbd_matrix_real_eigh(A, eigs, exc, src, vals_only)
    class(mbd_matrix_real), intent(inout) :: A
    type(mbd_matrix_real), intent(in), optional :: src
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_eigh(A, eigs, exc, src, vals_only)
    class(mbd_matrix_complex), intent(inout) :: A
    type(mbd_matrix_complex), intent(in), optional :: src
#endif
    real(dp), intent(out) :: eigs(:)
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: vals_only

#ifdef WITH_SCALAPACK
    if (.not. A%idx%parallel) then
        call eigh(A%val, eigs, exc, src%val, vals_only)
    else
#   if MBD_TYPE == 0
        call peigh(A%val, A%blacs, eigs, exc, src%val, vals_only)
#   elif MBD_TYPE == 1
        exc%code = MBD_EXC_UNIMPL
        exc%msg = 'Complex matrix diagonalization not implemented for scalapack'
#   endif
    endif
#else
    call eigh(A%val, eigs, exc, src%val, vals_only)
#endif
end subroutine

#if MBD_TYPE == 0
function mbd_matrix_real_eigvalsh(A, exc, destroy) result(eigs)
    class(mbd_matrix_real), target, intent(in) :: A
#elif MBD_TYPE == 1
function mbd_matrix_complex_eigvalsh(A, exc, destroy) result(eigs)
    class(mbd_matrix_complex), target, intent(in) :: A
#endif
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigs(3*A%idx%n_atoms)

#ifdef WITH_SCALAPACK
    if (.not. A%idx%parallel) then
        eigs = eigvalsh(A%val, exc, destroy)
    else
#   if MBD_TYPE == 0
        eigs = peigvalsh(A%val, A%blacs, exc, destroy)
#   elif MBD_TYPE == 1
        exc%code = MBD_EXC_UNIMPL
        exc%msg = 'Complex matrix diagonalization not implemented for scalapack'
#   endif
    end if
#else
    eigs = eigvalsh(A%val, exc, destroy)
#endif
end function

#if MBD_TYPE == 0
function mbd_matrix_real_eigvals(A, exc, destroy) result(eigs)
    class(mbd_matrix_real), target, intent(in) :: A
#elif MBD_TYPE == 1
function mbd_matrix_complex_eigvals(A, exc, destroy) result(eigs)
    class(mbd_matrix_complex), target, intent(in) :: A
#endif
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigs(3*A%idx%n_atoms)

#ifdef WITH_SCALAPACK
    if (A%idx%parallel) then
        exc%code = MBD_EXC_UNIMPL
        exc%msg = 'Complex matrix diagonalization not implemented for scalapack'
    else
        eigs = eigvals(A%val, exc, destroy)
    end if
#else
    eigs = eigvals(A%val, exc, destroy)
#endif
end function

#if MBD_TYPE == 0
subroutine mbd_matrix_real_contract_n_transp(this, dir, res)
    class(mbd_matrix_real), intent(in) :: this
    real(dp), intent(out), target :: res(:, :)
#elif MBD_TYPE == 1
subroutine mbd_matrix_complex_contract_n_transp(this, dir, res)
    class(mbd_matrix_complex), intent(in) :: this
    complex(dp), intent(out), target :: res(:, :)
#endif
    character(len=*), intent(in) :: dir

    integer :: my_i_atom, my_j_atom
#if MBD_TYPE == 0
    real(dp), pointer :: res_sub(:, :)
#elif MBD_TYPE == 1
    complex(dp), pointer :: res_sub(:, :)
#endif

    res(:, :) = 0d0
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            select case (dir(1:1))
            case ('R')
                res_sub => res(:, 3*(this%idx%i_atom(my_i_atom)-1)+1:)
            case ('C')
                res_sub => res(3*(this%idx%j_atom(my_j_atom)-1)+1:, :)
            end select
            associate ( &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                res_sub(:3, :3) = res_sub(:3, :3) + transpose(this_sub(:3, :3))
            end associate
        end do
    end do
#ifdef WITH_SCALAPACK
    if (this%idx%parallel) call all_reduce(res, this%blacs)
#endif
end subroutine

#if MBD_TYPE == 0
function contract_cross_33_real(k_atom, A, A_prime, B, B_prime) result(res)
    type(mbd_matrix_real), intent(in) :: A, B
    real(dp), intent(in) :: A_prime(:, :), B_prime(:, :)
    real(dp) :: res(A%idx%n_atoms)
#elif MBD_TYPE == 1
function contract_cross_33_complex(k_atom, A, A_prime, B, B_prime) result(res)
    type(mbd_matrix_complex), intent(in) :: A, B
    complex(dp), intent(in) :: A_prime(:, :), B_prime(:, :)
    complex(dp) :: res(A%idx%n_atoms)
#endif
    integer, intent(in) :: k_atom

    integer :: my_i_atom, my_j_atom, i_atom, j_atom

    res(:) = 0d0
    my_i_atom = findval(A%idx%i_atom, k_atom)
    if (my_i_atom > 0) then
        do my_j_atom = 1, size(A%idx%j_atom)
            j_atom = A%idx%j_atom(my_j_atom)
            associate ( &
                    A_sub => A%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:), &
                    A_prime_sub => A_prime(:, 3*(j_atom-1)+1:) &
            )
                res(j_atom) = -1d0/3*sum(A_sub(:3, :3)*A_prime_sub(:, :3))
            end associate
        end do
    end if
    my_j_atom = findval(A%idx%j_atom, k_atom)
    if (my_j_atom > 0) then
        do my_i_atom = 1, size(A%idx%i_atom)
            i_atom = A%idx%i_atom(my_i_atom)
            associate ( &
                    B_sub => B%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:), &
                    B_prime_sub => B_prime(3*(i_atom-1)+1:, :) &
            )
                res(i_atom) = res(i_atom) + &
                    (-1d0/3)*sum(B_prime_sub(:3, :)*B_sub(:3, :3))
            end associate
        end do
    end if
#ifdef WITH_SCALAPACK
    if (A%idx%parallel) call all_reduce(res, A%blacs)
#endif
end function

#if MBD_TYPE == 0
function mbd_matrix_real_contract_n33diag_cols(A) result(res)
    class(mbd_matrix_real), intent(in) :: A
    real(dp) :: res(A%idx%n_atoms)
#elif MBD_TYPE == 1
function mbd_matrix_complex_contract_n33diag_cols(A) result(res)
    class(mbd_matrix_complex), intent(in) :: A
    complex(dp) :: res(A%idx%n_atoms)
#endif

    integer :: i_xyz, my_j_atom, j_atom

    res(:) = 0d0
    do my_j_atom = 1, size(A%idx%j_atom)
        j_atom = A%idx%j_atom(my_j_atom)
        do i_xyz = 1, 3
            res(j_atom) = res(j_atom) + &
                sum(A%val(i_xyz::3, 3*(my_j_atom-1)+i_xyz))
        end do
    end do
    res = res/3
#ifdef WITH_SCALAPACK
    if (A%idx%parallel) call all_reduce(res, A%blacs)
#endif
end function

#if MBD_TYPE == 0
function mbd_matrix_real_contract_n33_rows(A) result(res)
    class(mbd_matrix_real), intent(in) :: A
    real(dp) :: res(A%idx%n_atoms)
#elif MBD_TYPE == 1
function mbd_matrix_complex_contract_n33_rows(A) result(res)
    class(mbd_matrix_complex), intent(in) :: A
    complex(dp) :: res(A%idx%n_atoms)
#endif

    integer :: my_i_atom, i_atom

    res(:) = 0d0
    do my_i_atom = 1, size(A%idx%i_atom)
        i_atom = A%idx%i_atom(my_i_atom)
        associate (A_sub => A%val(3*(my_i_atom-1)+1:, :))
            res(i_atom) = res(i_atom) + sum(A_sub(:3, :))
        end associate
    end do
#ifdef WITH_SCALAPACK
    if (A%idx%parallel) call all_reduce(res, A%blacs)
#endif
end function

#undef MBD_TYPE
#ifndef MBD_INCLUDED
#define MBD_INCLUDED
#define MBD_TYPE 1
#include "mbd_matrix_type.F90"
#undef MBD_INCLUDED

! #if MBD_TYPE == 0
type(mbd_matrix_real) function mbd_matrix_real_mmul( &
        A, B, transA, transB) result(C)
    class(mbd_matrix_real), intent(in) :: A
    type(mbd_matrix_real), intent(in) :: B
! #elif MBD_TYPE == 1
! type(mbd_matrix_complex) function mbd_matrix_complex_mmul( &
!         A, B, transA, transB) result(C)
!     class(mbd_matrix_complex), intent(in) :: A
!     type(mbd_matrix_complex), intent(in) :: B
! #endif
    logical, intent(in), optional :: transA, transB

    C%idx = A%idx
#ifdef WITH_SCALAPACK
    C%blacs = A%blacs
    if (.not. A%idx%parallel) then
        C%val = mmul(A%val, B%val, transA, transB)
    else
        C%val = pmmul(A%val, A%blacs, B%val, B%blacs, transA, transB, C%blacs)
    end if
#else
    C%val = mmul(A%val, B%val, transA, transB)
#endif
end function

subroutine mbd_matrix_real_invh(A, exc, src)
    class(mbd_matrix_real), intent(inout) :: A
    type(mbd_matrix_real), intent(in), optional :: src
    type(exception), intent(out), optional :: exc

#ifdef WITH_SCALAPACK
    if (.not. A%idx%parallel) then
        if (present(src)) then
            call invh(A%val, exc, src%val)
        else
            call invh(A%val, exc)
        end if
    else
        if (present(src)) then
            call pinvh(A%val, A%blacs, exc, src%val)
        else
            call pinvh(A%val, A%blacs, exc)
        end if
    end if
#else
    if (present(src)) then
        call invh(A%val, exc, src%val)
    else
        call invh(A%val, exc)
    end if
#endif
end subroutine

end module

#endif
