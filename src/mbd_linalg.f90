! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_linalg

use mbd_common, only: tostr, dp
use mbd_interface, only: print_log, print_error

implicit none

private
public :: diag, invert, inverted, diagonalize, sdiagonalize, diagonalized, &
    sdiagonalized, solve_lin_sys, eye, operator(.cprod.), sinvert

interface operator(.cprod.)
    module procedure cart_prod_
end interface

interface diag
    module procedure get_diag_
    module procedure get_diag_cmplx_
    module procedure make_diag_
end interface

interface invert
    module procedure invert_ge_dble_
    module procedure invert_ge_cmplx_
end interface

interface sinvert
    module procedure invert_sym_dble_
    ! module procedure invert_he_cmplx_
end interface

interface diagonalize
    module procedure diagonalize_ge_dble_
    module procedure diagonalize_ge_cmplx_
end interface

interface sdiagonalize
    module procedure diagonalize_sym_dble_
    module procedure diagonalize_he_cmplx_
end interface

interface diagonalized
    module procedure diagonalized_ge_dble_
end interface

interface sdiagonalized
    module procedure diagonalized_sym_dble_
end interface

external :: ZHEEV, DGEEV, DSYEV, DGETRF, DGETRI, DGESV, ZGETRF, ZGETRI, &
    ZGEEV, ZGEEB, DSYTRI, DSYTRF

contains


function eye(n) result(A)
    integer, intent(in) :: n
    real(dp) :: A(n, n)

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:n) A(i, i) = 1.d0
end function


subroutine invert_ge_dble_(A)
    real(dp), intent(inout) :: A(:, :)

    integer :: i_pivot(size(A, 1))
    real(dp), allocatable :: work_arr(:)
    integer :: n
    integer :: n_work_arr
    real(dp) :: n_work_arr_optim
    integer :: error_flag

    n = size(A, 1)
    if (n == 0) return
    call DGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    call DGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr))
    call DGETRI(n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


subroutine invert_ge_cmplx_(A)
    complex(dp), intent(inout) :: A(:, :)

    integer :: i_pivot(size(A, 1))
    complex(dp), allocatable :: work_arr(:)
    integer :: n
    integer :: n_work_arr
    complex(dp) :: n_work_arr_optim
    integer :: error_flag

    n = size(A, 1)
    call ZGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    call ZGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(dble(n_work_arr_optim))
    allocate (work_arr(n_work_arr))
    call ZGETRI(n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


subroutine invert_sym_dble_(A)
    real(dp), intent(inout) :: A(:, :)

    integer :: i_pivot(size(A, 1))
    real(dp), allocatable :: work_arr(:)
    integer :: n
    integer :: n_work_arr
    real(dp) :: n_work_arr_optim
    integer :: error_flag

    n = size(A, 1)
    if (n == 0) return
    call DSYTRF('U', n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr))
    call DSYTRF('U', n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    deallocate (work_arr)
    allocate (work_arr(n))
    call DSYTRI('U', n, A, n, i_pivot, work_arr, error_flag)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


function inverted(A) result(A_inv)
    real(dp), intent(in) :: A(:, :)
    real(dp) :: A_inv(size(A, 1), size(A, 2))

    A_inv = A
    call invert(A_inv)
end function


subroutine diagonalize_sym_dble_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    real(dp), intent(inout) :: A(:, :)
    real(dp), intent(out) :: eigs(size(A, 1))

    real(dp), allocatable :: work_arr(:)
    integer :: n
    real(dp) :: n_work_arr
    integer :: error_flag

    n = size(A, 1)
    call DSYEV(mode, "U", n, A, n, eigs, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr)))
    call DSYEV(mode, "U", n, A, n, eigs, work_arr, size(work_arr), error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_log( &
            "DSYEV failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


function diagonalized_sym_dble_(A, eigvecs) result(eigs)
    real(dp), intent(in) :: A(:, :)
    real(dp), intent(out), optional, target :: eigvecs(size(A, 1), size(A, 2))
    real(dp) :: eigs(size(A, 1))

    real(dp), pointer :: eigvecs_p(:, :)
    character(len=1) :: mode

    if (present(eigvecs)) then
        mode = 'V'
        eigvecs_p => eigvecs
    else
        mode = 'N'
        allocate (eigvecs_p(size(A, 1), size(A, 2)))
    end if
    eigvecs_p = A
    call sdiagonalize(mode, eigvecs_p, eigs)
    if (.not. present(eigvecs)) then
        deallocate (eigvecs_p)
    end if
end function


subroutine diagonalize_ge_dble_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    real(dp), intent(inout) :: A(:, :)
    complex(dp), intent(out) :: eigs(size(A, 1))

    real(dp), allocatable :: work_arr(:)
    integer :: n
    real(dp) :: n_work_arr
    integer :: error_flag
    real(dp) :: eigs_r(size(A, 1)), eigs_i(size(A, 1))
    real(dp) :: dummy
    real(dp) :: vectors(size(A, 1), size(A, 2))

    n = size(A, 1)
    call DGEEV('N', mode, n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr)))
    call DGEEV('N', mode, n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, work_arr, size(work_arr), error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_log( &
            "DGEEV failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    eigs = cmplx(eigs_r, eigs_i, 8)
    A = vectors
end subroutine


function diagonalized_ge_dble_(A, eigvecs) result(eigs)
    real(dp), intent(in) :: A(:, :)
    real(dp), intent(out), optional, target :: eigvecs(size(A, 1), size(A, 2))
    complex(dp) :: eigs(size(A, 1))

    real(dp), pointer :: eigvecs_p(:, :)
    character(len=1) :: mode

    if (present(eigvecs)) then
        mode = 'V'
        eigvecs_p => eigvecs
    else
        mode = 'N'
        allocate (eigvecs_p(size(A, 1), size(A, 2)))
    end if
    eigvecs_p = A
    call diagonalize(mode, eigvecs_p, eigs)
    if (.not. present(eigvecs)) then
        deallocate (eigvecs_p)
    end if
end function


subroutine diagonalize_he_cmplx_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    complex(dp), intent(inout) :: A(:, :)
    real(dp), intent(out) :: eigs(size(A, 1))

    complex(dp), allocatable :: work(:)
    complex(dp) :: lwork_cmplx
    real(dp), allocatable :: rwork(:)
    integer :: n, lwork
    integer :: error_flag
    integer, external :: ILAENV

    n = size(A, 1)
    allocate (rwork(max(1, 3*n-2)))
    call ZHEEV(mode, "U", n, A, n, eigs, lwork_cmplx, -1, rwork, error_flag)
    lwork = nint(dble(lwork_cmplx))
    allocate (work(lwork))
    call ZHEEV(mode, "U", n, A, n, eigs, work, lwork, rwork, error_flag)
    deallocate (rwork)
    deallocate (work)
    if (error_flag /= 0) then
        call print_error( &
            "ZHEEV failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


subroutine diagonalize_ge_cmplx_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    complex(dp), intent(inout) :: A(:, :)
    complex(dp), intent(out) :: eigs(size(A, 1))

    complex(dp), allocatable :: work(:)
    real(dp) :: rwork(2*size(A, 1))
    integer :: n, lwork
    complex(dp) :: lwork_arr
    integer :: error_flag
    complex(dp) :: dummy
    complex(dp) :: vectors(size(A, 1), size(A, 2))

    n = size(A, 1)
    call ZGEEV('N', mode, n, A, n, eigs, dummy, 1, &
        vectors, n, lwork_arr, -1, rwork, error_flag)
    lwork = nint(dble(lwork_arr))
    allocate (work(lwork))
    call ZGEEV('N', mode, n, A, n, eigs, dummy, 1, &
        vectors, n, work, lwork, rwork, error_flag)
    deallocate (work)
    if (error_flag /= 0) then
        call print_log( &
            "ZGEEV failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    A = vectors
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


function solve_lin_sys(A, b) result(x)
    real(dp), intent(in) :: A(:, :), b(size(A, 1))
    real(dp) :: x(size(b))

    real(dp) :: A_(size(b), size(b))
    integer :: i_pivot(size(b))
    integer :: n
    integer :: error_flag

    A_ = A
    x = b
    n = size(b)
    call DGESV(n, 1, A_, n, i_pivot, x, n, error_flag)
end function

end module mbd_linalg
