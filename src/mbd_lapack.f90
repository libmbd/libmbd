! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_lapack

use mbd_constants
use mbd_utils, only: exception_t, tostr

implicit none

private
public :: mmul, inv, invh, inverse, eig, eigh, eigvals, eigvalsh, det, mode

interface mmul
    module procedure mmul_real
    module procedure mmul_complex
end interface

interface inv
    module procedure inv_real
end interface

interface invh
    module procedure invh_real
end interface

interface eig
    module procedure eig_real
    module procedure eig_complex
end interface

interface eigh
    module procedure eigh_real
    module procedure eigh_complex
end interface

interface eigvals
    module procedure eigvals_real
    module procedure eigvals_complex
end interface

interface eigvalsh
    module procedure eigvalsh_real
    module procedure eigvalsh_complex
end interface

external :: ZHEEV, DGEEV, DSYEV, DGETRF, DGETRI, DGESV, ZGETRF, ZGETRI, &
    ZGEEV, DSYTRI, DSYTRF, DGEMM, ZGEMM

contains

function inverse(A, exc)
    real(dp), intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    real(dp) :: inverse(size(A, 1), size(A, 2))

    call inv_real(inverse, exc, src=A)
end function

function eigvalsh_real(A, exc, destroy) result(eigvals)
    real(dp), target, intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigvals(size(A, 1))

    real(dp), allocatable, target :: A_work(:, :)
    real(dp), pointer :: A_p(:, :)

    nullify (A_p)
    if (present(destroy)) then
        if (destroy) then
            A_p => A
        end if
    end if
    if (.not. associated(A_p)) then
        allocate (A_work(size(A, 1), size(A, 1)), source=A)
        A_p => A_work
    end if
    call eigh_real(A_p, eigvals, exc, vals_only=.true.)
end function

function eigvalsh_complex(A, exc, destroy) result(eigvals)
    complex(dp), target, intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigvals(size(A, 1))

    complex(dp), allocatable, target :: A_work(:, :)
    complex(dp), pointer :: A_p(:, :)

    nullify (A_p)
    if (present(destroy)) then
        if (destroy) then
            A_p => A
        end if
    end if
    if (.not. associated(A_p)) then
        allocate (A_work(size(A, 1), size(A, 1)), source=A)
        A_p => A_work
    end if
    call eigh_complex(A_p, eigvals, exc, vals_only=.true.)
end function

function eigvals_real(A, exc, destroy) result(eigvals)
    real(dp), target, intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigvals(size(A, 1))

    real(dp), allocatable, target :: A_work(:, :)
    real(dp), pointer :: A_p(:, :)

    nullify (A_p)
    if (present(destroy)) then
        if (destroy) then
            A_p => A
        end if
    end if
    if (.not. associated(A_p)) then
        allocate (A_work(size(A, 1), size(A, 1)), source=A)
        A_p => A_work
    end if
    call eig_real(A_p, eigvals, exc, vals_only=.true.)
end function

function eigvals_complex(A, exc, destroy) result(eigvals)
    complex(dp), target, intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigvals(size(A, 1))

    complex(dp), allocatable, target :: A_work(:, :)
    complex(dp), pointer :: A_p(:, :)

    nullify (A_p)
    if (present(destroy)) then
        if (destroy) then
            A_p => A
        end if
    end if
    if (.not. associated(A_p)) then
        allocate (A_work(size(A, 1), size(A, 1)), source=A)
        A_p => A_work
    end if
    call eig_complex(A_p, eigvals, exc, vals_only=.true.)
end function

function mmul_real(A, B, transA, transB) result(C)
    real(dp), intent(in) :: A(:, :), B(:, :)
    character, intent(in), optional :: transA, transB
    real(dp) :: C(size(A, 1), size(B, 2))

    character :: transA_, transB_
    integer :: n

    transA_= 'N'
    transB_ = 'N'
    if (present(transA)) transA_ = transA
    if (present(transB)) transB_ = transB
    n = size(A, 1)
    call DGEMM(transA_, transB_, n, n, n, 1d0, A, n, B, n, 0d0, C, n)
end function

function mmul_complex(A, B, transA, transB) result(C)
    complex(dp), intent(in) :: A(:, :), B(:, :)
    character, intent(in), optional :: transA, transB
    complex(dp) :: C(size(A, 1), size(B, 2))

    character :: transA_, transB_
    integer :: n

    transA_= 'N'
    transB_ = 'N'
    if (present(transA)) transA_ = transA
    if (present(transB)) transB_ = transB
    n = size(A, 1)
    call ZGEMM(transA_, transB_, n, n, n, (1d0, 0d0), A, n, B, n, (0d0, 0d0), C, n)
end function

subroutine inv_real(A, exc, src)
    real(dp), intent(inout) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)

    real(dp), allocatable :: work_arr(:)
    integer, allocatable :: i_pivot(:)
    integer :: n, n_work_arr, error_flag
    real(dp) :: n_work_arr_optim

    n = size(A, 1)
    if (n == 0) return
    if (present(src)) A = src
    allocate (i_pivot(n))
    call DGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DGETRF'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    call DGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr))
    call DGETRI(n, A, n, i_pivot, work_arr(1), n_work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DGETRI'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine

subroutine invh_real(A, exc, src)
    real(dp), intent(inout) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)

    integer, allocatable :: i_pivot(:)
    real(dp), allocatable :: work_arr(:)
    integer :: n, n_work_arr, error_flag
    real(dp) :: n_work_arr_optim

    n = size(A, 1)
    if (n == 0) return
    if (present(src)) A = src
    allocate (i_pivot(n))
    call DSYTRF('U', n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr))
    call DSYTRF('U', n, A, n, i_pivot, work_arr(1), n_work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DSYTRF'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    deallocate (work_arr)
    allocate (work_arr(n))
    call DSYTRI('U', n, A, n, i_pivot, work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DSYTRI'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    call fill_tril(A)
end subroutine

subroutine eigh_real(A, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only

    real(dp), allocatable :: work_arr(:)
    real(dp) :: n_work_arr
    integer :: error_flag, n

    n = size(A, 1)
    if (present(src)) A = src
    call DSYEV(mode(vals_only), 'U', n, A, n, eigs, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr)))
    call DSYEV(mode(vals_only), 'U', n, A, n, eigs, work_arr(1), size(work_arr), error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DSYEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine

subroutine eig_real(A, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    complex(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only

    real(dp) :: n_work_arr, dummy
    integer :: error_flag, n
    real(dp), allocatable :: eigs_r(:), eigs_i(:), vectors(:, :), work_arr(:)

    n = size(A, 1)
    if (present(src)) A = src
    allocate (eigs_r(n), eigs_i(n))
    if (mode(vals_only) == 'V') then
        allocate (vectors(n, n))
    else
        allocate (vectors(1, 1))
    end if
    call DGEEV( &
        'N', mode(vals_only), n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, n_work_arr, -1, error_flag &
    )
    allocate (work_arr(nint(n_work_arr)))
    call DGEEV( &
        'N', mode(vals_only), n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, work_arr(1), size(work_arr), error_flag &
    )
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DGEEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    eigs = cmplx(eigs_r, eigs_i, dp)
    if (mode(vals_only) == 'V') A = vectors
end subroutine

subroutine eigh_complex(A, eigs, exc, src, vals_only)
    complex(dp), intent(inout) :: A(:, :)
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    complex(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only

    complex(dp), allocatable :: work(:)
    complex(dp) :: lwork_cmplx
    real(dp), allocatable :: rwork(:)
    integer :: n, lwork, error_flag

    n = size(A, 1)
    if (present(src)) A = src
    allocate (rwork(max(1, 3*n-2)))
    call ZHEEV(mode(vals_only), 'U', n, A, n, eigs, lwork_cmplx, -1, rwork, error_flag)
    lwork = nint(dble(lwork_cmplx))
    allocate (work(lwork))
    call ZHEEV(mode(vals_only), 'U', n, A, n, eigs, work(1), lwork, rwork, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'ZHEEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine

subroutine eig_complex(A, eigs, exc, src, vals_only)
    complex(dp), intent(inout) :: A(:, :)
    complex(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    complex(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only

    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer :: n, lwork
    complex(dp) :: lwork_arr
    integer :: error_flag
    complex(dp) :: dummy
    complex(dp), allocatable :: vectors(:, :)

    n = size(A, 1)
    if (present(src)) A = src
    allocate (rwork(2*n))
    if (mode(vals_only) == 'V') then
        allocate (vectors(n, n))
    else
        allocate (vectors(1, 1))
    end if
    call ZGEEV( &
        'N', mode(vals_only), n, A, n, eigs, dummy, 1, &
        vectors, n, lwork_arr, -1, rwork, error_flag &
    )
    lwork = nint(dble(lwork_arr))
    allocate (work(lwork))
    call ZGEEV( &
        'N', mode(vals_only), n, A, n, eigs, dummy, 1, &
        vectors, n, work(1), lwork, rwork, error_flag &
    )
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'ZGEEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    if (mode(vals_only) == 'V') A = vectors
end subroutine

real(dp) function det(A) result(D)
    real(dp), intent(in) :: A(:, :)

    integer :: n, i, info
    real(dp), allocatable :: LU(:, :)
    integer, allocatable :: ipiv(:)

    n = size(A, 1)
    allocate (ipiv(n))
    LU = A
    call DGETRF(n, n, LU, n, ipiv, info)
    D = product([(LU(i, i), i = 1, n)])
end function

subroutine fill_tril(A)
    real(dp), intent(inout) :: A(:, :)

    integer :: i, j

    do i = 1, size(A, 1)
        do j = i+1, size(A, 1)
            A(j, i) = A(i, j)
        end do
    end do
end subroutine

character(len=1) function mode(vals_only)
    logical, intent(in), optional :: vals_only

    mode = 'V'
    if (present(vals_only)) then
        if (vals_only) mode = 'N'
    end if
end function

end module
