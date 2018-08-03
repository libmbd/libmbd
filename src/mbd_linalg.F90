! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_linalg

use mbd_common, only: tostr, dp, exception
use mbd_types, only: mat3n3n
use mbd_parallel, only: mbd_blacs

implicit none

private
public :: inv, invh, inverse, eig, eigh, eigvals, eigvalsh, solve, cprod

interface inv
    module procedure inv_re_
    module procedure inv_cplx_
end interface

interface invh
    module procedure invh_re_
    module procedure invh_mat3n3n_
end interface

interface eig
    module procedure eig_re_
    module procedure eig_cplx_
end interface

interface eigh
    module procedure eigh_re_
    module procedure eigh_cplx_
    module procedure eigh_mat3n3n_
end interface

interface eigvals
    module procedure eigvals_re_
    module procedure eigvals_cplx_
    module procedure eigvals_mat3n3n_
end interface

interface eigvalsh
    module procedure eigvalsh_re_
    module procedure eigvalsh_cplx_
    module procedure eigvalsh_mat3n3n_
end interface

external :: ZHEEV, DGEEV, DSYEV, DGETRF, DGETRI, DGESV, ZGETRF, ZGETRI, &
    ZGEEV, ZGEEB, DSYTRI, DSYTRF
#ifdef WITH_SCALAPACK
external :: PDSYEV, PDGETRF, PDGETRI
#endif

contains

subroutine inv_re_(A, exc, src)
    real(dp), intent(inout) :: A(:, :)
    type(exception), intent(out), optional :: exc
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
            exc%label = 'linalg'
            exc%origin = 'DGETRF'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    call DGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr))
    call DGETRI(n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'DGETRI'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine

subroutine inv_cplx_(A, exc, src)
    complex(dp), intent(inout) :: A(:, :)
    type(exception), intent(out), optional :: exc
    complex(dp), intent(in), optional :: src(:, :)

    integer, allocatable :: i_pivot(:)
    complex(dp), allocatable :: work_arr(:)
    integer :: n, n_work_arr, error_flag
    complex(dp) :: n_work_arr_optim

    n = size(A, 1)
    if (present(src)) A = src
    allocate (i_pivot(n))
    call ZGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'ZGETRF'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    call ZGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(dble(n_work_arr_optim))
    allocate (work_arr(n_work_arr))
    call ZGETRI(n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'ZGETRI'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine

subroutine invh_re_(A, exc, src)
    real(dp), intent(inout) :: A(:, :)
    type(exception), intent(out), optional :: exc
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
    call DSYTRF('U', n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
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
            exc%label = 'linalg'
            exc%origin = 'DSYTRI'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    call fill_tril(A)
end subroutine

#ifdef WITH_SCALAPACK
subroutine pinvh_re_(A, blacs, exc, src)
    real(dp), intent(inout) :: A(:, :)
    type(mbd_blacs), intent(in) :: blacs
    type(exception), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)

    integer, allocatable :: i_pivot(:), iwork_arr(:)
    real(dp), allocatable :: work_arr(:)
    integer :: n, n_work_arr, error_flag, n_iwork_arr
    real(dp) :: n_work_arr_optim

    n = 3*blacs%n_atoms
    if (n == 0) return
    if (present(src)) A = src
    allocate (i_pivot(n))
    call PDGETRF(n, n, A, 1, 1, blacs%desc, i_pivot, error_flag)
    ! call DSYTRF('U', n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'PDGETRF'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    call PDGETRI( &
        n, A, 1, 1, blacs%desc, i_pivot, &
        n_work_arr_optim, -1, n_iwork_arr, -1, error_flag &
    )
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr), iwork_arr(n_iwork_arr))
    call PDGETRI( &
        n, A, 1, 1, blacs%desc, i_pivot, &
        work_arr, n_work_arr, iwork_arr, n_iwork_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'PDSYTRI'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine
#endif

subroutine invh_mat3n3n_(A, exc, src)
    type(mat3n3n), intent(inout) :: A
    type(exception), intent(out), optional :: exc
    type(mat3n3n), intent(in), optional :: src

#ifndef WITH_SCALAPACK
    call invh_re_(A%re, exc, src%re)
#else
    call pinvh_re_(A%re, A%blacs, exc, src%re)
#endif
end subroutine

function inverse(A, exc)
    real(dp), intent(in) :: A(:, :)
    type(exception), intent(out), optional :: exc
    real(dp) :: inverse(size(A, 1), size(A, 2))

    call inv(inverse, exc, src=A)
end function

subroutine eigh_mat3n3n_(A, eigs, exc, src, vals_only)
    type(mat3n3n), intent(inout) :: A
    real(dp), intent(out) :: eigs(:)
    type(exception), intent(out), optional :: exc
    type(mat3n3n), intent(in), optional :: src
    logical, intent(in), optional :: vals_only

    if (allocated(A%re)) then
#ifndef WITH_SCALAPACK
        call eigh_re_(A%re, eigs, exc, src%re, vals_only)
#else
        call peigh_re_(A%re, A%blacs, eigs, exc, src%re, vals_only)
#endif
    else
        call eigh(A%cplx, eigs, exc, src%cplx, vals_only)
    end if
end subroutine

subroutine eigh_re_(A, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    real(dp), intent(out) :: eigs(:)
    type(exception), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only

    real(dp), allocatable :: work_arr(:)
    real(dp) :: n_work_arr
    integer :: error_flag, n

    n = size(A, 1)
    if (present(src)) A = src
    call DSYEV(mode(vals_only), 'U', n, A, n, eigs, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr)))
    call DSYEV(mode(vals_only), 'U', n, A, n, eigs, work_arr, size(work_arr), error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'DSYEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine

#ifdef WITH_SCALAPACK
subroutine peigh_re_(A, blacs, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    type(mbd_blacs), intent(in) :: blacs
    real(dp), intent(out) :: eigs(:)
    type(exception), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only

    real(dp), allocatable :: work_arr(:), vectors(:, :)
    real(dp) :: n_work_arr
    integer :: error_flag, n

    n = 3*blacs%n_atoms
    if (present(src)) A = src
    if (mode(vals_only) == 'V') then
        allocate (vectors(n, n))
    else
        allocate (vectors(1, 1))
    end if
    call PDSYEV( &
        mode(vals_only), 'U', n, A, 1, 1, blacs%desc, eigs, vectors, &
        1, 1, blacs%desc, n_work_arr, -1, error_flag &
    )
    allocate (work_arr(nint(n_work_arr)))
    call PDSYEV( &
        mode(vals_only), 'U', n, A, 1, 1, blacs%desc, eigs, vectors, &
        1, 1, blacs%desc, work_arr, size(work_arr), error_flag &
    )
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'PDSYEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    if (mode(vals_only) == 'V') A = vectors
end subroutine
#endif

subroutine eig_re_(A, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    complex(dp), intent(out) :: eigs(:)
    type(exception), intent(out), optional :: exc
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
        vectors, n, work_arr, size(work_arr), error_flag &
    )
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'DGEEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    eigs = cmplx(eigs_r, eigs_i, 8)
    if (mode(vals_only) == 'V') A = vectors
end subroutine

function eigvals_re_(A, exc, destroy)
    real(dp), target, intent(in) :: A(:, :)
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigvals_re_(size(A, 1))

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
    call eig_re_(A_p, eigvals_re_, exc, vals_only=.true.)
end function

function eigvals_cplx_(A, exc, destroy)
    complex(dp), target, intent(in) :: A(:, :)
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigvals_cplx_(size(A, 1))

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
    call eig_cplx_(A_p, eigvals_cplx_, exc, vals_only=.true.)
end function

subroutine eigh_cplx_(A, eigs, exc, src, vals_only)
    complex(dp), intent(inout) :: A(:, :)
    real(dp), intent(out) :: eigs(:)
    type(exception), intent(out), optional :: exc
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
    call ZHEEV(mode(vals_only), 'U', n, A, n, eigs, work, lwork, rwork, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'ZHEEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine

subroutine eig_cplx_(A, eigs, exc, src, vals_only)
    complex(dp), intent(inout) :: A(:, :)
    complex(dp), intent(out) :: eigs(:)
    type(exception), intent(out), optional :: exc
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
    if (mode(vals_only) == 'V') allocate (vectors(n, n))
    call ZGEEV( &
        'N', mode(vals_only), n, A, n, eigs, dummy, 1, &
        vectors, n, lwork_arr, -1, rwork, error_flag &
    )
    lwork = nint(dble(lwork_arr))
    allocate (work(lwork))
    call ZGEEV( &
        'N', mode(vals_only), n, A, n, eigs, dummy, 1, &
        vectors, n, work, lwork, rwork, error_flag &
    )
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'ZGEEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    if (mode(vals_only) == 'V') A = vectors
end subroutine

function eigvalsh_re_(A, exc, destroy)
    real(dp), target, intent(in) :: A(:, :)
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigvalsh_re_(size(A, 1))

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
    call eigh_re_(A_p, eigvalsh_re_, exc, vals_only=.true.)
end function

#ifdef WITH_SCALAPACK
function peigvalsh_re_(A, blacs, exc, destroy)
    real(dp), target, intent(in) :: A(:, :)
    type(mbd_blacs), intent(in) :: blacs
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: peigvalsh_re_(3*blacs%n_atoms)

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
    call peigh_re_(A_p, blacs, peigvalsh_re_, exc, vals_only=.true.)
end function
#endif

function eigvalsh_cplx_(A, exc, destroy)
    complex(dp), target, intent(in) :: A(:, :)
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigvalsh_cplx_(size(A, 1))

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
    call eigh_cplx_(A_p, eigvalsh_cplx_, exc, vals_only=.true.)
end function

function eigvalsh_mat3n3n_(A, exc, destroy)
    type(mat3n3n), target, intent(in) :: A
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigvalsh_mat3n3n_(3*A%blacs%n_atoms)

    if (allocated(A%re)) then
#ifndef WITH_SCALAPACK
        eigvalsh_mat3n3n_ = eigvalsh_re_(A%re, exc, destroy)
#else
        eigvalsh_mat3n3n_ = peigvalsh_re_(A%re, A%blacs, exc, destroy)
#endif
    else
        eigvalsh_mat3n3n_ = eigvalsh(A%cplx, exc, destroy)
    end if
end function

function eigvals_mat3n3n_(A, exc, destroy)
    type(mat3n3n), target, intent(in) :: A
    type(exception), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigvals_mat3n3n_(3*A%blacs%n_atoms)

    if (allocated(A%re)) then
        eigvals_mat3n3n_ = eigvals(A%re, exc, destroy)
    else
        eigvals_mat3n3n_ = eigvals(A%cplx, exc, destroy)
    end if
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

function solve(A, b, exc) result(x)
    real(dp), intent(in) :: A(:, :), b(:)
    real(dp) :: x(size(b))
    type(exception), intent(out), optional :: exc

    real(dp) :: A_(size(b), size(b))
    integer :: i_pivot(size(b))
    integer :: n
    integer :: error_flag

    A_ = A
    x = b
    n = size(b)
    call DGESV(n, 1, A_, n, i_pivot, x, n, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%label = 'linalg'
            exc%origin = 'DGESV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end function

function cprod(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))

    integer :: i, j

    do i = 1, size(a)
        do j = 1, size(b)
            c(i, j) = a(i)*b(j)
        end do
    end do
end function

character(len=1) function mode(vals_only)
    logical, intent(in), optional :: vals_only

    mode = 'V'
    if (present(vals_only)) then
        if (vals_only) mode = 'N'
    end if
end function

end module mbd_linalg
