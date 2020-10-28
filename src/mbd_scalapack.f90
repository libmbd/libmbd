! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_scalapack

use mbd_constants
use mbd_lapack, only: mode
use mbd_blacs, only: blacs_desc_t
use mbd_utils, only: exception_t, tostr

implicit none

private
public :: pmmul, pinvh, peigh, peigvalsh

interface pmmul
    module procedure pmmul_real
    module procedure pmmul_complex
end interface

interface pinvh
    module procedure pinvh_real
end interface

interface peigh
    module procedure peigh_real
    module procedure peigh_complex
end interface

interface peigvalsh
    module procedure peigvalsh_real
    module procedure peigvalsh_complex
end interface

external :: PDSYEV, PZHEEV, PDGETRF, PDGETRI, PDGEMM, PZGEMM

contains

subroutine pinvh_real(A, blacs, exc, src)
    real(dp), intent(inout) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)

    integer, allocatable :: i_pivot(:), iwork_arr(:)
    real(dp), allocatable :: work_arr(:)
    integer :: n, n_work_arr, error_flag, n_iwork_arr(1)
    real(dp) :: n_work_arr_optim(1)

    n = 3*blacs%n_atoms
    if (n == 0) return
    if (present(src)) A = src
    allocate (i_pivot(n))
    call PDGETRF(n, n, A, 1, 1, blacs%desc, i_pivot, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'PDGETRF'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    call PDGETRI( &
        n, A, 1, 1, blacs%desc, i_pivot, &
        n_work_arr_optim, -1, n_iwork_arr, -1, error_flag &
    )
    n_work_arr = nint(n_work_arr_optim(1))
    allocate (work_arr(n_work_arr), iwork_arr(n_iwork_arr(1)))
    call PDGETRI( &
        n, A, 1, 1, blacs%desc, i_pivot, &
        work_arr, n_work_arr, iwork_arr, n_iwork_arr(1), error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'PDSYTRI'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
end subroutine

function pmmul_real(A, blacsA, B, blacsB, transA, transB, blacsC) result(C)
    real(dp), intent(in) :: A(:, :), B(:, :)
    type(blacs_desc_t), intent(in) :: blacsA, blacsB, blacsC
    character, intent(in), optional :: transA, transB
    real(dp) :: C(size(A, 1), size(B, 2))

    character :: transA_, transB_
    integer :: n

    transA_= 'N'
    transB_ = 'N'
    if (present(transA)) transA_ = transA
    if (present(transB)) transB_ = transB
    n = 3*blacsA%n_atoms
    call PDGEMM( &
        transA_, transB_, n, n, n, 1d0, A, 1, 1, blacsA%desc, &
        B, 1, 1, blacsB%desc, 0d0, C, 1, 1, blacsC%desc &
    )
end function

function pmmul_complex(A, blacsA, B, blacsB, transA, transB, blacsC) result(C)
    complex(dp), intent(in) :: A(:, :), B(:, :)
    type(blacs_desc_t), intent(in) :: blacsA, blacsB, blacsC
    character, intent(in), optional :: transA, transB
    complex(dp) :: C(size(A, 1), size(B, 2))

    character :: transA_, transB_
    integer :: n

    transA_= 'N'
    transB_ = 'N'
    if (present(transA)) transA_ = transA
    if (present(transB)) transB_ = transB
    n = 3*blacsA%n_atoms
    call PZGEMM( &
        transA_, transB_, n, n, n, (1d0, 0d0), A, 1, 1, blacsA%desc, &
        B, 1, 1, blacsB%desc, (0d0, 0d0), C, 1, 1, blacsC%desc &
    )
end function

subroutine peigh_real(A, blacs, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only

    real(dp), allocatable :: work_arr(:), vectors(:, :)
    real(dp) :: n_work_arr(1)
    integer :: error_flag, n

    n = 3*blacs%n_atoms
    if (present(src)) A = src
    if (mode(vals_only) == 'V') then
        allocate (vectors(size(A, 1), size(A, 2)))
    else
        allocate (vectors(1, 1))
    end if
    call PDSYEV( &
        mode(vals_only), 'U', n, A, 1, 1, blacs%desc, eigs, vectors, &
        1, 1, blacs%desc, n_work_arr, -1, error_flag &
    )
    allocate (work_arr(nint(n_work_arr(1))))
    call PDSYEV( &
        mode(vals_only), 'U', n, A, 1, 1, blacs%desc, eigs, vectors, &
        1, 1, blacs%desc, work_arr(1), size(work_arr), error_flag &
    )
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'PDSYEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    if (mode(vals_only) == 'V') A = vectors
end subroutine

subroutine peigh_complex(A, blacs, eigs, exc, src, vals_only)
    complex(dp), intent(inout) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    complex(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only

    complex(dp), allocatable :: work_arr(:), vectors(:, :)
    integer :: n_work_arr, n_rwork_arr
    real(dp), allocatable :: rwork_arr(:)
    integer :: error_flag, n

    n = 3*blacs%n_atoms
    if (present(src)) A = src
    if (mode(vals_only) == 'V') then
        allocate (vectors(size(A, 1), size(A, 2)))
    else
        allocate (vectors(1, 1))
    end if
    allocate (work_arr(1), rwork_arr(1))
    call PZHEEV( &
        mode(vals_only), 'U', n, A, 1, 1, blacs%desc, eigs, vectors, &
        1, 1, blacs%desc, work_arr , -1, rwork_arr, -1, error_flag &
    )
    n_work_arr = nint(dble(work_arr(1)))
    n_rwork_arr = nint(rwork_arr(1))
    deallocate (work_arr, rwork_arr)
    if (mode(vals_only) == 'N') then
        n_rwork_arr = max(2*n, n_rwork_arr)
    else
        n_rwork_arr = max(4*n - 2, n_rwork_arr)
    end if
    allocate (work_arr(n_work_arr), source=(0d0, 0d0))
    allocate (rwork_arr(n_rwork_arr), source=0d0)
    call PZHEEV( &
        mode(vals_only), 'U', n, A, 1, 1, blacs%desc, eigs, vectors, &
        1, 1, blacs%desc, work_arr, n_work_arr, rwork_arr, n_rwork_arr, &
        error_flag &
    )
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'PZHEEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    if (mode(vals_only) == 'V') A = vectors
end subroutine

function peigvalsh_real(A, blacs, exc, destroy) result(eigs)
    real(dp), target, intent(in) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigs(3*blacs%n_atoms)

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
    call peigh_real(A_p, blacs, eigs, exc, vals_only=.true.)
end function

function peigvalsh_complex(A, blacs, exc, destroy) result(eigs)
    complex(dp), target, intent(in) :: A(:, :)
    type(blacs_desc_t), intent(in) :: blacs
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigs(3*blacs%n_atoms)

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
    call peigh_complex(A_p, blacs, eigs, exc, vals_only=.true.)
end function

end module
