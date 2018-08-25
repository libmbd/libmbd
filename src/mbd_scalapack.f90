! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_scalapack

use mbd_common, only: dp, exception => mbd_exc, MBD_EXC_LINALG, tostr
use mbd_parallel, only: mbd_blacs
use mbd_lapack, only: mode

implicit none

private
public :: pmmul, pinvh, peigh, peigvalsh

interface pmmul
    module procedure pmmul_real
end interface

interface pinvh
    module procedure pinvh_real
end interface

interface peigh
    module procedure peigh_real
end interface

interface peigvalsh
    module procedure peigvalsh_real
end interface

external :: PDSYEV, PDGETRF, PDGETRI, PDGEMM

contains

subroutine pinvh_real(A, blacs, exc, src)
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
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr), iwork_arr(n_iwork_arr))
    call PDGETRI( &
        n, A, 1, 1, blacs%desc, i_pivot, &
        work_arr, n_work_arr, iwork_arr, n_iwork_arr, error_flag)
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
    type(mbd_blacs), intent(in) :: blacsA, blacsB, blacsC
    logical, intent(in), optional :: transA, transB
    real(dp) :: C(size(A, 1), size(B, 2))

    character :: transA_, transB_
    integer :: n

    transA_= 'N'
    transB_ = 'N'
    if (present(transA)) then
        if (transA) transA_ = 'T'
    end if
    if (present(transB)) then
        if (transB) transB_ = 'T'
    end if
    n = 3*blacsA%n_atoms
    call PDGEMM( &
        transA_, transB_, n, n, n, 1d0, A, 1, 1, blacsA%desc, &
        B, 1, 1, blacsB%desc, 0d0, C, 1, 1, blacsC%desc &
    )
end function

subroutine peigh_real(A, blacs, eigs, exc, src, vals_only)
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
        allocate (vectors(size(A, 1), size(A, 2)))
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
            exc%code = MBD_EXC_LINALG
            exc%origin = 'PDSYEV'
            exc%msg = 'Failed with code ' // trim(tostr(error_flag))
        end if
        return
    endif
    if (mode(vals_only) == 'V') A = vectors
end subroutine

function peigvalsh_real(A, blacs, exc, destroy) result(eigs)
    real(dp), target, intent(in) :: A(:, :)
    type(mbd_blacs), intent(in) :: blacs
    type(exception), intent(out), optional :: exc
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

end module
