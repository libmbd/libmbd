! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DO_COMPLEX_TYPE

module mbd_elsi

use elsi, only: elsi_handle, elsi_init, elsi_set_mpi, elsi_set_blacs, &
    elsi_set_unit_ovlp, elsi_finalize, elsi_ev_real, elsi_ev_complex

use mbd_constants
use mbd_blacs, only: blacs_desc_t
use mbd_utils, only: exception_t, is_true

implicit none

private
public :: elsi_eigh, elsi_eigvalsh

interface elsi_eigh
    module procedure elsi_eigh_real
    module procedure elsi_eigh_complex
end interface

interface elsi_eigvalsh
    module procedure elsi_eigvalsh_real
    module procedure elsi_eigvalsh_complex
end interface

contains

#endif

#ifndef DO_COMPLEX_TYPE
subroutine elsi_eigh_real(A, blacs_desc, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    real(dp), intent(in), optional :: src(:, :)
#else
subroutine elsi_eigh_complex(A, blacs_desc, eigs, exc, src, vals_only)
    complex(dp), intent(inout) :: A(:, :)
    complex(dp), intent(in), optional :: src(:, :)
#endif
    type(blacs_desc_t), intent(in) :: blacs_desc
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: vals_only

    integer :: n_vecs, n
    type(elsi_handle) :: handle
#ifndef DO_COMPLEX_TYPE
    real(dp) :: DUMMY_MATRIX(1, 1)
    real(dp), allocatable :: vecs(:, :)
#else
    complex(dp) :: DUMMY_MATRIX(1, 1)
    complex(dp), allocatable :: vecs(:, :)
#endif

    n = 3*blacs_desc%n_atoms
    if (present(src)) A = src
    if (is_true(vals_only)) then
        n_vecs = 0
    else
        n_vecs = n
    end if
    call elsi_init(handle, 1, 1, 0, n, 0d0, n_vecs)
    call elsi_set_mpi(handle, blacs_desc%comm)
    call elsi_set_blacs(handle, blacs_desc%ctx, blacs_desc%blocksize)
    call elsi_set_unit_ovlp(handle, 1)
    allocate (vecs(size(A, 1), size(A, 2)))
#ifndef DO_COMPLEX_TYPE
    call elsi_ev_real(handle, A, DUMMY_MATRIX, eigs, vecs)
#else
    call elsi_ev_complex(handle, A, DUMMY_MATRIX, eigs, vecs)
#endif
    A = vecs
    call elsi_finalize(handle)
end subroutine

#ifndef DO_COMPLEX_TYPE
function elsi_eigvalsh_real(A, blacs_desc, exc, destroy) result(eigs)
    real(dp), allocatable, intent(inout) :: A(:, :)
#else
function elsi_eigvalsh_complex(A, blacs_desc, exc, destroy) result(eigs)
    complex(dp), allocatable, intent(inout) :: A(:, :)
#endif
    type(blacs_desc_t), intent(in) :: blacs_desc
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigs(3*blacs_desc%n_atoms)

#ifndef DO_COMPLEX_TYPE
    real(dp), allocatable :: A_(:, :)
#else
    complex(dp), allocatable :: A_(:, :)
#endif

    if (is_true(destroy)) then
        call move_alloc(A, A_)
    else
        A_ = A
    end if
    call elsi_eigh(A_, blacs_desc, eigs, exc, vals_only=.true.)
end function

#ifndef DO_COMPLEX_TYPE
#   define DO_COMPLEX_TYPE
#   include "mbd_elsi.F90"

end module

#endif
