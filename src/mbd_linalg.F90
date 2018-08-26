! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_linalg

use mbd_constants

implicit none

private
public :: outer, eye, diag

interface diag
    module procedure get_diag_real
    module procedure get_diag_complex
    module procedure make_diag_real
end interface

contains

function outer(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))

    integer :: i, j

    do i = 1, size(a)
        do j = 1, size(b)
            c(i, j) = a(i)*b(j)
        end do
    end do
end function

function eye(n) result(A)
    integer, intent(in) :: n
    real(dp) :: A(n, n)

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:n) A(i, i) = 1.d0
end function

function get_diag_real(A) result(d)
    real(8), intent(in) :: A(:, :)
    real(8) :: d(size(A, 1))

    integer :: i

    forall (i = 1:size(A, 1)) d(i) = A(i, i)
end function

function get_diag_complex(A) result(d)
    complex(dp), intent(in) :: A(:, :)
    complex(dp) :: d(size(A, 1))

    integer :: i

    forall (i = 1:size(A, 1)) d(i) = A(i, i)
end function

function make_diag_real(d) result(A)
    real(dp), intent(in) :: d(:)
    real(dp) :: A(size(d), size(d))

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:size(d)) A(i, i) = d(i)
end function

end module mbd_linalg
