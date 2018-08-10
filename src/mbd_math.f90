! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_math

use mbd_linalg, only: eye, inverse, diag, inv, eye
use mbd, only: dipole_matrix, mbd_system, mbd_damping, get_sigma_selfint
use mbd_common, only: dp, pi
use mbd_types, only: mat3n3n

implicit none

integer, parameter :: n_pts_coulomb = 500
real(dp), parameter :: L_coulomb = 10.d0
character(len=20) :: quadrature = 'simpson'

external :: DGETRF

contains

subroutine calc_coulomb_coupled_gauss(R1, R2, K, dip, coul)
    real(dp), intent(in) :: R1(3), R2(3), K(:, :)
    real(dp), intent(out), optional :: dip(3, 3), coul

    real(dp), allocatable :: u(:), w(:), x(:)
    real(dp) :: R(6), det_K
    integer :: i
    real(dp) :: det_K_plus_U2, coul_u, dot, dist, work(6, 6)
    real(dp), dimension(:, :), allocatable :: K11, K12, K22, dip_u

    dist = sqrt(sum((R1-R2)**2))
    allocate (x(n_pts_coulomb), w(n_pts_coulomb))
    select case (quadrature)
    case ('original')
        w = 1d0/n_pts_coulomb
        forall (i = 1:n_pts_coulomb) x(i) = w(1)/2+(i-1)*w(1)
        u = log(1d0/(1d0-x))*L_coulomb/dist
        w = 1d0/(1d0-x)*w*L_coulomb/dist
    case ('simpson')
        call simpson1by3(n_pts_coulomb, x, w)
        u = (1d0-x)/(1d0+x)
        w = 2*w/(1d0+x)**2
    end select
    R(1:3) = R1
    R(4:6) = R2
    if (present(coul)) coul = 0d0
    if (present(dip)) then
        dip(:, :) = 0d0
        allocate (K11(3, 3), K12(3, 3), K22(3, 3), dip_u(3, 3))
    end if
    do i = 1, n_pts_coulomb
        select case (size(K, 1))
        case (3)
            work(1:3, 1:3) = K
            work(4:6, 1:3) = eye(3)*u(i)**2
            work = -matmul( &
                work(:, 1:3), &
                matmul(inverse(K+eye(3)*u(i)**2), transpose(work(:, 1:3))) &
            )
            work(1:3, 1:3) = work(1:3, 1:3) + K
            work(4:6, 4:6) = work(4:6, 4:6) + eye(3)*u(i)**2
            det_K_plus_U2 = get_det(K+eye(3)*u(i)**2)
        case (6)
            work = K
            call add_U2(work, u(i)**2)  ! work is K+U2
            det_K_plus_U2 = get_det(work)
            call inv(work)  ! work is (K+U2)^-1
            work = K-matmul(K, matmul(work, K)) ! work is K-K*(K+U2)^-1*K
        end select
        dot = dot_product(R, matmul(work, R))
        coul_u = 1d0/sqrt(det_K_plus_U2)*exp(-dot)*w(i)
        if (present(coul)) coul = coul + coul_u
        if (present(dip)) then
            K11 = work(1:3, 1:3)
            K12 = work(1:3, 4:6)
            K22 = work(4:6, 4:6)
            dip_u = (-2*K12+4*get_outer( &
                matmul(K11, R1)+matmul(K12, R2), &
                matmul(K12, R1)+matmul(K22, R2) &
            ))*coul_u
            dip = dip + dip_u
        end if
    end do
    det_K = get_det(K)
    if (present(coul)) coul = 2.d0/sqrt(pi)*coul*sqrt(det_K)
    if (present(dip)) dip = 2.d0/sqrt(pi)*dip*sqrt(det_K)

    contains

    subroutine  add_U2(A, u_sq)
        real(dp), intent(inout) :: A(6, 6)
        real(dp), intent(in) :: u_sq

        integer :: i

        forall (i = 1:3)
            A(i, i) = A(i, i) + u_sq
            A(i, i+3) = A(i, i+3) - u_sq
            A(i+3, i) = A(i+3, i) - u_sq
            A(i+3, i+3) = A(i+3, i+3) + u_sq
        end forall
    end subroutine

end subroutine

real(dp) function get_coulomb_energy_coupled_osc(sys, q, m, w_t, C) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: q(:), m(:), w_t(:), C(:, :)

    real(dp) :: O(size(C, 1), size(C, 1))
    real(dp) :: Opp(size(C, 1)-6, size(C, 1)-6)
    real(dp) :: OAB(6, 6), OABm(6, 6), K(6, 6)
    real(dp) :: RA(3), RB(3)
    integer :: N, A, B, i, j
    integer :: AB(6), notAB(size(C, 1)-6)
    real(dp) :: ene_AB, ene_ABi(4)
    real(dp) :: prod_w_t, coul
    integer :: i2A(6)

    O = matmul(matmul(C, diag(w_t)), transpose(C))
    N = sys%siz()
    prod_w_t = product(w_t)
    ene = 0.d0
    do A = 1, N
        do B = A+1, N
            RA = sys%coords(:, A)
            RB = sys%coords(:, B)
            AB(:) = (/ (3*(A-1)+i, i = 1, 3),  (3*(B-1)+i, i = 1, 3) /)
            notAB(:) = (/ (i, i = 1, 3*(A-1)),  (i, i = 3*A+1, 3*(B-1)), (i, i = 3*B+1, 3*N) /)
            Opp(:, :) = O(notAB, notAB)
            OAB = O(AB, AB)-matmul(O(AB, notAB), matmul(inverse(O(notAB, notAB)), O(notAB, AB)))
            i2A = [(A, i = 1, 3), (B, i = 1, 3)]
            forall (i = 1:6, j = 1:6) OABm(i, j) = OAB(i, j)*sqrt(m(i2A(i))*m(i2A(j)))
            call calc_coulomb_coupled_gauss(RA, RB, OABm, coul=coul)
            ene_ABi(1) = coul
            K(1:3, 1:3) = m(B)*(OAB(4:6, 4:6)-matmul(OAB(4:6, 1:3), matmul(inverse(OAB(1:3, 1:3)), OAB(1:3, 4:6))))
            call calc_coulomb_coupled_gauss(RA, RB, K(1:3, 1:3), coul=coul)
            ene_ABi(2) = -coul
            K(1:3, 1:3) = m(A)*(OAB(1:3, 1:3)-matmul(OAB(1:3, 4:6), matmul(inverse(OAB(4:6, 4:6)), OAB(4:6, 1:3))))
            call calc_coulomb_coupled_gauss(RA, RB, K(1:3, 1:3), coul=coul)
            ene_ABi(3) = -coul
            ene_ABi(4) = 1d0/sqrt(sum((RA-RB)**2))
            ene_AB = q(A)*q(B)*sum(ene_ABi)
            ene = ene + ene_AB
        end do
    end do
end function

real(dp) function get_dipole_energy_coupled_osc(sys, a0, w, w_t, C, damp) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: a0(:), w(:), w_t(:), C(:, :)
    type(mbd_damping), intent(in) :: damp

    integer :: A, B, i, j, N
    type(mat3n3n) :: T

    N = sys%siz()
    T = dipole_matrix(sys, damp, grad=.false.)
    do  A = 1, N
        do B = 1, N
            i = 3*(A-1)
            j = 3*(B-1)
            T%re(i+1:i+3, j+1:j+3) = &
                w(A)*w(B)*sqrt(a0(A)*a0(B))*T%re(i+1:i+3, j+1:j+3)
        end do
    end do
    T%re = matmul(matmul(transpose(C), T%re), C)
    ene = sum(diag(T%re)/(4*w_t))
end function

real(dp) function get_det(A) result(D)
    real(dp), intent(in) :: A(:, :)

    integer :: n, i, info
    real(dp) :: LU(size(A, 1), size(A, 1))
    integer :: ipiv(size(A, 1))

    n = size(A, 1)
    LU = A
    call DGETRF(n, n, LU, n, ipiv, info)
    D = product((/ (LU(i, i), i = 1, n) /))
end function

function get_outer(a, b) result(C)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: C(size(a), size(b))

    integer :: i, j

    forall (i = 1:size(a), j = 1:size(b)) C(i, j) = a(i)*b(j)
end function

subroutine simpson1by3(n, x, w)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n), w(n)

    integer :: i
    real(dp) :: h, delta

    delta = 1d-6
    h = 2*(1d0-delta)/(n-1)
    do i = 1, n
        x(i) = -(1d0-delta)+h*(i-1)
        if (2*(i-(i/2)) == i) then
            w(i) = 2*h/3
        else
            w(i) = 4*h/3
        end if
    end do
    w(1) = h/3
    w(n)= w(1)
end subroutine

end module
