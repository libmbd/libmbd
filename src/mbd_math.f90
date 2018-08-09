! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_math

use mbd_linalg, only: eye, inverse, diag, inv
use mbd, only: dipole_matrix, mbd_system, mbd_damping, get_sigma_selfint
use mbd_common, only: dp, pi
use mbd_types, only: mat3n3n

implicit none

integer, parameter :: n_pts_coulomb = 50
real(dp), parameter :: L_coulomb = 10.d0
real(dp), parameter :: point_charge = 100.d0

external :: DGETRF

contains

subroutine calc_coulomb_coupled_gauss(R1, R2, K, dip, coul)
    real(dp), intent(in) :: R1(3), R2(3), K(6, 6)
    real(dp), intent(out), optional :: dip(3, 3), coul

    real(dp), dimension(n_pts_coulomb) :: u, w, x
    real(dp) :: R(6), s
    integer :: i
    real(dp) :: det_K_plus_U2, coul_u, dot, dist
    real(dp) :: work(6, 6), Ks(6, 6)
    real(dp), dimension(3, 3) :: K11, K12, K22, dip_u

    ! print *, "det(K)", get_det(K)
    s = get_det(K)**(-1.d0/6)
    Ks = s*K
    w(:) = 1.d0/n_pts_coulomb
    forall (i = 1:n_pts_coulomb) x(i) = w(1)/2+(i-1)*w(1)
    ! u = x/(1.d0-x)
    ! w = 1.d0/(1.d0-x)**2*w
    dist = sqrt(sum((R1-R2)**2))
    u = log(1.d0/(1.d0-x))*sqrt(s)*L_coulomb/dist
    w = 1.d0/(1.d0-x)*w*sqrt(s)*L_coulomb/dist
    R(1:3) = R1/sqrt(s)
    R(4:6) = R2/sqrt(s)
    if (present(coul)) coul = 0.d0
    if (present(dip)) dip(:, :) = 0.d0
    do i = 1, n_pts_coulomb
        work = Ks
        call add_U2(work, u(i)**2)  ! work is K+U2
        det_K_plus_U2 = get_det(work)
        ! call print_matrix('K', K)
        ! call print_matrix('K+U2', work)
        call inv(work)  ! work is (K+U2)^-1
        ! call print_matrix('(K+U2)^-1', work)
        work = Ks-matmul(Ks, matmul(work, Ks)) ! work is K-K*(K+U2)^-1*K
        ! call print_matrix('K-K*(K+U2)^-1*K', work)
        dot = dot_product(R, matmul(work, R))
        coul_u = 1.d0/sqrt(det_K_plus_U2)*exp(-dot)*w(i)
        if (present(coul)) coul = coul + coul_u
        if (present(dip)) then
            K11 = work(1:3, 1:3)
            K12 = work(1:3, 4:6)
            K22 = work(4:6, 4:6)
            dip_u = (-2*K12+4*get_outer( &
                matmul(K11, R1)+matmul(K12, R2), &
                matmul(K12, R1)+matmul(K22, R2) &
            )/s)*coul_u
            dip = dip + dip_u
        end if
        ! print *, "u =", u(i)**2, "w =", w(i), "1/sqrt(det(K+U2)) =", 1.d0/sqrt(det_K_plus_U2), &
        !     "dot =", dot, "exp =", exp(-dot), "add =", coul_u
    end do
    if (present(coul)) coul = 2.d0/sqrt(pi)*coul/sqrt(s)
    if (present(dip)) dip = 2.d0/sqrt(pi)*dip*s**(-3.d0/2)

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

real(dp) function get_coulomb_energy_coupled_osc(R, q, m, w_t, C) result(ene)
    real(dp), intent(in) :: R(:, :), q(size(R, 1)), m(size(R, 1)), w_t(3*size(R, 1))
    real(dp), intent(in) :: C(3*size(R, 1), 3*size(R, 1))

    real(dp) :: O(size(C, 1), size(C, 1))
    real(dp) :: Opp(size(C, 1)-6, size(C, 1)-6)
    real(dp) :: OAB(6, 6), OABm(6, 6), K(6, 6)
    real(dp) :: RA(3), RB(3)
    integer :: N, A, B, i, j
    integer :: AB(6), notAB(size(C, 1)-6)
    real(dp) :: ene_AB, ene_ABi(4)
    real(dp) :: prod_w_t, coul
    integer :: i2A(6) = (/ 1, 1, 1, 2, 2, 2 /)

    O = matmul(matmul(C, diag(w_t)), transpose(C))
    N = size(R, 1)
    prod_w_t = product(w_t)
    ene = 0.d0
    do A = 1, N
        do B = A+1, N
            RA = R(A, :)
            RB = R(B, :)
            AB(:) = (/ (3*(A-1)+i, i = 1, 3),  (3*(B-1)+i, i = 1, 3) /)
            notAB(:) = (/ (i, i = 1, 3*(A-1)),  (i, i = 3*A+1, 3*(B-1)), (i, i = 3*B+1, 3*N) /)
            Opp(:, :) = O(notAB, notAB)
            OAB = O(AB, AB)-matmul(O(AB, notAB), matmul(inverse(O(notAB, notAB)), O(notAB, AB)))
            forall (i = 1:6, j = 1:6) OABm(i, j) = OAB(i, j)*sqrt(m(i2A(i))*m(i2A(j)))
            call calc_coulomb_coupled_gauss(RA, RB, OABm, coul=coul)
            ene_ABi(1) = 1.d0/sqrt(get_det(OAB))*coul
            K(:, :) = 0.d0
            K(1:3, 1:3) = point_charge*eye(3)
            K(4:6, 4:6) = m(B)*(OAB(4:6, 4:6)-matmul(OAB(4:6, 1:3), matmul(inverse(OAB(1:3, 1:3)), OAB(1:3, 4:6))))
            call calc_coulomb_coupled_gauss(RA, RB, K, coul=coul)
            ene_ABi(2) = -1.d0/sqrt(get_det(OAB(1:3, 1:3))*get_det(K(4:6, 4:6))/m(B)**3)*coul
            K(:, :) = 0.d0
            K(1:3, 1:3) = m(A)*(OAB(1:3, 1:3)-matmul(OAB(1:3, 4:6), matmul(inverse(OAB(4:6, 4:6)), OAB(4:6, 1:3))))
            K(4:6, 4:6) = point_charge*eye(3)
            call calc_coulomb_coupled_gauss(RA, RB, K, coul=coul)
            ene_ABi(3) = -1.d0/sqrt(get_det(OAB(4:6, 4:6))*get_det(K(1:3, 1:3))/m(A)**3)*coul
            K(:, :) = 0.d0
            K(1:3, 1:3) = point_charge*eye(3)
            K(4:6, 4:6) = point_charge*eye(3)
            call calc_coulomb_coupled_gauss(RA, RB, K, coul=coul)
            ene_ABi(4) = 1.d0/sqrt(get_det(OAB))*coul
            ene_AB = q(A)*q(B)*sqrt(prod_w_t)*sum(ene_ABi)
            ! print *, q(A)*q(B)*sqrt(prod_w_t)*ene_ABi
            ene = ene + ene_AB
        end do
    end do
end function

real(dp) function get_dipole_energy_coupled_osc(sys, a0, w, w_t, C) result(ene)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: a0(size(sys%coords, 1)), w(size(sys%coords, 1)), w_t(3*size(sys%coords, 1))
    real(dp), intent(in) :: C(3*size(sys%coords, 1), 3*size(sys%coords, 1))

    integer :: A, B, i, j, N
    type(mat3n3n) :: T
    type(mbd_damping) :: damp
    real(dp), allocatable :: dummy(:)

    N = size(sys%coords, 1)
    damp%version = 'dip,gg'
    damp%sigma = get_sigma_selfint(a0, dummy)
    T = dipole_matrix(sys, damp, grad=.false.)
    do  A = 1, N
        do B = 1, N
            i = 3*(A-1)
            j = 3*(B-1)
            T%re(i+1:i+3, j+1:j+3) = w(A)*w(B)*sqrt(a0(A)*a0(B))*T%re(i+1:i+3, j+1:j+3)
        end do
    end do
    ! call print_matrix('T', T)
    T%re = matmul(matmul(transpose(C), T%re), C)
    ! call print_matrix('T_t', T)
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

subroutine swap_ints(a, b)
    integer, intent(inout) :: a, b

    integer :: c

    c = a
    a = b
    b = c
end subroutine

function get_outer(a, b) result(C)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: C(size(a), size(b))

    integer :: i, j

    forall (i = 1:size(a), j = 1:size(b)) C(i, j) = a(i)*b(j)
end function

end module
