! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_coulomb

use mbd_constants
use mbd_linalg, only: eye, outer, diag
use mbd_lapack, only: inverse, det, inv
use mbd_system_type, only: mbd_system
use mbd_dipole, only: dipole_matrix
use mbd_damping_type, only: mbd_damping, damping_fermi
use mbd_matrix_type, only: mbd_matrix_real

implicit none

private

public :: coulomb_energy, dipole_energy

integer :: n_pts_coulomb = 500
real(dp) :: L_coulomb = 10d0
character(len=20) :: quadrature = 'simpson'

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
            det_K_plus_U2 = det(K+eye(3)*u(i)**2)
        case (6)
            work = K
            call add_U2(work, u(i)**2)  ! work is K+U2
            det_K_plus_U2 = det(work)
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
            dip_u = (-2*K12+4*outer( &
                matmul(K11, R1)+matmul(K12, R2), &
                matmul(K12, R1)+matmul(K22, R2) &
            ))*coul_u
            dip = dip + dip_u
        end if
    end do
    det_K = det(K)
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

real(dp) function coulomb_energy(sys, q, m, w_t, C, damp)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: q(:), m(:), w_t(:), C(:, :)
    type(mbd_damping), intent(in) :: damp

    real(dp), allocatable :: O(:, :)
    real(dp) :: OAB(6, 6), OABm(6, 6), RA(3), RB(3), ene_ABi(4), prod_w_t, &
        K(3, 3), s_vdw, f_damp
    integer, allocatable :: notAB(:)
    integer :: N, A, B, i, j, AB(6), i2A(6)

    allocate (notAB(size(C, 1)-6))
    O = matmul(matmul(C, diag(w_t)), transpose(C))
    N = sys%siz()
    prod_w_t = product(w_t)
    coulomb_energy = 0.d0
    do A = 1, N
        do B = A+1, N
            RA = sys%coords(:, A)
            RB = sys%coords(:, B)
            AB(:) = [(3*(A-1)+i, i = 1, 3), (3*(B-1)+i, i = 1, 3)]
            notAB(:) = [ &
                (i, i = 1, 3*(A-1)), &
                (i, i = 3*A+1, 3*(B-1)), &
                (i, i = 3*B+1, 3*N) &
            ]
            OAB = O(AB, AB) - matmul( &
                O(AB, notAB), &
                matmul(inverse(O(notAB, notAB)), O(notAB, AB)) &
            )
            i2A = [(A, i = 1, 3), (B, i = 1, 3)]
            forall (i = 1:6, j = 1:6)
                OABm(i, j) = OAB(i, j)*sqrt(m(i2A(i))*m(i2A(j)))
            end forall
            call calc_coulomb_coupled_gauss(RA, RB, OABm, coul=ene_ABi(1))
            K = m(B)*(OAB(4:6, 4:6) - matmul( &
                OAB(4:6, 1:3), matmul(inverse(OAB(1:3, 1:3)), OAB(1:3, 4:6)) &
            ))
            call calc_coulomb_coupled_gauss(RA, RB, K, coul=ene_ABi(2))
            K = m(A)*(OAB(1:3, 1:3) - matmul( &
                OAB(1:3, 4:6), matmul(inverse(OAB(4:6, 4:6)), OAB(4:6, 1:3)) &
            ))
            call calc_coulomb_coupled_gauss(RA, RB, K, coul=ene_ABi(3))
            ene_ABi(2:3) = -ene_ABi(2:3)
            ene_ABi(4) = 1d0/sqrt(sum((RA-RB)**2))
            select case (damp%version)
            case ('fermi')
                s_vdw = damp%ts_sr*sum(damp%r_vdw([A, B]))
                f_damp = damping_fermi(RA-RB, s_vdw, damp%ts_d)
            case default
                f_damp = 1d0
            end select
            coulomb_energy = coulomb_energy + f_damp*q(A)*q(B)*sum(ene_ABi)
        end do
    end do
end function

real(dp) function dipole_energy(sys, a0, w, w_t, C, damp)
    type(mbd_system), intent(inout) :: sys
    real(dp), intent(in) :: a0(:), w(:), w_t(:), C(:, :)
    type(mbd_damping), intent(in) :: damp

    integer :: A, B, i, j, N
    type(mbd_matrix_real) :: T

    N = sys%siz()
    T = dipole_matrix(sys, damp)
    do  A = 1, N
        do B = 1, N
            i = 3*(A-1)
            j = 3*(B-1)
            T%val(i+1:i+3, j+1:j+3) = &
                w(A)*w(B)*sqrt(a0(A)*a0(B))*T%val(i+1:i+3, j+1:j+3)
        end do
    end do
    T%val = matmul(matmul(transpose(C), T%val), C)
    dipole_energy = sum(diag(T%val)/(4*w_t))
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
