module mbd_math

use mbd, only: invert, pi

implicit none

integer, parameter :: n_pts_coulomb = 50
real(8), parameter :: L_coulomb = 10.d0

contains

subroutine calc_coulomb_coupled_gauss(R1, R2, K, dip, coul)
    real(8), intent(in) :: R1(3), R2(3), K(6, 6)
    real(8), intent(out), optional :: dip(3, 3), coul

    real(8), dimension(n_pts_coulomb) :: u, w, x
    real(8) :: R(6), s
    integer :: i
    real(8) :: det_K_plus_U2, coul_u, dot, dist
    real(8) :: work(6, 6), Ks(6, 6)
    real(8), dimension(3, 3) :: K11, K12, K22, dip_u

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
        call invert(work)  ! work is (K+U2)^-1
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
        real(8), intent(inout) :: A(6, 6)
        real(8), intent(in) :: u_sq

        integer :: i

        forall (i = 1:3)
            A(i, i) = A(i, i) + u_sq
            A(i, i+3) = A(i, i+3) - u_sq
            A(i+3, i) = A(i+3, i) - u_sq
            A(i+3, i+3) = A(i+3, i+3) + u_sq
        end forall
    end subroutine

end subroutine


real(8) function get_det(A) result(D)
    real(8), intent(in) :: A(:, :)

    integer :: n, i, info
    real(8) :: LU(size(A, 1), size(A, 1))
    integer :: ipiv(size(A, 1))

    n = size(A, 1)
    LU = A
    call DGETRF(n, n, LU, n, ipiv, info)
    D = product((/ (LU(i, i), i = 1, n) /))
end function

function get_outer(a, b) result(C)
    real(8), intent(in) :: a(:), b(:)
    real(8) :: C(size(a), size(b))

    integer :: i, j

    forall (i = 1:size(a), j = 1:size(b)) C(i, j) = a(i)*b(j)
end function

subroutine print_matrix(label, A)
    character(len=*), intent(in) :: label
    real(8), intent(in) :: A(:, :)

    integer :: m, n, i, j

    m = size(A, 1)
    n = size(A, 2)

    write (6, '(A,":")') label
    do i = 1, m
        do j = 1, n
            write (6, "(g10.3)", advance="no") A(i, j)
        end do
        write (6, *)
    end do
end subroutine

end module
