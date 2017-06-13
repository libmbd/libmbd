!\file integralPack.f90
!\author Mainak Sadhukhan
!\brief a collection of numerical quadrature codes

module quadratures
use constants, only: pi
use ofunction

contains
!subroutine integration_gL(order, farray, ifzero2inf, res)
!!the Gauss-Legendre quadrature routine
!implicit none
!integer, intent(in) :: order
!double precision, dimension(order, 2), intent(in) :: farray
!logical, intent(in) :: ifzero2inf
!double precision, intent(out) :: res
!integer :: n, i
!double precision, dimension(order):: x, w
!
!res= 10.d0
!
!if(ifzero2inf .eqv. .true.)then
! print*, "Zero - Infty"
! do i = 1, order
!print*, farray(i, :)
! end do
!
!else
! print*, "Ordinary"
!end if
!
!
!end subroutine integration_gL
!!
subroutine gl_points(n, x, w)
implicit none
! provides Gauss-Legendre abcissa and weights
integer, parameter ::  m = 10
integer, intent(in) :: n
double precision, dimension(n), intent(out) :: x, w
integer :: i
do i = 1, n
 x(i) = root_Legendre(n, i)
 w(i) = 2.d0/((1.d0-x(i)**2)*(derivlegendre(n,x(i)))**2)
end do
end subroutine gl_points
!
function root_Legendre(n, k)result(res)
! provides the kth root for P_n(n)
implicit none
double precision, parameter :: lim_err=1.d-12
integer, intent(in) :: n, k
double precision :: err, res,root_guess_l,root_guess_u, x

err = 1.d0
! guess from A&S 22.16.6
root_guess_l= cos(pi()*(2*k-1)/(2.d0*n+1))
root_guess_u= cos((pi()*2*k)/(2.d0*n+1))
! Newton-Raphson
x =root_guess_u
!print*, "Guess", x
do while (err >= dabs(lim_err))
!print*, "P = ",legendre(n, x),"P' = ",derivlegendre(n, x)
x = x - legendre(n, x)/derivlegendre(n, x)
!Print*, "err = ", abs(legendre(n, x))
err =abs(legendre(n, x))
end do
!
res = x
end function root_Legendre
end module quadratures




