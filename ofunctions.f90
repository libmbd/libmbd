!\file ofunctions.f90
!\author M. Sadhukhan
!\brief Provides values for a given order of an orthogoal polynomial

module ofunction
contains
recursive function legendre(n, x)result(res)
! Provides P_n(x)
implicit none
integer, intent(in):: n
double precision, intent(in) :: x
double precision :: res
integer :: i

if (n == 0) then
 res = 1.d0
else if (n == 1) then
 res = x
else
 res = ((2*n-1)*x*legendre(n-1, x)-(n-1)*legendre(n-2, x))/(n*1.d0)
end if
!
end function legendre
! 
recursive function derivlegendre(n, x)result(res)
! Provides first derivative of P_n(x) wrt x
implicit none
integer, intent(in):: n
double precision, intent(in) :: x
double precision :: res
integer :: i

if (n == 0) then
res = 0.d0
else if (n == 1) then
res = 1.d0
else
res = ((2*n-1)*legendre(n-1, x)+ (2*n-1)*x*derivlegendre(n-1, x)-(n-1)*derivlegendre(n-2, x))/(n*1.d0)
end if
!
end function derivlegendre
!
integer recursive function factorial(n)result(res)
implicit none
integer :: n

if (n == 0)then
 res = 1
else
 res= n*factorial(n-1)
end if
end function factorial
!
end module ofunction
