!\file constant.f90
!\author M. Sadhukhan
!\brief Provides numerical and other constants
module constants
contains
double precision function pi()
pi = 4.d0*datan(1.d0)
!return pi
end function pi
!
end module constants
