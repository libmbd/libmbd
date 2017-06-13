!\file io.f90
!\author M. Sadhukhan
!\brief contains custom-designed
module ios
contains
subroutine printmat(A, message)
implicit none
double precision, dimension(:, :), intent(in) :: A
character(len=*), optional :: message
character(len=20) :: fmt
integer :: i, n

if(present(message))print*, message
n = size(A, 1)
write (fmt, *) n
fmt = "(1x,"//trim(adjustl(fmt))//"f8.4)"
!print*, fmt
!stop
do i = 1, n
 write(*, trim(fmt)) A(i, :)
end do

end subroutine printmat
end module ios
