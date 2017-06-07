!\file io.f90
!\author M. Sadhukhan
!\brief contains custom-designed
module ios
contains
subroutine printmat(A, message)
implicit none
double precision, dimension(:, :), intent(in) :: A
character(len=*), optional :: message
integer :: i, n

if(present(message))print*, message
n = size(A, 1)
do i = 1, n
print*, A(i, :)
end do

end subroutine printmat
end module ios
