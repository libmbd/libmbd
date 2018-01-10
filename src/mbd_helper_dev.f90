!\file mbd_helper_dev.f90
!\author M. Sadhukhan
module mbd_helper_dev

use mbd_common, only: pi

external :: DGETRI, DGETRF

contains

recursive function legendre(n, x)result(res)
! Provides P_n(x)
implicit none
integer, intent(in):: n
double precision, intent(in) :: x
double precision :: res

if (n == 0) then
 res = 1.d0
else if (n == 1) then
 res = x
else
 res = ((2.d0*n-1.d0)*x*legendre(n-1, x)-(n-1)*legendre(n-2, x))/(n*1.d0)
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

if (n == 0) then
res = 0.d0
else if (n == 1) then
res = 1.d0
else
res = ((2.d0*n-1.d0)*legendre(n-1, x)+ (2.d0*n-1.d0)*x*derivlegendre(n-1, x)-(n-1)*derivlegendre(n-2, x))/(n*1.d0)
end if
!
end function derivlegendre
!
subroutine gl_points(n, x, w)
implicit none
! provides Gauss-Legendre abcissa and weights
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
root_guess_l= cos(pi*(2*k-1)/(2.d0*n+1))
root_guess_u= cos((pi*2*k)/(2.d0*n+1))
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
integer recursive function factorial(n)result(res)
implicit none
integer :: n

if (n == 0)then
 res = 1
else
 res= n*factorial(n-1)
end if
end function factorial

subroutine diag_double(vector, matrix)
! provides matrix with the diagonal containing vector
! suitable for double precision
implicit none
integer :: i
double precision, dimension(:), intent(in) :: vector
double precision, dimension(size(vector), size(vector)), intent(out):: matrix
matrix = 0.d0

do i = 1, size(vector)
 matrix(i,i) = vector(i)
end do
end subroutine diag_double
!
subroutine det(A, dt)
! provides determinant dt of matrix A
implicit none
integer :: n, info, i
double precision, dimension(:, :) :: A
double precision :: dt
integer, dimension(size(A, 1)) :: ipiv
dt = 1.d0
n = size(A, 1)
!do i = 1, n
!print*, A(i, :)
!end do
call dgetrf(n, n, A, n, ipiv, info)

!print*, "ipiv :", ipiv

if (info == 0)then
 do i = 1, n
   if(ipiv(i) == i)then
    dt = dt*A(i, i)
   else
    dt = -dt*A(i, i)
   end if
 end do
else
 print*, "dgertf did not exited sucessfully"
 print*, "info =", info
 stop
end if
!print*, "dt = ",dt
end subroutine det
!
function i_matrix(n) result(matrix)
! provides identity matrix of dimension n
implicit none
integer, intent(in) :: n
double precision, dimension(n, n) :: matrix
double precision, dimension(n) :: diag1
diag1 = 1.d0
call diag_double(diag1, matrix)
end function i_matrix
!
function matrix_combine( n, matA, matB)result(matAB)
implicit none
! combine two nxn matrices
integer, intent(in) ::  n
double precision, dimension(n, n), intent(in) :: matA, matB
double precision, dimension(2*n, n) :: matAB

matAB(1:n, :) = matA
matAB(n+1:2*n, : ) = matB
end function matrix_combine
!
function matrix_embed(n, matA, Bign, rpos, cpos)result(Bigmatrix)
implicit none
integer, intent(in) :: n, Bign, rpos, cpos
double precision, dimension(n, n), intent(in) :: matA
double precision, dimension(Bign, Bign) :: Bigmatrix
if (rpos+n-1 > Bign .or. cpos+n-1 > Bign) then
 print*, "Matrix sizes are not fit for embedding ... aborting execusion"
 stop
end if
Bigmatrix = 0.d0
Bigmatrix(rpos:, cpos:) =matA

end function matrix_embed
!
function Matrix_invert(A)result(Ainv)
! Provides matrix inversion of A
implicit none
double precision, dimension(:,:), intent(in) :: A
double precision, dimension(size(A, 1),size(A, 1)):: Ainv

integer :: n, info
integer, dimension(size(A,1)) :: ipiv
double precision, dimension(size(A,1)**2) :: work
ipiv= 1
work = 0.d0
n = size(A, 1)
Ainv = 0.d0
if (n>1)then
call DGETRI(n, A, n, ipiv, work, n**2, info)
if(info /=0)then
 print*, "Matrix inversion failed ... Aborting, info =", info
 stop
end if
 Ainv= A
else if (n == 1)then
 Ainv= A
end if

end function Matrix_invert
!
function my_repeat(n, number)result(vector)
implicit none
integer, intent(in) :: n
double precision,intent(in) :: number
double precision, dimension(n) :: vector
vector = number
end function my_repeat
!
function combine_vector(nA, vA, nB, vB)result(vAB)
implicit none
integer, intent(in) ::  nA, nB
double precision, dimension(na), intent(in) :: vA
double precision, dimension(nb), intent(in) :: vB
double precision, dimension(na+nb) :: vAB
vAB = 0.d0
vAB(1:na) = vA
vAB(na+1:nb+na) = vB
end function combine_vector


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

end module
