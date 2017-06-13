!\file matrix_utils.f90
!\author M. Sadhukhan
!\brief provides non-standard matrix operations
module matrix_util
contains
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
integer :: i
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

integer :: n, info, i
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
function repeat(n, number)result(vector)
implicit none
integer, intent(in) :: n
double precision,intent(in) :: number
double precision, dimension(n) :: vector
vector = number
end function repeat
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
end module matrix_util
