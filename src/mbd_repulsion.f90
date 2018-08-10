module mbd_repulsion
use mbd_linalg, only: inverse
use mbd_math, only: get_det
use mbd_types, only: scalar, mat3n3n
use mbd, only: dipole_matrix, damping_fermi, mbd_system, mbd_damping
use mbd_common, only: pi, dp

external :: DGETRI, DGETRF

contains
!
subroutine E2body(largeomega, RA, RB, ee)
implicit none
double precision, dimension(6,6),intent(in) :: largeomega
double precision,dimension(3),intent(in) :: RA, RB
double precision, intent(out) :: ee
integer :: i, nint
double precision, dimension(:), allocatable ::x, w
double precision, dimension(6,6) :: BigU2, intermedom
double precision, dimension(6) :: RAB
double precision :: den, u, exponent
!
ee = 0.d0
if(allocated(x))deallocate(x)
if(allocated(w))deallocate(w)
rab(1:3)= ra
rab(4:6)= rb
nint = 500
allocate(x(nint),w(nint))
!call gl_points(nint, x, w)
call simpson1by3(nint, x, w)

do i = 1, nint
 u = (1.d0-x(i))/(1.d0+x(i))
 call U2(u, BigU2)
 den = 1.d0/dsqrt(get_det(largeomega+BigU2))/(1.d0+x(i))**2.d0
 intermedom=largeomega-matmul(largeomega,matmul(inverse(largeomega+BigU2), largeomega))
 exponent =dot_product(rab, matmul(intermedom,rab))
 !print*, "u", u, "denometer", den, "exponent", exponent, w(i), exp(-exponent), get_det(largeomega+BigU2),&
 !&2.d0*w(i)*exp(-exponent)*den
  ee = ee + 2.d0*w(i)*exp(-exponent)*den
end do

! do i = 1, nint
!   h=1.d0/(nint-1)
!   xx= -1.d0+h*i
!   u = (1.d0-xx)/(1.d0+xx)
!   call U2(u, BigU2)
!   den = 1.d0/dsqrt(get_det(largeomega+BigU2))/(1.d0+xx)**2.d0
!   intermedom=largeomega-matmul(largeomega,matmul(inverse(largeomega+BigU2), largeomega))
!   exponent =dot_product(rab, matmul(intermedom,rab))
!   ! ! print*, 'den =', den, "exponential", exp(-exponent)
!
!   if((2*(i-i/2)== i).and. i>1 .and. i< (nint-1))then
!     ee = ee + 2.d0*exp(-exponent)*den*2.d0*h/3.d0
!  else if ((2*(i-i/2)== i+1).and. i>1 .and. i< (nint-1)) then
!    ee = ee + 2.d0*exp(-exponent)*den*4.d0*h/3.d0
!  else if ( i==1 .or. i ==(nint-1)) then
!    ee = ee + 2.d0*exp(-exponent)*den*1.d0*h/3.d0
!  end if
! end do

if(allocated(x))deallocate(x)
if(allocated(w))deallocate(w)
end subroutine E2body
!
subroutine E1body(largeomega1, RA, RB, en)
implicit none
double precision, dimension(3,3), intent(in) :: largeomega1
double precision,dimension(3), intent(in) :: RA, RB
double precision, intent(out) :: en
integer :: i, nint
double precision, dimension(:), allocatable ::x, w
double precision, dimension(6,6) :: intermedom, intermedom1
double precision, dimension(6) :: RAB
double precision :: den, u, exponent
double precision, dimension(6,3) :: cmatrix

en = 0.d0
if(allocated(x))deallocate(x)
if(allocated(w))deallocate(w)
! if(allocated(x1))deallocate(x1)
! if(allocated(w1))deallocate(w1)
rab(1:3)= ra
rab(4:6)= rb
!nint = 14
nint =500
allocate(x(nint),w(nint))
 !allocate(x1(nint1),w1(nint1))
 !call gl_points(nint, x, w)
call simpson1by3(nint, x, w)
 do i = 1, nint
  u = (1.d0-x(i))/(1.d0+x(i))
  intermedom =matrix_embed(3, largeomega1, 6, 1, 1)
 forall(i=4:6)intermedom(i, i) = u**2.d0
 cmatrix = matrix_combine( 3, largeomega1, u*u*i_matrix(3))
 den = 1.d0/dsqrt(get_det(largeomega1+u*u*i_matrix(3)))/(1.d0+x(i))**2.d0
 intermedom1=intermedom - matmul(cmatrix,matmul(inverse(largeomega1+u*u*i_matrix(3)),transpose(cmatrix)))
 exponent =dot_product(rab, matmul(intermedom1,rab))
 !print*, 'den =', den, "exponential", exp(-exponent)
 en = en + 2.d0*w(i)*den*exp(-exponent)
end do
!write(*, *)'en 1 =', en
! en = 0.d0
! do i = 1, nint1-1
!   h=1.d0/(nint1-1)
!   xx= -1.d0+h*i
!   u = (1.d0-xx)/(1.d0+xx)
! !  write(*,*)i, xx, h, u
!   intermedom =matrix_embed(3, largeomega1, 6, 1, 1)
!   forall(i=4:6)intermedom(i, i) = u**2.d0
!   cmatrix = matrix_combine( 3, largeomega1, u*u*i_matrix(3))
!   den = 1.d0/dsqrt(get_det(largeomega1+u*u*i_matrix(3)))
!   intermedom1=intermedom - matmul(cmatrix,matmul(inverse(largeomega1+u*u*i_matrix(3)),transpose(cmatrix)))
!   exponent =dot_product(rab, matmul(intermedom1,rab))
!    print*, 'den =', den, "exponential", exp(-exponent)
!
!   if((2*(i-i/2)== i).and. i>1 .and. i< (nint1-1))then
!    en = en + 2.d0*den*exp(-exponent)*2.d0*h/(3.d0*(1.d0+xx)**2.d0)
!  else if ((2*(i-i/2)== i+1).and. i>1 .and. i< (nint1-1)) then
!    en = en + 2.d0*den*exp(-exponent)*4.d0*h/(3.d0*(1.d0+xx)**2.d0)
!  else if ( i==1 .or. i ==(nint1-1)) then
!    en = en + 2.d0*den*exp(-exponent)*1.d0*h/(3.d0*(1.d0+xx)**2.d0)
!  end if
! end do
! write(*, *)'en 2 =', en
! stop

if(allocated(x))deallocate(x)
if(allocated(w))deallocate(w)
! if(allocated(x1))deallocate(x1)
! if(allocated(w1))deallocate(w1)
end subroutine E1body
!
subroutine full_coulomb(n, coords, C, w, w0, a0, rvdw0, alpha, beta, version,dampswitch, ecoul, en, ee,nn)
implicit none
integer, intent(in) :: n
character(len=*), intent(in) :: version
double precision, dimension(n), intent(in) :: w0, a0, rvdw0
double precision, dimension(n, 3), intent(in) :: coords
double precision, dimension(3*n,3*n), intent(in) :: C
double precision, dimension(3*n), intent(in) :: w
double precision, intent(in) :: alpha, beta, dampswitch
double precision, intent(out) :: ecoul, ee, en, nn
double precision,dimension(3) :: RA, RB
double precision :: MA, MB, dt2b, dt1b, erep, eatt1, eatt2, nrep, sigma
type(scalar) :: fdamp
integer :: A, B, i, j
double precision, dimension(3*n,3*n) :: wmatrix, bigomega
double precision, dimension(6,6) :: largeomega
double precision, dimension(3,3) :: largeomega1, largeomega2

double precision, dimension(6) :: mab
integer, dimension(6) ::AB
integer, dimension(3*(n-2)) ::notAB
integer, dimension(3) ::AA, BB
integer, dimension(3*(n-1)) ::notAA, notBB

erep = 0.d0
eatt1 = 0.d0
eatt2 = 0.d0
nrep =0.d0
ecoul = 0.d0
ee = 0.d0
en = 0.d0
nn = 0.d0
!print*,"w0 = ",w0
call diag_double(w, wmatrix)
!call printmat(C,"C Matrix")
bigomega=matmul(C,matmul(wmatrix, transpose(C)))
!call printmat(bigomega,"bigomega")
!Stop
do A = 1, n
 do B = A+1, n
  RA = coords(A,:)
  RB = coords(B,:)
  mA = 1.d0/(a0(A)*w0(A)*w0(A))
  mB = 1.d0/(a0(B)*w0(B)*w0(B))
  mab = combine_vector(3, my_repeat(3, ma), 3, my_repeat(3, mb))
  nn = nn + 1.d0/dsqrt(dot_product(RA-RB, RA-RB))
  !write(*, *)"A : ", A, " B : ",B," nn = ", 1.d0/dsqrt(dot_product(RA-RB, RA-RB))

  !! <AB|(1/|rA-rB|)|AB>
  AB(:)= (/ (3*(A-1)+i, i = 1, 3),  (3*(B-1)+i, i = 1, 3) /)
  if(n == 2)then
   largeomega = bigomega(AB, AB)
   dt2b = 1.d0
  else
    notAB(:) = (/(i, i = 1, 3*(A-1)),(i, i = 3*A+1, 3*(B-1)),(i, i = 3*B+1, 3*n)/)
!    print*,"notAB = ",notAB
    largeomega = bigomega(AB, AB)-matmul(bigomega(AB,notAb),matmul(inverse(bigomega(notAB, notAB)),bigomega(notAB, Ab)))
!    write(*, *)"bigomega(notAB, notAB) before", bigomega(notAB, notAB)
!call printmat(largeomega,"largeomega")

!read(*, *)dummy
    dt2b= get_det(bigomega(notAB, notAB))
  end if
!
!    write(*, *)"bigomega(notAB, notAB) after", bigomega(notAB, notAB)
  do i = 1, 6
   do j = 1, 6
    largeomega(i,j)= dsqrt(mAb(i)*mAb(j))*largeomega(i,j)
   end do
  end do
!  call printmat(largeomega,"largeomega")
!  read(*, *)dummy

   call E2body(largeomega, RA, RB, ee)
  ! print*,'ee == ',ee
!   read(*, *)dummy
   ee=ee*(mA*mb)**1.5d0
!**1.5d0
!print*,'ee == ',ee
!read(*, *)dummy
   ee = ee/dsqrt(dt2b)
!   print*,'ee == ',ee, 'erep =', erep
  ! read(*, *)dummy
!   stop
   erep = erep+ee
!   print*,'ee == ',erep, RA, RB, A, B
!   read(*, *)dummy
   ee = 0.d0
  !! <A|rA-RB|^-1|A>
   AA(:) = (/(3*(A-1)+i, i =1, 3)/)
   notAA(:) = (/(i, i = 1, 3*(A-1)),(i, i = 3*A+1, 3*n)/)
!   write(*, *)"shape of bigomega(notAA, notAA) :",shape(bigomega(notAA, notAA))
!   do idx=1, size(bigomega(notAA, notAA), 1)
!    write(*, *) bigomega(notAA(idx), notAA)
!   end do

   dt1b = get_det(bigomega(notAA, notAA))

!   write(*, *)"Det", dt1b
!   stop
   largeomega1 = bigomega(AA, AA) - matmul(bigomega(AA, notAA),matmul(inverse(bigomega(notAA, notAA)),bigomega(notAA, AA)))
   largeomega1 = ma*largeomega1
   call E1body(largeomega1, RA, RB, en)
   eatt1 = eatt1 + en*(mA**1.5d0)/dsqrt(dt1b)
!   write(*, *)"added eatt1 = ", eatt1, "determinant = ", dt1b
!*(mA**1.5d0)/dsqrt(dt1b)
   en = 0.d0
  !! <B|rB-RA|^-1|B>
   BB(:) = (/(3*(B-1)+i, i =1, 3)/)
   notBB(:) = (/(i, i = 1, 3*(B-1)),(i, i = 3*B+1, 3*n)/)
   dt1b = get_det(bigomega(notBB, notBB))
   largeomega2 = bigomega(BB, BB) - matmul(bigomega(BB, notBB),matmul(inverse(bigomega(notBB, notBB)),bigomega(notBB, BB)))
   largeomega2 = mb*largeomega2
   call E1body(largeomega2, RB, RA, en)
   eatt2 = eatt2 + en*(mB**1.5d0)/dsqrt(dt1b)
!(mB**1.5d0)/dsqrt(dt1b)
   en = 0.d0
  sigma = beta*(rvdw0(A)+rvdw0(B))
 select case(version)
  case ('fermi')
   fdamp=(damping_fermi(RA-RB, sigma, alpha, .false.))
  case ('')
    fdamp%val = 1.d0
  case default
    fdamp%val = 1.d0
 end select
! write(*, *)"A : ", A, " B : ",B," eatt1 = ", eatt1, " eatt2 = ", eatt2," erep = ", erep, " nn = ", nn
!  & alpha*((sqrt(dot_product(RA-RB, RA-RB)))/sigma-1.d0)
  eatt1= (eatt1+eatt1*(fdamp%val-1.d0)*dampswitch)
  eatt2= (eatt2+eatt2*(fdamp%val-1.d0)*dampswitch)
  erep = (erep+erep*(fdamp%val-1.d0)*dampswitch)
  nn = nn + nn*(fdamp%val-1.d0)*dampswitch
!  write(*, *)"A =", A,"B =", B, "fdamp = ", fdamp
 end do
end do

en = eatt1+eatt2
en = -en*2.d0*sqrt(product(w))/dsqrt(pi)
ee= erep*2.d0*sqrt(product(w))/dsqrt(pi)
ecoul = nn + en + ee
!write(*, *)"ecoul = ",ecoul,"e_n = ",en,"ee = ",ee,"nn = ",nn

!print*,"sigma", sigma,"damping =", damping_fermi(sqrt(dot_product(RA-RB, RA-RB)),sigma , alpha)
!ecoul = ecoul!*(1.d0-damping_fermi(sqrt(dot_product(RA-RB, RA-RB)),sigma , alpha))
!
end subroutine full_coulomb
!
subroutine U2(u, BigU2)
implicit none
double precision, intent(in) :: u
!double precision, dimension(3, 3) :: matA, matB
double precision, dimension(6, 6), intent(out) :: BigU2
double precision, dimension(6, 6) :: TempU2
TempU2 = 0.d0
BigU2=0.d0
TempU2 = matrix_embed(3, i_matrix(3),6,1, 1)
BigU2= BigU2+ TempU2
TempU2 = 0.d0
TempU2 = matrix_embed(3, i_matrix(3),6,4, 4)
BigU2= BigU2+ TempU2
TempU2 = 0.d0
TempU2 = matrix_embed(3, -i_matrix(3),6,4, 1)
BigU2= BigU2+ TempU2
TempU2 = 0.d0
TempU2 = matrix_embed(3, -i_matrix(3),6,1, 4)
BigU2= BigU2+ TempU2
BigU2= (u**2.d0)*BigU2
end subroutine U2


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
subroutine simpson1by3(n, x, w)
implicit none
integer :: i, n
double precision, dimension(n) :: x, w
double precision:: h, delta

delta = 1.d-6
h=2.0*(1.d0-delta)/((n-1)*1.d0)
do i = 1, n
 x(i)=-(1.d0-delta)+h*(i-1)
 if(2*(i-(i/2))==i)then
  w(i) = 2.d0*h/3.d0
 else
  w(i) = 4.d0*h/3.d0
 end if
end do
w(1)= 1.d0*h/3.d0
w(n)= w(1)
end subroutine simpson1by3




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
double precision, dimension(:, :), intent(in) :: A
double precision, intent(out) :: dt
integer, dimension(size(A, 1)**2) :: ipiv
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
Bigmatrix(rpos:rpos+(n-1), cpos:cpos+(n-1)) =matA

end function matrix_embed
!
function Matrix_invert(A)result(Ainv)
! Provides matrix inversion of A
implicit none
double precision, dimension(:,:), intent(in) :: A
double precision, dimension(size(A, 1),size(A, 1)):: Ainv, Testmat
integer :: n, info
integer, dimension(size(A,1)) :: ipiv
double precision, dimension(size(A,1)) :: work

ipiv= 1
work = 0.d0
n = size(A, 1)
Ainv = 0.d0
!if (n>1)then
!call dgetrf(n, n, A, n, ipiv, info)
!if(info /=0)then
!print*, "Matrix inversion failed ... Aborting, info =", info
!stop
!end if
!call DGETRI(n, A, n, ipiv, work, n**2, info)
!if(info /=0)then
! print*, "Matrix inversion failed ... Aborting, info =", info
! stop
!end if
! Ainv= A
!else if (n == 1)then
! Ainv= A
!end if

Testmat = A

call DGETRF(N,N,A,N,IPIV,info)
if(info .eq. 0) then
!write(*,*)"succeded"
else
 write(*,*)"failed"
end if

call DGETRI(N,A,N,IPIV,WORK,N,info)
if(info .eq. 0) then
!write(*,*)"succeded"
else
 write(*,*)"Inversion failed ... aborting!"
 stop
end if
Ainv=A
!Testmat = matmul(Ainv, Testmat)
!print*, "Shape of A", shape(A)
!print*, "(Testmat)^{-1}"
!print*, "======================"
!do k = 1, size(A, 1)
! write(*,'(1x,6(E10.4,1x))')Ainv(k, :)
!end do
!call printmat(Ainv,"A inverse")
!print*, "======================"

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
fmt = "(1x,"//trim(adjustl(fmt))//"(E13.5E4 ,1x))"
!print*, fmt
!stop
do i = 1, n
 write(*, trim(fmt)) A(i, :)
end do
end subroutine printmat

function is_symmetric(matrix)result(symm)
implicit none
double precision, dimension(:, :), intent(in) :: matrix
logical :: symm
integer :: i, j
symm = .true.
do i = 1, size(matrix, 1)
 do j = i+1, size(matrix, 1)
  if(dabs(matrix(i, j)-matrix(j, i))>=1.d-8)then
   symm = .false.
  end if
 end do
end do
end function is_symmetric

end module mbd_repulsion
