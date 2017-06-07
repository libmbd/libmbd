!\file mbd_e1.f90
!\notes Takes C, frequency, coords, mass and charge to calculate the E1
!\todo Write the whole goddamn code!

module mbd_repulsion
use matrix_util
use quadratures
use ios
use mbd
contains
subroutine E1_twobody(BigOmAB, RA, RB, Coul_en_2b)
implicit none
double precision, dimension(6, 6), intent(in) :: BigOmAB
double precision, dimension(3), intent(in) :: RA, RB
double precision, dimension(6):: RAB
double precision, dimension(6, 6) :: BigU2, HugeOmega, BigMat
double precision :: u, dt, exponent, Coul_en_2b
double precision, dimension(:), allocatable :: x, w
integer :: i, ninteg
!print*, "in E1_twobody"

Coul_en_2b = 0.d0
ninteg = 6
u = 1.d0

allocate(x(ninteg), w(ninteg))
RAB(1:3) =RA
RAB(4:6) = RB
call gl_points(ninteg, x, w)
do i = 1, ninteg
!get the roots of legendre
 u = (1.d0-x(i))/(1.d0+x(i))
 !print*, "u = ", u
 call U2(u, BigU2)
 BigMat=BigOmAB + BigU2
 call det(BigMat, dt) ! Add the U2 and get det
 !print*,"det in 2 body= ", dt
 Bigmat = Matrix_invert(BigMat)
 HugeOmega = BigOmAB - matmul(BigOmAB, matmul(BigMat, BigOmAB))
 exponent= dot_product(RAB,matmul(HugeOmega, RAB))
 coul_en_2b = coul_en_2b + w(i)*2.d0*(exp(-exponent)/dsqrt(dt))/(1.d0+u)**2
end do
deallocate(x, w)
print*, "Coulomb 2-body = ",coul_en_2b
end subroutine E1_twobody
!
subroutine E1_onebody(BigOmA, RA, RB, Coul_en_1b)
implicit none
double precision, dimension(3, 3), intent(in) :: BigOmA
double precision, dimension(3), intent(in) :: RA, RB
double precision, dimension(6):: RAB
double precision, dimension(3, 3) :: BigMat
double precision, dimension(6, 6) :: HugeOmega, zero_matrix
double precision, dimension(6, 3) :: om_u_matrix
double precision :: u, dt, exponent, Coul_en_1b
double precision, dimension(:), allocatable :: x, w
integer :: i, ninteg
!print*, "in E1_onebody"

Coul_en_1b = 0.d0
ninteg = 6
u = 1.d0
zero_matrix = 0.d0
HugeOmega = 0.d0
allocate(x(ninteg), w(ninteg))
RAB(1:3) =RA
RAB(4:6) = RB
call gl_points(ninteg, x, w)
do i = 1, ninteg
!get the roots of legendre
u = (1.d0-x(i))/(1.d0+x(i))

BigMat=BigOmA + (u**2)*i_matrix(3)
om_u_matrix = matrix_combine( 3, BigOmA, (u**2)*i_matrix(3))
call det(BigMat, dt)
!print*, "det in 1 body= ", dt
Bigmat = Matrix_invert(BigMat)
HugeOmega= matrix_embed(3, BigOmA, 6, 1, 1)+matrix_embed(3,(u**2)*i_matrix(3), 6, 4, 4)
HugeOmega = HugeOmega - matmul(om_u_matrix, matmul(BigMat, transpose(om_u_matrix)))
exponent= dot_product(RAB,matmul(HugeOmega, RAB))
coul_en_1b = coul_en_1b + w(i)*2.d0*(exp(-exponent)/dsqrt(dt))/(1.d0+u)**2
end do
deallocate(x, w)
print*, "Coulomb 1-body = ",coul_en_1b
end subroutine E1_onebody
!
subroutine fullcoulomb(natom, C, coords, charge, mass, omega, omzero, e1)
implicit none
integer, intent(in) :: natom
double precision, dimension(natom), intent(in) :: charge, mass, omzero
double precision, dimension(3*natom, 3*natom), intent(in) :: C
double precision, dimension(3*natom), intent(in) :: omega
double precision, dimension(natom,3) :: coords
double precision, intent(out) :: e1
integer :: i, j, A, B
double precision :: ninv, coul_en_2b, coul_en_1b, coul_en
integer, dimension(6) :: AB
integer, dimension(:),allocatable :: notAB
double precision, dimension(3*natom, 3*natom) :: BigOmega
double precision, dimension(6, 6) :: m_matrix, BigOmAB
double precision, dimension(3) :: RA, RB
double precision, dimension(3, 3) :: BigOmA
integer, dimension(3) :: AA
integer, dimension(3*(natom-1)) :: notAA

coul_en_2b = 0.d0
coul_en_1b = 0.d0
coul_en = 0.d0

call diag_double(omega, BigOmega)
Bigomega = matmul(transpose(C),matmul(BigOmega,C))
if (natom > 2)then
 allocate(notAB(3*(natom-2)))
end if
do A = 1, natom
AA(:) = (/ (3*(A-1)+i, i = 1, 3)/)
notAA(:) = (/ (i, i = 1, 3*(A-1)),(i, i = 3*A+1, 3*natom)/)
BigOmA = BigOmega(AA, AA)-matmul(BigOmega(AA, notAA),matmul(matrix_invert(BigOmega(notAA, notAA)),BigOmega(notAA, AA)))
BigOmA = sqrt(mass(A))*BigOmA
 do B = A+1, natom
  ra = coords(A, :)
  rb = coords(B, :)
  AB(:)= (/ (3*(A-1)+i, i = 1, 3),  (3*(B-1)+i, i = 1, 3) /)
  if(natom > 2)then
notAB(:) = (/ (i, i = 1, 3*(A-1)),  (i, i = 3*A+1, 3*(B-1)), (i, i = 3*B+1, 3*natom) /)
BigOmAB = BigOmega(AB, AB)-matmul(BigOmega(AB, notAB),matmul(matrix_invert(BigOmega(notAB, notAB)),BigOmega(notAB, AB)))
call det(BigOmega(notAB, notAB), ninv)
  else
BigOmAB = BigOmega(AB, AB)
  end if
call diag_double(combine_vector(3,repeat(3, sqrt(mass(A))), 3, repeat(3, sqrt(mass(B)))), m_matrix)
BigOmAB=matmul(m_matrix, BigOmAB)
!print*, "calling E1_twobody"
call E1_twobody(BigOmAB, RA, RB, coul_en_2b)
coul_en = coul_en + coul_en_2b/1.d0
coul_en_2b = 0.d0
!print*, "calling E1_onebody"
call E1_onebody(BigOmA, RA, RB, coul_en_1b)
coul_en = coul_en -2.d0*coul_en_1b
coul_en_1b = 0.d0
coul_en= coul_en+ 1.d0/sqrt(dot_product(RA-RB, RA-RB))
 end do
! onebody
end do
print*, "coul_en = ", coul_en
e1 = coul_en
if(allocated(notAB)) deallocate (notAB)
end subroutine fullcoulomb
!
subroutine U2(u, BigU2)
implicit none
double precision, intent(in) :: u
!double precision, dimension(3, 3) :: matA, matB
double precision, dimension(6, 6), intent(out) :: BigU2
double precision, dimension(6, 6) :: TempU2
integer :: i
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
BigU2= (u**2)*BigU2
end subroutine U2
!
! Jan Herman's Dipole energy code 
! copied as it is
real(8) function get_dipole_energy_coupled_osc(R, a0, w, w_t, C) result(ene)
real(8), intent(in) :: R(:, :), a0(size(R, 1)), w(size(R, 1)), w_t(3*size(R, 1))
real(8), intent(in) :: C(3*size(R, 1), 3*size(R, 1))

real(8) :: T(size(C, 1), size(C, 1))
integer :: A, B, i, j, N

T(:, :) = 0.d0
N = size(R, 1)
call add_dipole_matrix('', 'dip,gg', R, a0, w, relay=T)
do  A = 1, N
do B = 1, N
i = 3*(A-1)
j = 3*(B-1)
T(i+1:i+3, j+1:j+3) = w(A)*w(B)*sqrt(a0(A)*a0(B))*T(i+1:i+3, j+1:j+3)
end do
end do
! call print_matrix('T', T)
T = matmul(matmul(transpose(C), T), C)
! call print_matrix('T_t', T)
ene = sum(diag(T)/(4*w_t))
end function

end module mbd_repulsion
