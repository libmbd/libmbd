program main
!use ofunction
use quadratures
use matrix_util
use mbd_repulsion
implicit none
integer, parameter :: n = 2
integer :: i, j
double precision :: r, dt, a, e1, dip_en, erep, eatt
double precision, dimension(n) :: vector
double precision, dimension(n,n) :: matA, matB
double precision, dimension(3*n,3*n) :: C
double precision, dimension (3) :: RA, RB, RC
double precision, dimension(n, 3) :: coords

open(1, file='e1.txt')
r=0.d0
a = 1.d0/sqrt(2.d0)
!C(1, :)= [ 1.d0, 0.d0, 0.d0,0.d0,0.d0,0.d0]!,0.d0,0.d0,0.d0]
!
!C(2, :)= [ 0.d0, 1.d0, 0.d0,0.d0,0.d0,0.d0]!,0.d0,0.d0,0.d0]
!
!C(3, :)= [ 0.d0, 0.d0, 1.d0,0.d0,0.d0,0.d0]!,0.d0,0.d0,0.d0]
!
!C(4, :)= [ 0.d0, 0.d0, 0.d0,1.d0,0.d0,0.d0]!,0.d0,0.d0,0.d0]
!
!C(5, :)= [ 0.d0, 0.d0, 0.d0,0.d0,1.d0,0.d0]!,0.d0,0.d0,0.d0]
!
!C(6, :)= [ 0.d0, 0.d0, 0.d0,0.d0,0.d0,1.d0]!,0.d0,0.d0,0.d0]
!!!C(7, :)= [ 0.d0, 0.d0, 0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0]
!!
!!C(8, :)= [ 0.d0, 0.d0, 0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0]
!!
!!C(9, :)= [ 0.d0, 0.d0, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0]
C = 0.d0
C(1, 3) = -a
C(1, 5) = -a
C(2, 2) = a
C(2, 4) = a
C(3, 1) = -a
C(3, 6) = -a
C(4, 3) = a
C(4, 5) = -a
C(5, 2) = -a
C(5, 4) = a
C(6, 1) = -a
C(6, 6) = a
call printmat(matmul(transpose(C), C), "C:")
RA = [0.0d0, 0.0d0, 0.0d0]
RC =[0.0d0, 1.0d0, 0.0d0]
do i = 1, 100
r = 0.d0 +i*0.1d0
RB = [0.d0, 0.0d0, r]

coords = reshape([RA, RB], shape(coords))
dip_en = get_dipole_energy_coupled_osc(coords, repeat(n, 1.d0), repeat(n, 2.d0), repeat(3*n, 2.d0), C)
call fullcoulomb(n, C, coords, repeat(n,1.d0), repeat(n, 1.d0), repeat(3*n, 2.d0), repeat(n, 2.d0), e1, eatt, erep)
print*, "shape of coords = ", shape(coords)

write(1, *)r, e1, eatt, erep, dip_en, 1.d0/r
end do
close(1)
end program
