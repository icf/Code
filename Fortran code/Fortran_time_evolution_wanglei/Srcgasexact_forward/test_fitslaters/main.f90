program main
use trans4calfun
implicit none

  
   integer::N                     ! num of variables 
   integer::M                     ! number of functions
   
   integer:: i, j

  complex(kind=8):: overlap1, overlap2, overlap12
  real(kind=8), external:: rannyu
  complex(kind=8),allocatable:: phi(:,:)

Norb=11
Npar=2

N=2*Norb*Npar ! real/imag part* Matrix size
M=1

allocate(phi_1(Norb, Npar))
allocate(phi_2(Norb, Npar))
allocate(phi(Norb, Npar))

do i=1,Norb
  do j=1, Npar
phi_1(i,j)=cmplx(rannyu(),rannyu())
phi_2(i,j)=cmplx(rannyu(),rannyu())
phi(i,j)=0.5d0*(phi_1(i,j)+phi_2(i,j))
  enddo 
enddo


!do i=1, Norb
!  do j=1, Npar
!    print * ,i, j, phi_1(i,j)
!  enddo 
!enddo 


!---get the offset of the target function (the one done not dependent on variables)
call  cal_overlap(Norb,Npar,phi_1,phi_1,overlap1)
call  cal_overlap(Norb,Npar,phi_2,phi_2,overlap2)
call  cal_overlap(Norb,Npar,phi_1,phi_2,overlap12)

offset=dble(overlap1)/4d0+dble(overlap2)/4d0+dble(overlap12)/2d0


call  slater_det_fitting(N,M,phi,Norb,Npar)

!print *, 'phi_target'
!do i=1, Norb
!  do j=1, Npar
!    print * ,i, j, phi_target(i,j)
!  enddo 
!enddo 


deallocate(phi_1)
deallocate(phi_2)
deallocate(phi)


end program 
