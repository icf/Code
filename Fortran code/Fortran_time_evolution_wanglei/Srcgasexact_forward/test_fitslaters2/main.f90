program main
use trans4calfun
implicit none

  
   integer::N                     ! num of variables 
   integer::M                     ! number of functions
   
   integer:: i, j, spin

  complex(kind=8):: overlap1(2), overlap2(2), overlap12(2)
  real(kind=8), external:: rannyu
  complex(kind=8),allocatable:: phi_target(:,:,:)
  complex(kind=8),allocatable:: phi_1(:,:,:), phi_2(:,:,:)
  real(kind=8):: offset

Nsite=10
nup=2
ndo=1
numupdo=(/nup,ndo/)

N=2*Nsite*(nup+ndo) ! real/imag part* Matrix size
M=1

allocate(phi_1(Nsite,max(nup,ndo),2))
allocate(phi_2(Nsite,max(nup,ndo),2))
allocate(phi_target(Nsite,max(nup,ndo),2))

do spin=1,2
  do j=1, numupdo(spin)
    do i=1,Nsite
phi_1(i,j,spin)=cmplx(rannyu(),rannyu())
phi_2(i,j,spin)=cmplx(rannyu(),rannyu())
phi_target(i,j,spin)=0.5d0*(phi_1(i,j,spin)+phi_2(i,j,spin))
  enddo 
 enddo 
enddo


!do i=1, Norb
!  do j=1, Npar
!    print * ,i, j, phi_1(i,j)
!  enddo 
!enddo 


!---get the offset of the target function (the one done not dependent on variables)
do spin=1,2
call  cal_overlap(Nsite,numupdo(spin),phi_1(1:Nsite,1:numupdo(spin),spin),phi_1(1:Nsite,1:numupdo(spin),spin),overlap1(spin))
call  cal_overlap(Nsite,numupdo(spin),phi_2(1:Nsite,1:numupdo(spin),spin),phi_2(1:Nsite,1:numupdo(spin),spin),overlap2(spin))
call  cal_overlap(Nsite,numupdo(spin),phi_1(1:Nsite,1:numupdo(spin),spin),phi_2(1:Nsite,1:numupdo(spin),spin),overlap12(spin))
enddo 


offset=dble(product(overlap1(:)))/4d0+dble(product(overlap2(:)))/4d0+dble(product(overlap12(:)))/2d0

call slater_det_fitting(N,M,phi_target,phi_1,phi_2,offset)

!print *, 'phi_target'
!do i=1, Norb
!  do j=1, Npar
!    print * ,i, j, phi_target(i,j)
!  enddo 
!enddo 


deallocate(phi_1)
deallocate(phi_2)
deallocate(phi_target)

end program 
