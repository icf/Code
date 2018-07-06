module trans4calfun
implicit none

integer:: Nsite,nup,ndo
integer:: numupdo(2)
complex(kind=8), allocatable:: phi_ref1(:,:,:), phi_ref2(:,:,:) ! refenence state 
real(kind=8):: offset_ref

end module

!======================================================================
subroutine  slater_det_fitting(N,M,phi,phi1,phi2,offset)
use trans4calfun  
implicit none
  
   integer,intent(IN)::N                     ! num of variables 
   integer,intent(IN)::M                     ! number of functions =1 here
   
   complex(kind=8),intent(IN):: phi1(Nsite,max(nup,ndo),2), phi2(Nsite,max(nup,ndo),2)
   complex(kind=8),intent(INOUT):: phi(Nsite,max(nup,ndo),2)
   real(kind=8), intent(IN):: offset

   real(kind=8)::x(N),f(M)

   real(kind=8)::h=1d-6
   real(kind=8)::dmax=10d0
   real(kind=8)::acc=1d-16
  
  
   integer:: maxfun=10**6
   integer::iprint


   real(kind=8)::W(2*m*n+2*n*n+2*m+5*n+10)

  integer:: spin,i, j, counter
!-------------------------------------------------------------
allocate(phi_ref1(Nsite,max(nup,ndo),2))
allocate(phi_ref2(Nsite,max(nup,ndo),2))

phi_ref1=phi1
phi_ref2=phi2
offset_ref=offset
! print *, 'm, n', m, n

iprint=0

counter=1
do spin=1,2
 do j=1, numupdo(spin)
  do i=1,Nsite
 x(2*counter-1)=dble(phi(i,j,spin))
 x(2*counter)=dimag(phi(i,j,spin))
 counter=counter+1   
  enddo 
 enddo
enddo 



call va05ad(m,n,f,x,h,dmax, acc, maxfun, iprint,w)

counter=1
do spin=1,2
 do j=1, numupdo(spin)
  do i=1,Nsite
 phi(i,j,spin)=cmplx(x(2*counter-1),x(2*counter))
 counter=counter+1   
  enddo 
 enddo
enddo 



print *, 'fitting error', sum(abs(f(:))), iprint

deallocate(phi_ref1)
deallocate(phi_ref2)


end subroutine slater_det_fitting

!***************************************************************************

subroutine calfun(m, n, f, x) ! for va05
use trans4calfun
implicit none
   integer,intent(IN)::m,n
   real(kind=8)::x(n),f(m)
   integer:: spin,i, j, counter 


   complex(kind=8):: overlap1(2), overlap2(2), overlap_target(2)
   complex(kind=8):: phi_target(Nsite,max(nup,ndo),2)

counter=1
do spin=1,2
 do j=1, numupdo(spin)
  do i=1,Nsite
 phi_target(i,j,spin)=cmplx(x(2*counter-1),x(2*counter))
 counter=counter+1   
  enddo 
 enddo
enddo 

do spin=1,2

call  cal_overlap(Nsite,numupdo(spin),phi_ref1(1:Nsite,1:numupdo(spin),spin),phi_target(1:Nsite,1:numupdo(spin),spin),overlap1(spin))
!print *,'1t'

call  cal_overlap(Nsite,numupdo(spin),phi_ref2(1:Nsite,1:numupdo(spin),spin),phi_target(1:Nsite,1:numupdo(spin),spin),overlap2(spin))


call  cal_overlap(Nsite,numupdo(spin),phi_target(1:Nsite,1:numupdo(spin),spin),phi_target(1:Nsite,1:numupdo(spin),spin),overlap_target(spin))

enddo 


f(1)=dble(product(overlap_target(:))) - dble(product(overlap1(:))) - dble(product(overlap2(:)))+ offset_ref


!print *, 'x',x
!print *, 'error',sum(abs(f(:)))
!pause

end subroutine calfun
!=======================================================================

