module trans4calfun
implicit none

integer:: Norb, Npar
complex(kind=8), allocatable:: phi_1(:,:), phi_2(:,:)
real(kind=8):: offset

end module

!======================================================================
subroutine  slater_det_fitting(N,M,phi,Norb,Npar)
!use trans4calfun  
implicit none
  
   integer,intent(IN)::N                     ! num of variables 
   integer,intent(IN)::M                     ! number of functions =1 here
   
   integer,intent(IN):: Norb, Npar
   complex(kind=8),intent(INOUT):: phi(Norb,Npar)

   real(kind=8)::x(N),f(M)

   real(kind=8)::h=1d-4
   real(kind=8)::dmax=10d0
   real(kind=8)::acc=1d-16
  
  
   integer:: maxfun=10**6
   integer::iprint


   real(kind=8)::W(2*m*n+2*n*n+2*m+5*n+10)

  integer:: i, j, counter
!-------------------------------------------------------------
! print *, 'm, n', m, n

iprint=0

counter=1
do i=1,Norb
  do j=1, Npar
 x(2*counter-1)=dble(phi(i,j))
 x(2*counter)=dimag(phi(i,j))
!  x(counter)=dble(phi_target(i,j))
   counter=counter+1   
  enddo 
enddo


call va05ad(m,n,f,x,h,dmax, acc, maxfun, iprint,w)


counter=1
do i=1,Norb
  do j=1, Npar
 phi(i,j)= cmplx(x(2*counter-1), x(2*counter))
 !phi_target(i,j)= x(counter) 
 counter=counter+1   
  enddo 
enddo

print *, 'fitting error', sum(abs(f(:))), iprint

end subroutine slater_det_fitting

!***************************************************************************

subroutine calfun(m, n, f, x) ! for va05
use trans4calfun
implicit none
   integer,intent(IN)::m,n
   real(kind=8)::x(n),f(m)
   integer:: i, j, counter 


   complex(kind=8):: overlap1, overlap2, overlap_target
   complex(kind=8):: phi_target(Norb,Npar)

counter=1
do i=1,Norb
  do j=1, Npar
   phi_target(i,j)= cmplx(x(2*counter-1), x(2*counter))
   !phi_target(i,j)= x(counter)
   counter=counter+1   
!print *, i, j, phi_target(i,j)
  enddo 
enddo


call  cal_overlap(Norb,Npar,phi_1,phi_target,overlap1)
!print *,'1t'

call  cal_overlap(Norb,Npar,phi_2,phi_target,overlap2)
!print *,'2t'

call  cal_overlap(Norb,Npar,phi_target,phi_target,overlap_target)
!print *,'tt'


f(1)=dble(overlap_target) - dble(overlap1) - dble(overlap2)+ offset


!print *, 'x',x
!print *, 'error',sum(abs(f(:)))
!pause

end subroutine calfun
!=======================================================================

