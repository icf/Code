module trans4calfun
implicit none

!integer:: Nsite,nup,ndo
!integer:: numupdo(2)
complex(kind=8), allocatable:: phi_ref(:)! refenence state 
real(kind=8):: offset_ref1, offset_ref2

end module

!======================================================================
subroutine  slater_det_fitting(N,M,dim,phi_target,phiref,offset1,offset2,f)
use system_parameters ! use nup, ndo , numupdo there
use hubbard_param ! use Nsite there
use trans4calfun  
implicit none
  
   integer,intent(IN)::N                     ! num of variables 
   integer,intent(IN)::M                     ! number of functions =1 here
   integer,intent(IN):: dim
   
   complex(kind=8),intent(IN):: phiref(dim)
   complex(kind=8),intent(INOUT):: phi_target(Nsite,max(nup,ndo),2)
   real(kind=8), intent(IN):: offset1, offset2

   real(kind=8)::x(N),f(M)

   real(kind=8)::h=1d-8
   real(kind=8)::dmax=10d0
   real(kind=8)::acc=1d-10
  
  
   integer:: maxfun=10**6
   integer::iprint


   real(kind=8)::W(2*m*n+2*n*n+2*m+5*n+10)

  integer:: spin,i, j, counter
!-------------------------------------------------------------
allocate(phi_ref(dim))

phi_ref=phiref
offset_ref1=offset1
offset_ref2=offset2

 !print *, 'm, n', m, n
 !pause

iprint=0

counter=1
do spin=1,2
 do j=1, numupdo(spin)
  do i=1,Nsite
 x(2*counter-1)=dble(phi_target(i,j,spin))
 x(2*counter)=dimag(phi_target(i,j,spin))
 counter=counter+1   
  enddo 
 enddo
enddo 



call va05ad(m,n,f,x,h,dmax, acc, maxfun, iprint,w)

counter=1
do spin=1,2
 do j=1, numupdo(spin)
  do i=1,Nsite
 phi_target(i,j,spin)=cmplx(x(2*counter-1),x(2*counter))
 counter=counter+1   
  enddo 
 enddo
enddo 

!print *, 'fitting error', sum(abs(f(:))), iprint

deallocate(phi_ref)


end subroutine slater_det_fitting

!***************************************************************************

subroutine calfun(m, n, f, x) ! for va05
use system_parameters ! use nup, ndo , numupdo there
use hubbard_param ! use Nsite there
use trans4calfun
use calphy_module
implicit none
   integer,intent(IN)::m,n
   real(kind=8)::x(n),f(m)
   integer:: spin,i, j, counter 


   complex(kind=8):: phi_target(Nsite,max(nup,ndo),2)
   complex(kind=8):: phi2hf(4)
complex(kind=8):: overlap_mix, dbocc_mix, kin_mix
complex(kind=8):: occ_mix(Nsite,2)




counter=1
do spin=1,2
 do j=1, numupdo(spin)
  do i=1,Nsite
 phi_target(i,j,spin)=cmplx(x(2*counter-1),x(2*counter))
 counter=counter+1   
  enddo 
 enddo
enddo 


call slater2fock(nup,ndo,4,phi_target(1:Nsite,1:nup,1),phi_target(1:Nsite,1:ndo,2),phi2hf)

call cal_den_mix(4, nup, ndo,phi_ref,phi2hf,occ_mix,dbocc_mix, overlap_mix,kin_mix)


!f(1)=offset_ref1-dble(overlap_mix)
!f(2)=offset_ref2-dimag(overlap_mix)


f(1)=1d0-dble(overlap_mix)
f(2)=0d0-dimag(overlap_mix)




!print *, 'x',x
!print *, 'error',sum(abs(f(:)))
!pause

end subroutine calfun
!=======================================================================

