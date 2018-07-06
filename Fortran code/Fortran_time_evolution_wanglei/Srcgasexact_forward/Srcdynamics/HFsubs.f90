subroutine initialize_expK
use Hubbard_param
use HF_module
use flags_param
implicit none

real(kind=8):: KKmat(Nsite,Nsite), Eigenvalue(Nsite)

integer:: i,j,k

Kmat=0d0
do i=1,Nsite-1
 Kmat(i,i+1)=-Thop
enddo 

!PBC
Kmat(1,Nsite)=-Thop

 do i=1,Nsite
   do j=i+1,Nsite
     Kmat(j,i)=Kmat(i,j)
   enddo 
 enddo 

! do i=1,Nsite
! print *, 'Kmat',i,Kmat(i,:)
! enddo 
 
KKmat=Kmat
call reigen(KKmat,Nsite,Eigenvalue)

 ExpK=0d0
 DO i=1,Nsite
   DO j=1,Nsite
         DO k=1,Nsite
ExpK(i,j) = ExpK(i,j) + KKmat(i,k)*exp(-Xi*dt/dble(TrotterOrder)*EigenValue(k))*KKmat(j,k)
         END DO
   END DO
 END DO

end subroutine 

subroutine initialize_phi
use HF_module
implicit none

integer:: i,site,num, spin
phiT=0d0
     do spin=1,2
      do num=1,numupdo(spin)
       phiT(num,num,spin) = 1d0
      enddo 
     enddo 

end subroutine 





subroutine propagationHF_onestep(expKHF) ! propagation HF state one step
use Hubbard_param
use HF_module
use flags_param
implicit none

complex(kind=8), intent(IN):: expKHF(Nsite,Nsite,2)
integer:: spin

   do spin=1,2
   phiT(1:Nsite,1:numupdo(spin),spin)=matmul(expKHF(1:Nsite,1:Nsite,spin),phiT(1:Nsite,1:numupdo(spin),spin)) 
   enddo 

end subroutine propagationHF_onestep



subroutine cal_expKHF(occ,expKHF)
use Hubbard_param
use HF_module
use flags_param
implicit none
real(kind=8),intent(IN):: occ(Nsite,2)
complex(kind=8),intent(OUT):: expKHF(Nsite,Nsite,2)

real(kind=8):: KHFmat(Nsite,Nsite), Eigenvalue(Nsite)
integer:: i,j,k, spin

ExpKHF=0d0

  do spin=1,2
     KHFmat=Kmat
    do i=1,Nsite
      KHFmat(i,i)=Kmat(i,i)+U*occ(i,3-spin)
    enddo 

call reigen(KHFmat,Nsite,Eigenvalue)

 DO i=1,Nsite
   DO j=1,Nsite
         DO k=1,Nsite
ExpKHF(i,j,spin) = ExpKHF(i,j,spin) + KHFmat(i,k)*exp(-Xi*dt*EigenValue(k))*KHFmat(j,k)
         END DO
   END DO
 END DO
 enddo 
end subroutine 


