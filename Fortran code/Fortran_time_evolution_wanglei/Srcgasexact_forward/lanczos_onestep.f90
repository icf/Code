subroutine lanczos_onestep(step,direction,matsize,dim,nup,ndo,dt,Hlanc,Eigenvalue, phi,phi_time)
use lanc_param
implicit none

real(kind=8):: dt
integer:: step, matsize, dim, nup,ndo, direction
complex(kind=8):: phi_time(dim)

complex(kind=8), parameter:: Xi=cmplx(0d0,1d0)

integer:: i, j,k

integer:: diff

complex(kind=8)::  phi(lanc_size,dim)

real(kind=8):: Hlanc(lanc_size,lanc_size)
real(kind=8):: Eigenvalue(lanc_size)

complex(kind=8):: ExpHlanc(lanc_size,lanc_size)
real(kind=8):: res



if (mod(step,stepforcycle)==0) then 
diff=stepforcycle
else
diff=step-(step/stepforcycle)*stepforcycle
endif

 DO i=1,matsize
   DO j=1,matsize
       ExpHlanc(i,j)=0.0_rKind
         DO k=1,matsize
ExpHlanc(i,j) = ExpHlanc(i,j) + Hlanc(i,k)*exp(-dble(direction)*Xi*dt*(diff)*EigenValue(k))*Hlanc(j,k)
         END DO
   END DO
 END DO


res=0d0
phi_time=0d0
do k=1, dim
 do j=1, matsize
    phi_time(k)=phi_time(k)+Phi(j,k)*ExpHlanc(j,1)
 enddo 
!print *, 'k', phi_time(k)
res=res+abs(phi_time(k))**2
enddo 

if (mod(step,stepforcycle)==0) then 
!----perform lanczos on a give target state phi_time
call sub_lancdynamics(nup, ndo, dim,phi_time,matsize, Hlanc,Phi)
call reigen(Hlanc(1:matsize,1:matsize),matsize,eigenvalue(1:matsize))
endif

end subroutine 

