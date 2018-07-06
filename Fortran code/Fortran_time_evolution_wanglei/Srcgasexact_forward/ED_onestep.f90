subroutine ed_onestep(direction,dim,phiio)
use system_parameters
implicit none

integer:: dim,direction
complex(kind=8):: phiio(dim)

complex(kind=8):: phi0(dim)


integer::  j,k

real(kind=8):: res

phi0=phiio
phiio=0d0

if (direction==1) then 

do j=1,dim
  do k=1,dim
phiio(j)=phiio(j)+ExpHam(j,k)*phi0(k)
 enddo 
enddo 

else 

do j=1,dim
  do k=1,dim
phiio(j)=phiio(j)+conjg(ExpHam(j,k))*phi0(k)
 enddo 
enddo 

endif

res=0d0
do j=1, dim 
res=res+ abs(phiio(j))**2
enddo 
!print *, 'norm', res

if (abs(res-1d0)>1d-8) then 
        print *, 'the norm of the exact phi is wrong'
        stop
endif


end subroutine 

