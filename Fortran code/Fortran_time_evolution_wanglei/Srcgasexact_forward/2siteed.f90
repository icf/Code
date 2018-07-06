program main
implicit none

complex(kind=8),parameter:: Xi=cmplx(0d0,1d0)
real(kind=8),parameter:: dt=0.01d0
real(kind=8),parameter:: Thop=-1d0
real(kind=8),parameter:: U=2d0

real(kind=8):: Ham(4,4), Eigenvalue(4)
complex(kind=8)::  phi0(4), phi(4),ExpHam(4,4) 

real(kind=8):: time,den

integer:: i,j,k

Ham(1,:)=(/U, Thop, -Thop, 0d0/)
Ham(2,:)=(/Thop,0d0,0d0,Thop/)
Ham(3,:)=(/-Thop,0d0,0d0,-Thop/)
Ham(4,:)=(/0d0, Thop, -Thop, U/)

do i=1,4
print *, Ham(i,:)
enddo 

call reigen(Ham,4,eigenvalue)

ExpHam=0d0
DO i=1,4
   DO j=1,4
         DO k=1,4
ExpHam(i,j) = ExpHam(i,j) + Ham(i,k)*exp(-Xi*dt*EigenValue(k))*Ham(j,k)
         END DO
   END DO
 END DO


!phi0=(/1d0,0d0,0d0,0d0/)
phi0=(/1d0/sqrt(2d0),0d0,0d0,1d0/sqrt(2d0)/)

do k=1,4
print *, 'phi0',phi0(k)
enddo 

open(11,file='time_den1.dat')

do time=dt,10d0,dt
phi=0d0
do j=1,4
  do k=1,4
phi(j)=phi(j)+ExpHam(j,k)*phi0(k)
 enddo 
enddo 

den=abs(phi(1))**2+abs(phi(4))**2
!print *, sum(abs(phi(:))**2)
write(11,*) time, den

phi0=phi
enddo 

close(11)

end program
