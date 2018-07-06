program main
implicit none

complex(kind=8),parameter:: Xi=cmplx(0d0,1d0)
real(kind=8),parameter:: dt=0.1d0
real(kind=8),parameter:: tmax=10d0
real(kind=8),parameter:: Thop=-1d0
real(kind=8),parameter:: U=2d0

real(kind=8):: Ham1(4,4), Eigenvalue1(4)
real(kind=8):: Ham2(4,4), Eigenvalue2(4)

complex(kind=8)::  phi0(4), phi_plus(4), phi_minus(4)
complex(kind=8):: ExpHam1(4,4),ExpHam2(4,4) 

real(kind=8):: time,den_plus, den_minus

integer:: i,j,k

real(kind=8):: res


!---------------------------------
Ham1(1,:)=(/U, Thop, -Thop, 0d0/)
Ham1(2,:)=(/Thop,0d0,0d0,Thop/)
Ham1(3,:)=(/-Thop,0d0,0d0,-Thop/)
Ham1(4,:)=(/0d0, Thop, -Thop, U/)

!---------------------------------
Ham2(1,:)=(/-U, Thop, -Thop, 0d0/)
Ham2(2,:)=(/Thop,0d0,0d0,Thop/)
Ham2(3,:)=(/-Thop,0d0,0d0,-Thop/)
Ham2(4,:)=(/0d0, Thop, -Thop, -U/)


!do i=1,4
!print *, Ham1(i,:)-Ham2(i,:)
!enddo 

call reigen(Ham1,4,eigenvalue1)
call reigen(Ham2,4,eigenvalue2)


!phi0=(/0d0,0d0,1d0,0d0/)

phi0=(/1d0/sqrt(2d0),0d0,0d0,1d0/sqrt(2d0)/)

!do k=1,4
!print *, 'phi0',phi0(k)
!enddo 

open(11,file='time_den1.dat')

do time=dt,tmax,dt

ExpHam1=0d0
ExpHam2=0d0

DO i=1,4
   DO j=1,4
         DO k=1,4
ExpHam1(i,j) = ExpHam1(i,j) + Ham1(i,k)*exp(-Xi*time*EigenValue1(k))*Ham1(j,k)
ExpHam2(i,j) = ExpHam2(i,j) + Ham2(i,k)*exp(-Xi*time*EigenValue2(k))*Ham2(j,k)
         END DO
   END DO
 END DO


phi_plus=0d0
phi_minus=0d0
do j=1,4
  do k=1,4

phi_plus(j)=phi_plus(j)+0.5d0*(ExpHam1(j,k)+ExpHam2(j,k)) *phi0(k)
phi_minus(j)=phi_minus(j)+0.5d0*(ExpHam1(j,k)-ExpHam2(j,k)) *phi0(k)

!phi(j)=phi(j)+ExpHam2(j,k) *phi0(k)
 enddo 
enddo 


!res=0d0
!do j=1,4
!res=res+abs(phi(j))**2
!enddo
!phi=phi/sqrt(res)


!den_plus=abs(phi_plus(1))**2+abs(phi_plus(2))**2
!den_minus=abs(phi_minus(1))**2+abs(phi_minus(2))**2

den_plus=abs(phi_plus(1))**2+abs(phi_plus(4))**2
den_minus=abs(phi_minus(1))**2+abs(phi_minus(4))**2



!print *, sum(abs(phi(:))**2)
write(11,*) time, den_plus, den_minus


!phi0=phi !only when we do the iteration

enddo 

close(11)

end program
