subroutine initialize
use system_parameters
use hubbard_param
implicit none


NAMELIST /SystemSettings/ Nwalkers,TrotterOrder,HS,dt,Maxsteps,Dir,IniState, Nsamples,Shift,Phaseless,Local_energy,PopContrlSteps, Maxshift, UpperBound, LowerBound
NAMELIST /FHParams/  Nsite, nup, ndo,Thop,U, Mu

OPEN(138,file='input')
READ(138,SystemSettings)
READ(138,FHParams)
CLOSE(138)

write(*, NML=SystemSettings)
write(*, NML=FHParams)

call ini_EDdynamics()

numupdo(1)=nup
numupdo(2)=ndo

allocate(phi_ini(Nsite,max(nup,ndo),2))
!allocate(phiT(Nsite,max(nup,ndo),2))
allocate(phi(Nsite,max(nup,ndo),2,Nwalkers))

!allocate(occT(Nsite,2))

allocate(weight(Nwalkers))
!allocate(wscaling(Nwalkers))
allocate(wphase(Nwalkers))
allocate(impfunc0(Nwalkers))
allocate(impfunc1(Nwalkers))


allocate(Kmat(Nsite,Nsite))
allocate(expK(Nsite,Nsite))

call initialize_phi_ini()
!call initialize_weight()

call initialize_expK()
call initialize_gamma_cont()

!call initialize_r250()

end subroutine 

subroutine initialize_weight()
use system_parameters
implicit none
integer:: i
do i=1,Nwalkers
weight(i)=1d0
!wscaling(i)=1d0
wphase(i)=0d0
impfunc0(i)=1d0
impfunc1(i)=1d0
enddo 
end subroutine

subroutine initialize_phi_ini
use system_parameters
use hubbard_param
implicit none

integer:: site,num, spin

phi_ini=0d0

select case(IniState)

case(1)

     do spin=1,2
      do num=1,numupdo(spin)
       phi_ini(num,num,spin) = 1d0
      enddo 
     enddo 
!------------------------------------------
case(2)

     do spin=1,2
      do site=1, Nsite
      do num=1,numupdo(spin)
       phi_ini(site,num,spin) = 1d0/sqrt(dble(Nsite))
      enddo 
      enddo 
     enddo 

case(3)
!-----------------------------------------
     spin=2
      do num=1,numupdo(spin)
       phi_ini(num,num,spin) = 1d0
      enddo 
   
     spin=1 
      do num=1,numupdo(spin)
       phi_ini(Nsite-num+1,num,spin) = 1d0
      enddo 
end select

end subroutine


subroutine initialize_phi
use system_parameters
use hubbard_param
implicit none

integer:: i,dim

!   phiT=phi_ini ! reset the trial wf to the initial wf

   do i=1,Nwalkers
       phi(:,:,:,i) = phi_ini
    enddo 


dim=nstate(nup,ndo)
call slater2fock(nup,ndo,dim,phi_ini(1:Nsite,1:nup,1),phi_ini(1:Nsite,1:ndo,2),phi_time)

end subroutine


subroutine initialize_expK
use system_parameters
use hubbard_param
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

!subroutine initialize_gamma
!use system_parameters
!real(kind=8):: y, res

!y=exp(dt*U/2d0)

!res=y**2-1d0 

!     if (res<0d0 ) then 
!        print *, 'error in initialize gamma'
!      endif

!      gamma=log(y+sqrt(res))
!print *, 'gamma', gamma
!end subroutine 

subroutine initialize_gamma_cont
use system_parameters
use hubbard_param
implicit none

real(kind=8):: rdummy

rdummy=sqrt(abs(U)*dt/2d0) 

if (U>0d0) then 
gamma=rdummy -dble(HS)*Xi*rdummy
else
gamma=rdummy +dble(HS)*Xi*rdummy
endif 

print *, 'Cont. HS tansformation, gamma=', gamma
end subroutine 


subroutine initialize_gamma_disc
use system_parameters
use hubbard_param
implicit none

complex(kind=8):: y
complex(kind=8),external:: csqrt, clog

complex(kind=8):: cdummy


      if (abs(U)<1d-10) then 
        gamma=0d0
        return
       endif

y=cos(dt*U/2d0)+Xi*sin(dt*U/2d0)

gamma=clog(y-csqrt(y**2-1d0))

!print *, 'gamma', gamma !,clog(y+csqrt(y**2-1d0))
!print *, y**2-csqrt(y**2-1d0)**2
!print *, (exp(gamma)+exp(-gamma))/2d0, exp(Xi*dt*U/2d0)

cdummy=(exp(gamma)+exp(-gamma))/2d0
      if (abs(cdummy-y)<1d-10) then
        print *, 'gamma checked!', gamma
      else
        print *, 'wrong gamma!'
        stop
      endif 

end subroutine 


function csqrt(z)
implicit none

complex(kind=8),parameter:: Xi=cmplx(0d0,1d0)

complex(kind=8),intent(IN):: z
complex(kind=8):: csqrt

real(kind=8):: rho, theta
real(kind=8):: re, im
real(kind=8):: pi
pi=dacos(-1d0)

re=dble(z)
im=dimag(z)

rho=dsqrt(re**2+im**2)
theta=dacos(re/rho)

rho=sqrt(rho)
theta=theta/2d0

csqrt=rho*dcos(theta)+rho*dsin(theta)*Xi

end function 


function clog(z)
implicit none

complex(kind=8),parameter:: Xi=cmplx(0d0,1d0)

complex(kind=8),intent(IN):: z
complex(kind=8):: clog

real(kind=8):: rho, theta
real(kind=8):: re, im

re=dble(z)
im=dimag(z)

rho=abs(z)
theta=datan(im/re)

if (isnan(theta)) then 
print *, 'wrong theta in clog',z
stop
endif 

clog=log(rho)+Xi*theta
!print *, rho, theta, clog
end function 



