
subroutine get_xbar(i,direction, xbar,x2bar)
use system_parameters
use hubbard_param
use calphy_module
implicit none

integer,intent(IN):: i , direction
complex(kind=8),intent(OUT):: xbar(Nsite) , x2bar(Nsite)
complex(kind=8):: overlap(2)

integer:: spin, site
complex(kind=8):: gf(Nsite, Nsite, 2)
complex(kind=8):: occup, occdo
complex(kind=8):: cdummy

integer:: dim
complex(kind=8), allocatable:: phi2hf(:)
complex(kind=8):: overlap_mix, dbocc_mix, kin_mix
complex(kind=8):: occ_mix(Nsite,2)

real(kind=8), external:: rannyu

dim=nstate(nup,ndo)
allocate(phi2hf(dim))

if (direction==1) then 
cdummy=gamma
else if (direction==-1) then 
cdummy=conjg(gamma)
else 
print *, 'wrong gamma'
stop
endif

select case(Shift)

case(0) ! noshift
do site=1, Nsite
xbar(site) = 0d0
enddo 
x2bar=0d0

case(3) ! random shift
do site=1, Nsite
xbar(site) = (rannyu() + rannyu() *Xi ) * cdummy
enddo 
!xbar = (rannyu() + rannyu() *Xi ) * cdummy
x2bar=0d0


case(1) ! mixed shift 

call slater2fock(nup,ndo,dim,phi(1:Nsite,1:nup,1,i),phi(1:Nsite,1:ndo,2,i),phi2hf)
call  cal_den_mix(dim, nup, ndo,phi_time,phi2hf,occ_mix,dbocc_mix, overlap_mix, kin_mix)

do site=1, Nsite
xbar(site)= -cdummy*(occ_mix(site,1)-occ_mix(site,2))
x2bar(site)=0d0
enddo 

case(4) ! mixed shift without the denominator

call slater2fock(nup,ndo,dim,phi(1:Nsite,1:nup,1,i),phi(1:Nsite,1:ndo,2,i),phi2hf)
call  cal_den_mix(dim, nup, ndo,phi_time,phi2hf,occ_mix,dbocc_mix, overlap_mix, kin_mix)

do site=1, Nsite
xbar(site)= -cdummy*(occ_mix(site,1)-occ_mix(site,2)) *overlap_mix
x2bar(site)=0d0
enddo 

case(5) ! shift with dt correct, only for TO=1

call slater2fock(nup,ndo,dim,phi(1:Nsite,1:nup,1,i),phi(1:Nsite,1:ndo,2,i),phi2hf)
call  cal_den_mix(dim, nup, ndo,phi_time,phi2hf,occ_mix,dbocc_mix, overlap_mix, kin_mix)

do site=1, Nsite
xbar(site)= -cdummy*(occ_mix(site,1)-occ_mix(site,2)) *overlap_mix/impfunc0(i)
x2bar(site)=0d0
enddo 


case(2) !self shift
!----Calculate HF Green's function------
do spin=1,2
!    call cal_overlap(Nsite,numupdo(spin),phiT(1:Nsite,1:numupdo(spin),spin),phi(1:Nsite,1:numupdo(spin),spin,i),overlap(spin))
!    call cal_green(Nsite,numupdo(spin),phiT(1:Nsite,1:numupdo(spin),spin),phi(1:Nsite,1:numupdo(spin),spin,i),gf(1:Nsite,1:Nsite,spin))

    call cal_overlap(Nsite,numupdo(spin),phi(1:Nsite,1:numupdo(spin),spin,i),phi(1:Nsite,1:numupdo(spin),spin,i),overlap(spin))
    call cal_green(Nsite,numupdo(spin),phi(1:Nsite,1:numupdo(spin),spin,i),phi(1:Nsite,1:numupdo(spin),spin,i),gf(1:Nsite,1:Nsite,spin))

enddo 

do site=1, Nsite
occup=1d0-gf(site,site,1)
occdo=1d0-gf(site,site,2)

xbar(site)=-cdummy*(occup*overlap(1)+dble(HS)*occdo*overlap(2))
x2bar(site)=0d0
enddo 

end select

do site=1, Nsite
if (abs(xbar(site)) > Maxshift  ) then 
   xbar(site)= xbar(site)/abs(xbar(site)) * Maxshift
endif
end do 


if (isnan(dble(xbar(1)))) then 
print *, 'xbar nan'
stop
endif 

end subroutine  get_xbar



subroutine cal_expKHF(direction, occ,expKHF)
use system_parameters
use hubbard_param
implicit none

integer:: direction 
real(kind=8),intent(IN):: occ(Nsite,2)
complex(kind=8),intent(OUT):: expKHF(Nsite,Nsite,2)

real(kind=8):: KHFmat(Nsite,Nsite), Eigenvalue(Nsite)
integer:: i,j,k, spin

ExpKHF=0d0

  do spin=1,2
     KHFmat=Kmat
    do i=1,Nsite
      KHFmat(i,i)=Kmat(i,i) !+U*occ(i,3-spin)
    enddo 

call reigen(KHFmat,Nsite,Eigenvalue)

 DO i=1,Nsite
   DO j=1,Nsite
         DO k=1,Nsite
ExpKHF(i,j,spin) = ExpKHF(i,j,spin) + KHFmat(i,k)*exp(-Xi*dble(direction)*dt*EigenValue(k))*KHFmat(j,k)
         END DO
   END DO
 END DO
 enddo 
end subroutine 


subroutine cal_overlap(n,m,phi_l,phi_r,overlap)
implicit none

integer,intent(IN):: n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: overlap

complex(kind=8):: tempmat(m,m), phi_l_dag(m,n)

if (isnan(dble(phi_l(1,1)))) then 
print *, 'wrong mat to caldet'
stop
endif 

phi_l_dag=transpose(conjg(phi_l))

tempmat= matmul(phi_l_dag, phi_r)
call caldet_c(m,tempmat,overlap)

end subroutine 

subroutine cal_green(n,m,phi_l,phi_r,gf)
implicit none

integer,intent(IN):: n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: gf(n,n)

integer:: i,j
integer:: k,l
complex(kind=8):: tempmat(m,m)

tempmat= matmul(transpose(conjg(phi_l)), phi_r)
call  inverse(tempmat,m)


gf=0d0
    do i=1,n
      do j=1,n
         if (i==j) then 
            gf(i,j)=1d0
           else
            gf(i,j)=0d0
         endif

          do k=1,m
            do l=1,m
          gf(i,j)=gf(i,j) - phi_r(i,k)*tempmat(k,l) *conjg(phi_l(j,l))
            enddo 
         enddo 

       enddo 
      enddo 

end subroutine 


subroutine cal_occ(n,gf,occ)
implicit none
integer,intent(IN):: n
complex(kind=8),intent(IN):: gf(n,n)
real(kind=8), intent(OUT):: occ(n)

integer:: i

do i=1,n
occ(i)=1d0-dble(gf(i,i))
enddo 

end subroutine 


subroutine population_control(step)
use system_parameters
use hubbard_param
implicit none

integer,intent(IN):: step
integer:: itable(Nwalkers)

complex(kind=8),allocatable:: phi_temp(:,:,:,:) ! Nsite,num,spin,walker
complex(kind=8):: impfunc0_temp(Nwalkers)
real(kind=8):: wphase_temp(Nwalkers)
!real(kind=8):: wscaling_temp(Nwalkers)
integer:: i,j

allocate(phi_temp(Nsite,max(nup,ndo),2,Nwalkers))

phi_temp=phi
wphase_temp=wphase
!wscaling_temp=wscaling
impfunc0_temp=impfunc0

call  reconfiguration(Nwalkers,weight,itable)

do i=1,Nwalkers
j=itable(i)
phi(:,:,:,i)=phi_temp(:,:,:,j)
wphase(i)=wphase_temp(j)
!wscaling(i)=wscaling_temp(j)
impfunc0(i)=impfunc0_temp(j)
weight(i)=1d0
enddo 

deallocate(phi_temp)

end subroutine population_control 


subroutine assess_walkers(reconfig,meanw,maxw,minw)
use system_parameters
implicit none

real(kind=8)::  maxw, minw, meanw
integer:: i

logical,intent(OUT):: reconfig

meanw=0d0
do i=1,Nwalkers

if (weight(i)<0d0) then 
print *, 'negative weight', weight(i)
stop
endif
meanw=meanw+weight(i)
enddo
meanw=meanw/dble(Nwalkers)

maxw=maxval(weight(:),1)
minw=minval(weight(:),1)

!if (maxw-minw>1d0*meanw)  then 
if (maxw>2d0 .or. minw<0.5d0)  then 
!if (meanw>2d0 .or. meanw<0.5d0) then 
reconfig=.true.
else
reconfig=.false.
endif

!if (reconfig) then 
!print *, maxw,minw,meanw, reconfig
!pause
!endif

end subroutine assess_walkers
