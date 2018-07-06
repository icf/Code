program main
use system_parameters
use hubbard_param
use io_module
use lanc_param
use calphy_module
implicit none

complex(kind=8),allocatable:: gf(:,:,:), expKHF(:,:,:)
real(kind=8),allocatable:: occ_exact(:,:)
real(kind=8):: dbocc_exact, kin_exact
integer:: step,sample,spin, i
!integer:: output_counter
real(kind=8):: time, ln_wscaling, meanw,maxw, minw
character(len=200):: BaseName,OccName, OverlapName, ShiftName
!logical:: alive
!complex(kind=8):: dbocc
!complex(kind=8):: denominator

complex(kind=8):: cdummy


integer:: direction 
logical:: reconfig

real(kind=8),external:: phaseless_proj
complex(kind=8),allocatable:: xbar(:,:),x2bar(:,:) ! Nsite, Nwalker
complex(kind=8),allocatable:: expfactor(:) ! Nwalker
real(kind=8):: rdummy

complex(kind=8), external :: clog

!real(kind=8):: Hlanc(lanc_size,lanc_size)
!complex(kind=8):: ExpHlanc(lanc_size,lanc_size)
!real(kind=8):: Eigenvalue(lanc_size)
integer:: dim
complex(kind=8), allocatable:: phimat(:,:), phi2hf(:)
complex(kind=8):: overlap_mix, dbocc_mix, kin_mix
complex(kind=8),allocatable:: occ_mix(:,:)

complex(kind=8):: obsvr(3), deno !, obsvl(2)

!complex(kind=8),allocatable::Q(:,:),R(:,:)

real(kind=8):: fit_error

!complex(kind=8),allocatable:: phi_last(:)
!real(kind=8):: delta_theta


call initialize()

allocate(expKHF(Nsite,Nsite,2))
allocate(gf(Nsite,Nsite,2))

allocate(xbar(Nsite,Nwalkers))
allocate(x2bar(Nsite,Nwalkers))
allocate(expfactor(Nwalkers))
allocate(occ_exact(Nsite,2))
allocate(occ_mix(Nsite,2))


dim=nstate(nup, ndo)
allocate(phimat(lanc_size,dim))
allocate(phi2hf(dim))
!allocate(phi_last(dim))


call createFileName(BaseName,Dir)
call appendBaseName(BaseName,'tcpqmc_')

call appendBaseName(BaseName,'L',Nsite)
call appendBaseName(BaseName,'nup',nup)
call appendBaseName(BaseName,'ndo',ndo)

call appendBaseName(BaseName,'IniState',IniState)


call appendBaseName(BaseName,'Thop',2,Thop)
call appendBaseName(BaseName,'U',2,U)
call appendBaseName(BaseName,'Mu',2,Mu)

if (Phaseless) then 
call appendBaseName(BaseName,'Phaseless')
endif

call appendBaseName(BaseName,'Shift', Shift)


if (Local_energy) then 
call appendBaseName(BaseName,'LocalE')
endif


call appendBaseName(BaseName,'dt',4,dt)
call appendBaseName(BaseName,'TO',TrotterOrder)
call appendBaseName(BaseName,'Nwalkers',Nwalkers)
call appendBaseName(BaseName,'Maxshift',2,Maxshift)
call appendBaseName(BaseName,'UpperBound',2,UpperBound)
call appendBaseName(BaseName,'LowerBound',2,LowerBound)
call appendBaseName(BaseName,'PopCtrlSteps',PopContrlSteps)


call copyName(BaseName,OccName)
call copyName(BaseName,OverlapName)
call copyName(BaseName,ShiftName)


call appendBaseName(OccName,'_Occ.dat')
call appendBaseName(OverlapName,'_Overlap.dat')
call appendBaseName(ShiftName, '_Shift.dat')

!inquire(file=OccName, exist=alive)
!if (alive) then 
!call openUnit(OccName,11,'A')
!print *, 'file exist, append'
!else
!call openUnit(OccName,11,'N')
!print *, 'file not exist, create new'
!endif

call openUnit(OccName,11)
!call openUnit(OverlapName,13)

do sample=1,Nsamples
print *, 'sample', sample


!call openUnit(ShiftName,99)
do Totalsteps=2*Maxsteps,2*Maxsteps,2


call initialize_phi()
call initialize_weight()
!call initialize_mats()

ln_wscaling=0d0
do step=1,Totalsteps

time= (step)*dt

!print *, step, Totalsteps

!if (step<=TotalSteps/2) then 
direction=1 ! forward proprogation 
!else 
!direction=-1 ! backward proprogation 
!end if 

!care about the order , first perform the mc step 

call propagation_onestep1(step,direction,xbar,x2bar, expfactor,fit_error)

call ed_onestep(direction, dim, phi_time)


!phi_last=phi_time
!call lanczos_onestep(step,direction,matsize,dim,nup,ndo,dt,Hlanc,Eigenvalue, phimat,phi_time)

!overlap_mix=0d0
!do i=1,dim
!overlap_mix=overlap_mix+ conjg(phi_last(i))*phi_time(i)
!enddo 

!delta_theta= dimag(clog(overlap_mix))
!print *, phi_time
!print *, phi_last
!print *, overlap_mix
!print *, delta_theta
!pause

!--calculate physical quantities on phi_time
call cal_den(dim, nup, ndo,phi_time, occ_exact(:,1),occ_exact(:,2),dbocc_exact,kin_exact)


     do i=1, Nwalkers

call slater2fock(nup,ndo,dim,phi(1:Nsite,1:nup,1,i),phi(1:Nsite,1:ndo,2,i),phi2hf)

call cal_den_mix(dim, nup, ndo,phi_time,phi2hf,occ_mix,dbocc_mix, overlap_mix,kin_mix)

!print *, 'phi_time', phi_time
!print *, 'phi2hf', phi2hf

!if (step==Totalsteps) then 
!        print *, 'phi_time', Totalsteps, phi_time 
!print *, step,Totalsteps, abs(overlap_mix)
!pause
!pause
!endif 


 
!        impfunc1(i)=overlap(1)*overlap(2)
         impfunc1(i)=overlap_mix
         !print *, phi_time
         !print *, phi2hf
         !print *, overlap_mix
         !pause

!          cdummy=(impfunc1(i)/impfunc0(i)) 
!          write(99,'(i5,20f15.6)') step, dble(clog(cdummy)), dimag(clog(cdummy)), dble(expfactor(i)),  dimag(expfactor(i)), abs(xbar(1,i)) , abs((impfunc1(i)/impfunc0(i))*exp(expfactor(i))), weight(i), abs(impfunc1(i)), dble(impfunc0(i)), dimag(impfunc0(i)), dble(impfunc1(i)), dimag(impfunc1(i)), abs(impfunc0(i)), abs(impfunc1(i))


!write(99,'(i4,10f15.6)') step, occ_exact(1,1),dble(overlap_mix),dimag(overlap_mix), abs(overlap_mix)

          if (Local_energy) then 
           cdummy= exp(expfactor(i))
          else
           cdummy=(impfunc1(i)/impfunc0(i)) * exp(expfactor(i))
!          cdummy=dble(impfunc1(i))/dble(impfunc0(i)) * exp(expfactor(i))
          endif

           rdummy=abs(cdummy)
           if (isnan(rdummy)) then 
            print *, 'rdummy nan', cdummy, abs(cdummy)
            print *, expfactor(i), impfunc0(i), impfunc1(i)
            stop
           endif
!--------phaseless-------------------------------
!if (Phaseless) then 
!     rdummy=rdummy*max(cos(dimag(clog(impfunc1(i)/impfunc0(i)))), 0d0)
!endif
!-----------------------------------------
          impfunc0(i)=impfunc1(i)
!-----------------------------------------
!control over the amplitude fluctuation  
!          if (rdummy>Upperbound) then 
          !do spin=1,2
          !phi(1:Nsite,1:numupdo(spin),spin,i)=phi(1:Nsite,1:numupdo(spin),spin,i)*(sqrt(Upperbound/rdummy))**(1/dble(numupdo(spin)))
          !enddo 
!         impfunc0(i)=impfunc0(i) * Upperbound/rdummy
          ! wscaling(i)=wscaling(i)* rdummy/Upperbound
!          rdummy=Upperbound
!          endif

!          if (rdummy<Lowerbound) then 
          !wscaling(i)=wscaling(i)* rdummy/Lowerbound
          if (rdummy<1d-10) then  ! for zero W factor, the phase is no defined 
                 cdummy=Lowerbound
          endif
!          rdummy=Lowerbound
!          endif

    weight(i)=weight(i)*rdummy


     if (.not. Phaseless) then 
          wphase(i)=wphase(i)+ dimag(clog(cdummy))
     endif

!-----------------------------------------
!control over the amplitude fluctuation  
          if (weight(i)>Upperbound) then 
          weight(i)=Upperbound
          endif

          if (weight(i)<Lowerbound) then 
          weight(i)=Lowerbound
          endif

    enddo 


if (mod(step, PopContrlSteps)==0)then 
call assess_walkers(reconfig,meanw,maxw,minw)
!if (reconfig) then 
ln_wscaling=ln_wscaling + log(meanw)
call population_control(step)
!endif
endif



!do spin=1,2
!n=Nsite
!m=numupdo(spin)
!allocate(Q(N,N))
!allocate(R(N,M))

!  do i=1, Nwalkers
!  call zqr(N,M,phi(1:N,1:M,spin,i),Q,R)
 
! phi(1:N,1:M,spin,i)=Q(1:N,1:M)

! cdummy=1d0
! do j=1,min(m,n)
!  cdummy=cdummy*R(j,j)
! enddo 

! impfunc0(i)=impfunc0(i)/abs(cdummy)
!enddo 
!deallocate(Q)
!deallocate(R)
!enddo 


obsvr=0d0
deno=0d0
!------Normalize the wave function for each walker--------
do i=1,Nwalkers
    do spin=1,2
   call cal_overlap(Nsite,numupdo(spin), phi(1:Nsite,1:numupdo(spin),spin,i),phi(1:Nsite,1:numupdo(spin),spin,i),cdummy)

!weight(i)=weight(i)*sqrt(abs(cdummy))

impfunc0(i)=impfunc0(i)/(sqrt(abs(cdummy)))**(1/dble(numupdo(spin)))

phi(1:Nsite,1:numupdo(spin),spin,i)=phi(1:Nsite,1:numupdo(spin),spin,i)/(sqrt(abs(cdummy)))**(1/dble(numupdo(spin)))

 enddo 

call slater2fock(nup,ndo,dim,phi(1:Nsite,1:nup,1,i),phi(1:Nsite,1:ndo,2,i),phi2hf)
call cal_den_mix(dim, nup, ndo,phi_time,phi2hf,occ_mix,dbocc_mix, overlap_mix,kin_mix)
obsvr(1)=obsvr(1)+(occ_mix(1,1))*weight(i)* exp(Xi*wphase(i))
obsvr(2)=obsvr(2)+dbocc_mix*weight(i)* exp(Xi*wphase(i))
obsvr(3)=obsvr(3)+kin_mix*weight(i)* exp(Xi*wphase(i))
deno=deno+weight(i)* exp(Xi*wphase(i))

!print *, 'energy', dbocc_mix*U+kin_mix
!pause
enddo 


!print *, 'TotalSteps',Totalsteps
!pause

!do spin=1,2
!call cal_green(Nsite,numupdo(spin),phiHF(1:Nsite,1:numupdo(spin),spin,TotalSteps/2),phiHF(1:Nsite,1:numupdo(spin),spin,TotalSteps/2),gf(1:Nsite,1:Nsite,spin))
!call cal_occ(Nsite,gf(1:Nsite,1:Nsite,spin),occT(:,spin))
!enddo 


!deallocate(Bmats)
!deallocate(PhiHF)
!deallocate(occHF)
!deallocate(PhiHF_persite)

obsvr=obsvr/deno

!print *, time, abs(denominator)
write(11,'(i4,20f17.6)') step,time,occ_exact(1,1),dbocc_exact,kin_exact,dble(obsvr(:)),abs(deno),abs(deno)*dexp(ln_wscaling)/dble(Nwalkers), maxw, minw, meanw, fit_error


enddo  ! step
enddo  ! Totalsteps

enddo  ! Samples
close(11)
!close(99)

end program 
