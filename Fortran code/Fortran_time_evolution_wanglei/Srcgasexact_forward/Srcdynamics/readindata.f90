program main
use prec
use hubbard_param
use lanc_param
use flags_param
use calphy_module
use io_module

implicit none

integer:: nup, ndo 
integer:: dim 
real(kind=RKind):: gsenergy(noflow)
integer::ns, ms 
integer:: alloc_error

complex(kind=8),allocatable:: phi_time(:)
real(kind=8):: time,  res

real(kind=8),allocatable:: denup(:), dendo(:)
real(kind=8),allocatable:: occ_hf(:,:),occ_mix(:,:)
real(kind=8):: dbocc_mix

integer,allocatable:: code(:)
real(kind=8):: comup,comdo, dbocc
character(len=132):: BaseName,RadName, WfName
integer:: i,j,k
integer:: matsize
integer:: step, diff, idummy

integer(kind=4),external:: large_pickup, array2integer
integer:: maxnum, spin


NAMELIST /SystemSettings/ Nsite,nup,ndo,Maxsteps,stepforcycle,dt,Dir
NAMELIST /FHParams/  Thop,U,Mu,Vext, stagger_h


OPEN(138,file='input')
READ(138,SystemSettings)
READ(138,FHParams)
CLOSE(138)

write(*, NML=SystemSettings)
write(*, NML=FHParams)


Nsmax=large_pickup(Nsite/2,Nsite)**2
print *, 'Maximum subspace size', Nsmax


if(invTable_flag==1) then 
alloc_error=0
allocate( invTable(0:2**(2*Nsite)-1), stat=alloc_error)
if(alloc_error/=0) then
print *, 'wrong in allocate invTable', alloc_error
stop
endif
endif

allocate(Tmn(Nsite,Nsite))
allocate(Nstate(0:Nsite,0:Nsite))
allocate(denup(Nsite))
allocate(dendo(Nsite))
allocate(occ_hf(Nsite,2))
allocate(code(2*Nsite))


call numstate1()    !produce Nstate(nup, ndo)

dim=nstate(nup, ndo)
allocate(phi_time(dim))



call createFileName(BaseName,Dir)
call appendBaseName(BaseName,'FHED_')
call appendBaseName(BaseName,'L',Nsite)
call appendBaseName(BaseName,'nup',nup)
call appendBaseName(BaseName,'ndo',ndo)

call appendBaseName(BaseName,'NNt',2,Thop)
call appendBaseName(BaseName,'U',2,U)
!call appendBaseName(BaseName,'Vext',2,Vext)
!if (abs(stagger_h)>1d-10) then 
!call appendBaseName(BaseName,'staggerh',2,stagger_h)
!endif 
!call appendBaseName(BaseName,'cycle',stepforcycle)


call copyName(BaseName,RadName)
call copyName(BaseName,WfName)


call appendBaseName(RadName,'Rad.dat')
call appendBaseName(WfName,'Wf.dat')


call openUnit(RadName,11)
call openUnit(WfName,12)


!--calculate physical quantities on phi_time
call cal_den(dim, nup, ndo,phi_time, denup,dendo,dbocc)

do ns=1,Nsite
print *, 'ns,nup,ndo', ns,denup(ns), dendo(ns)
enddo 
pause

do step =1,2*Maxsteps
time=(step-1)*dt

!--calculate physical quantities on phi_time
call cal_den(dim, nup, ndo,phi_time, denup,dendo,dbocc)

print *, 'time,norm', time, res
write(11,'(4f15.6)') time, denup(1),dendo(1),dbocc


!call slater2fock(nup,ndo,dim,phiT(1:Nsite,1:nup,2),phiT(1:Nsite,1:ndo,2),phiHF)

!call  cal_den_mix(dim, nup, ndo,phiHF,phi_time,occ_mix,dbocc_mix)


read(12,*) idummy, phi_time
print *, idummy, phi_time


enddo  ! end of time evolution 
close(11)
close(12)


deallocate(phi_time)

deallocate(Tmn)
deallocate(Nstate)
deallocate(denup)
deallocate(dendo)

end program main
