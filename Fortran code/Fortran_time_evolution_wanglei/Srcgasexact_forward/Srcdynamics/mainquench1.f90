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

real(kind=8),allocatable:: gsphi(:)
complex(kind=8),allocatable:: phi_time(:), phi(:,:)

real(kind=8):: Hlanc(lanc_size,lanc_size)
complex(kind=8):: ExpHlanc(lanc_size,lanc_size)
real(kind=8):: Eigenvalue(lanc_size)
real(kind=8):: time,  res

real(kind=8),allocatable:: denup(:), dendo(:)
real(kind=8),allocatable:: occ_hf(:,:),occ_mix(:,:)
real(kind=8):: dbocc_mix

integer,allocatable:: code(:)
real(kind=8):: comup,comdo, dbocc
character(len=132):: BaseName,RadName, WfName
integer:: i,j,k
integer:: matsize
real(kind=8)::  Uin
integer:: step, diff

integer(kind=4),external:: large_pickup, array2integer
integer:: maxnum, spin


NAMELIST /SystemSettings/ Nsite,nup,ndo,Maxsteps,stepforcycle,dt,Dir
NAMELIST /FHParams/  Thop,Uin,Mu,Vext, stagger_h

U=0d0

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


 !--- assign the hopping matrix ----
     Tmn=0.0_Rkind
      do ns=1, Nsite-1 
      Tmn(ns, ns+1)=1.0_RKind
      enddo
      !---------PBC----------------
             Tmn(1,Nsite)=1.0_RKind
       
  do ms=1, Nsite
      do ns=1, ms-1   ! ns< ms 
     Tmn(ms,ns)=Tmn(ns,ms)
      enddo
  enddo

  Tmn=Tmn*Thop
  

call numstate1()    !produce Nstate(nup, ndo)
!print *, nstate
!=======================openmp initialize ==========================
!$ call omp_set_num_threads(MAX_threads)

!$omp parallel private(myid)
!$  myid = OMP_GET_THREAD_NUM()
!$  num_threads = OMP_GET_NUM_THREADS()
!$  print *,'threads', myid, 'of',  num_threads,  'is alive'
!$omp end parallel  

!$ print *,'============================================='
!=========================================================


dim=nstate(nup, ndo)
allocate(gsphi(dim))
allocate(phi_time(dim))
allocate(phi(lanc_size,dim))



!--setup the initial state 
!--It is a ground state 
!call sub_lanc(nup,ndo,dim,gsenergy,gsphi)
!print * ,nup,ndo,nstate(nup,ndo), gsenergy

!--Or it is set up by hand

allocate(Table(Nstate(nup, ndo)))

call buildbasis(nup, ndo) ! build Table and invTable 

code=0
do i=1,nup
code(2*i-1)=1
enddo 
do i=1,ndo
!code(2*i)=1
code(2*Nsite-2*(i-1))=1
enddo 
!do j=1,dim
! call findbas(code, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej
! print * ,code
!enddo 
j=invTable(array2integer(code))   ! array rep-->demical rep-----> real inv table--->sequentail rep in subspace
print *, 'j',invTable(array2integer(code))

gsphi=0d0
gsphi(j)=1d0
deallocate(Table)
!--------

U=Uin
!---Copy initial state----------
do k=1,dim
phi_time(k)=gsphi(k)
enddo 
!----perform lanczos on a give target state phi_time
call sub_lancdynamics(nup, ndo, dim,phi_time,matsize, Hlanc,Phi)
call reigen(Hlanc(1:matsize,1:matsize),matsize,eigenvalue(1:matsize))


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

if (step<=Maxsteps) then 
call lanczos_onestep(step,1,stepforcycle,matsize,dim,nup,ndo,dt,Hlanc,Eigenvalue, phi,phi_time)
else
call lanczos_onestep(step,-1,stepforcycle,matsize,dim,nup,ndo,dt,Hlanc,Eigenvalue, phi,phi_time)
endif 


write(12,*) step, phi_time


enddo  ! end of time evolution 
close(11)
close(12)

print *, 'phi_time', phi_time
print *, 'gsphi', gsphi

deallocate(gsphi)
deallocate(phi_time)
deallocate(phi)

deallocate(Tmn)
deallocate(Nstate)
deallocate(denup)
deallocate(dendo)

end program main
