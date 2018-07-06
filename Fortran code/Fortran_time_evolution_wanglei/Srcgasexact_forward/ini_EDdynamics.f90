subroutine ini_EDdynamics
use prec
use hubbard_param
use lanc_param
use flags_param
use calphy_module
use io_module
use system_parameters
use hopt_module

implicit none

integer:: dim 
integer::ns, ms 
integer:: alloc_error

real(kind=8),allocatable:: gsphi(:) !, phimat(:,:)


real(kind=8):: Hlanc(lanc_size,lanc_size)
complex(kind=8):: ExpHlanc(lanc_size,lanc_size)
real(kind=8):: Eigenvalue(lanc_size)

real(kind=8),allocatable:: denup(:), dendo(:)
real(kind=8),allocatable:: occ_hf(:,:),occ_mix(:,:)

integer,allocatable:: code(:)
real(kind=8):: comup,comdo, dbocc
integer:: i,j,k
integer:: matsize
integer:: step, diff

integer(kind=4),external:: large_pickup, array2integer
integer:: maxnum, spin

real(kind=8), allocatable:: phi_old(:), phi_new(:)


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
!allocate(gsphi(dim))
allocate(phi_time(dim))
!allocate(phiexact(dim))
!allocate(phimat(lanc_size,dim))


!--setup the initial state 
!--It is a ground state 
!call sub_lanc(nup,ndo,dim,gsenergy,gsphi)
!print * ,nup,ndo,nstate(nup,ndo), gsenergy

!--Or it is set up by hand

!allocate(Table(Nstate(nup, ndo)))

!call buildbasis(nup, ndo) ! build Table and invTable 

!code=0
!do i=1,nup
!code(2*i-1)=1
!enddo 
!do i=1,ndo
!code(2*i)=1
!code(2*Nsite-2*(i-1))=1
!enddo 
!do j=1,dim
! call findbas(code, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej
! print * ,code
!enddo 
!j=invTable(array2integer(code))   ! array rep-->demical rep-----> real inv table--->sequentail rep in subspace
!print *, 'j',invTable(array2integer(code))

!gsphi=0d0
!gsphi(j)=1d0
!deallocate(Table)
!--------

!---Copy initial state----------
!do k=1,dim
!phi_time(k)=gsphi(k)
!enddo 
!----perform lanczos on a give target state phi_time
!call sub_lancdynamics(nup, ndo, dim,phi_time,matsize, Hlanc,phimat)
!call reigen(Hlanc(1:matsize,1:matsize),matsize,eigenvalue(1:matsize))


allocate(Ham(dim,dim))
allocate(ExpHam(dim,dim))
allocate(EigenHam(dim))

!Ham(1,:)=(/U, Thop, -Thop, 0d0/)
!Ham(2,:)=(/Thop,0d0,0d0,Thop/)
!Ham(3,:)=(/-Thop,0d0,0d0,-Thop/)
!Ham(4,:)=(/0d0, Thop, -Thop, U/)


allocate(phi_old(dim))
allocate(phi_new(dim))

do j=1,dim
phi_old=0d0
phi_old(j)=1d0
call Hopt(dim, nup, ndo, Phi_old, Phi_new)
Ham(:,j)=phi_new(:)
print *,j, Ham(:,j)
enddo 

call reigen(Ham,dim,EigenHam)

ExpHam=0d0
DO i=1,dim
   DO j=1,dim
         DO k=1,dim
ExpHam(i,j) = ExpHam(i,j) + Ham(i,k)*exp(-Xi*dt*EigenHam(k))*Ham(j,k)
         END DO
   END DO
 END DO


end subroutine ini_EDdynamics
