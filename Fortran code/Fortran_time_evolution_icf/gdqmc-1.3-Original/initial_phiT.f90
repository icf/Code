
subroutine get_phiT()
use phiT_param
implicit none
if(PT.eq.0) then
  call read_phiT()
else if(PT.eq.1) then
  !Set from the FS wave function.
  call set_FS_phiT()
end if
end subroutine get_phiT


!Read PhiT from file. 1.Read the Dtot 2. Read the phiT and coe_multi 
subroutine read_phiT()
use phiT_param
use lattice_param
use model_param
use io_module
use mpi_serial_param
implicit none
integer::i,j,k
character(len=300)::filen
!Get Dtot,coe_multi:
open(unit=10,file='phiT_coe.dat',status='old')
  read(10,*) Dtot
  call allocate_phiT()
  do i=1,Dtot,1
     read(10,*) coe_multi(i)
  end do
close(10)
call norm_array(Dtot,coe_multi(1))
!Read phiT
open(unit=10,file='phiT.dat',status='old')
do i=1,Dtot,1
  do k=1,Ntot,1
     do j=1,2*Nsite,1
        read(10,*) phiT(j,k,i)
     end do
  end do
end do
close(10)

!call createFileName(filen,'./')
!call appendBaseName(filen,'l',Nsite)
!call appendBaseName(filen,'nt',Ntot)
!call appendBaseName(filen,'.dat')
!call openUnit(filen,10,'B')
!     read(10) Dtot
!     call allocate_phiT()
!     read(10) coe_multi(1:Dtot)
!     do i=1,Dtot,1
!        read(10) phiT(1:2*Nsite,1:Ntot,i)
!     end do
!close(10)
if(rank.eq.0) then  
  write(*,*) "Dtot of phiT:",Dtot
end if
!write(*,*) Dtot,coe_multi;call mystop
end subroutine read_phiT


!Set the phiT from free electrons wave function.
subroutine set_FS_phiT()
use phiT_param
use lattice_param
use mpi_serial_param
use model_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

Dtot=1
call allocate_phiT()
coe_multi(1)=1.d0
if(rank.eq.0) then 
  if(dtype.EQ.'c') then
    call input_phiT(2*Nsite,Ntot,Hzero,phiT(1,1,1))
  else if(dtype.EQ.'d') then
    call input_phiT_d()
  end if
end if
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call MPI_BCAST(phiT(1,1,1),1,phtype,0,MPI_COMM_WORLD,IERR)
#endif

end subroutine set_FS_phiT



!Diagonialize h0, give the lowest eigenstate to ph
subroutine input_phiT(nl,nt,h0,ph)
implicit none
integer,intent(IN)::nl
integer,intent(IN)::nt
complex(kind=8),intent(IN)::h0(nl,nl)
complex(kind=8),intent(OUT)::ph(nl,nt)
complex(kind=8)::hu(nl,nl)
real(kind=8)::ev(nl)


call check_Hermite_c(h0,nl)
call zcopy(nl*nl,h0,1,hu,1)
call eigen(Hu,nl,ev)
call zcopy(nl*nt,hu,1,ph,1)
end subroutine input_phiT


!Used to call input_phiT in decouple condition
subroutine input_phiT_d()
use phiT_param
use lattice_param
use mpi_serial_param
use model_param
implicit none
complex(kind=8)::h0(Nsite,Nsite)
complex(kind=8)::ph(Nsite,Ntot)
!complex(kind=8),allocatable::h0(:,:)
!complex(kind=8),allocatable::ph(:,:)
integer::i,j

!allocate(h0(Nsite,Nsite),ph(Nsite,Ntot))

h0(1:Nsite,1:Nsite)=Hzero(1:Nsite,1:Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),Hzero(1,1),2*Nsite,h0(1,1),Nsite)
call input_phiT(Nsite,Nspin(1),h0,ph(1,1))
phiT(1:Nsite,1:Nspin(1),1)=ph(1:Nsite,1:Nspin(1))
!call mkl_zomatcopy('C','N',Nsite,Nspin(1),(1.d0,0.d0),ph(1,1),Nsite,phiT(1,1,1),2*Nsite)


h0(1:Nsite,1:Nsite)=Hzero((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),Hzero(Nsite+1,Nsite+1),2*Nsite,h0(1,1),Nsite)
call input_phiT(Nsite,Nspin(2),h0,ph(1,1))
phiT((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,1)=ph(1:Nsite,1:Nspin(2))
!call mkl_zomatcopy('C','N',Nsite,Nspin(2),(1.d0,0.d0),ph(1,1),Nsite,phiT(Nsite+1,Nspin(1)+1,1),2*Nsite)

!do i=1,Nsite,1
!   do j=1,Nspin(2),1
!      write(*,*) ph(i,j),phiT(Nsite+i,Nspin(1)+j,1)
!   end do
!end do
!deallocate(h0,ph)
end subroutine input_phiT_d
