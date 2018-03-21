
subroutine get_phiT()
use phiT_param
use sc_loop_param
implicit none
if(PT.eq.0) then
  call read_phiT()
elseif(PT.eq.1) then
     if(sc_loop_flag .EQ. 0) then
        !Set from the FS wave function.
        call set_FS_phiT()
     elseif(sc_loop_flag .EQ. 1)then
        write(*,*)'get phiT from cc'
        call get_cc_phiT()
    endif
endif
    
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

!--------------------------------
!get phiT from cc
!--------------------------------
!Set the phiT from free electrons wave function.
subroutine get_cc_phiT()
use phiT_param
use lattice_param
use mpi_serial_param
use model_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

call deallocate_phiT()
call allocate_phiT()

if(rank.eq.0) then 
  if(dtype.EQ.'c') then
    write(*,*)'Error: sc_loop only apply to dtype=d'
    call mystop
  else if(dtype.EQ.'d') then
    call input_phiT_d_cc()
  end if
end if
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call MPI_BCAST(phiT(1,1,1),1,phtype,0,MPI_COMM_WORLD,IERR)
#endif

end subroutine get_cc_phiT
    
    
subroutine input_phiT_d_cc()
use phiT_param
use lattice_param
use mpi_serial_param
use model_param
use sc_loop_param
implicit none

call cc_decomposition
phiT(1:Nsite,1:Nspin(1),1:Dtot)=phi_sc(1:Nsite,1:Nspin(1),1:Dtot)

phiT((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,1:Dtot)=phi_sc(1:Nsite,Nspin(1)+1:Ntot,1:Dtot)

end subroutine input_phiT_d_cc
    
!------------------------------
!get phi_sc from |phi_sc><phi_sc|=cicj_sc_global
!------------------------------
subroutine cc_decomposition()
use phiT_param
use lattice_param
use model_param
use sc_loop_param
implicit none
complex(kind=8)::cicj_sc_global_up(Nsite,Nsite)
complex(kind=8)::cicj_sc_global_dn(Nsite,Nsite)
real(kind=8)::ev_up(Nsite),ev_dn(Nsite)

integer::i

cicj_sc_global_up(1:Nsite,1:Nsite)=cicj_sc_global(1:Nsite,1:Nsite)
cicj_sc_global_up=(conjg(transpose(cicj_sc_global_up))+cicj_sc_global_up)/2.0
cicj_sc_global_dn(1:Nsite,1:Nsite)=cicj_sc_global(Nsite+1:2*Nsite,Nsite+1:2*Nsite)
cicj_sc_global_dn=(conjg(transpose(cicj_sc_global_dn))+cicj_sc_global_dn)/2.0

call eigen(cicj_sc_global_up,Nsite,ev_up)
call eigen(cicj_sc_global_dn,Nsite,ev_dn)

!check ev_up and ev_dn
do i=1,Nsite,1
    write(*,*)'This is ev_up from cicj_sc_global_up, the order is',i,'the value is',ev_up(i)
enddo
do i=1,Nsite,1
    write(*,*)'This is ev_dn from cicj_sc_global_dn, the order is',i,'the value is',ev_dn(i)
enddo


do i=1,Dtot,1
   phi_sc(1:Nsite,1:Nspin(1),i)=cicj_sc_global_up(1:Nsite,Nsite-i*Nspin(1)+1:Nsite-(i-1)*Nspin(1))
   phi_sc(1:Nsite,Nspin(1)+1:Ntot,i)=cicj_sc_global_dn(1:Nsite,Nsite-i*Nspin(2)+1:Nsite-(i-1)*Nspin(2)) 
   coe_multi(i)=abs(sum(ev_up(Nsite-i*Nspin(1)+1:Nsite-(i-1)*Nspin(1)))+sum(ev_dn(Nsite-i*Nspin(2)+1:Nsite-(i-1)*Nspin(2))))
enddo

call norm_array(Dtot,coe_multi(1))
write(*,*)'coe_muti:',coe_multi

end subroutine cc_decomposition
