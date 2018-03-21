
subroutine get_phiT()
use phiT_param
use lattice_param
use model_param
use io_module
use mpi_serial_param
implicit none
integer::i,j,k
if(PT.eq.0) then
  call read_phiT()
  if(I_wavefun.eq.2)then
    call read_F_matrix()
    if(Nzeta.gt.0)then
      call read_D_unpaired()
    endif
  endif
else if(PT.eq.1) then
  !Set from the FS wave function.
  call set_FS_phiT()
  if(I_wavefun.eq.2)then
    call read_F_matrix()
    if(Nzeta.gt.0)then
      call read_D_unpaired()
    endif
  endif
end if
!IN THE BCS CASE phiT is used only to build the initial walker


if(rank.eq.0)then
  open(2,file='psi_trial.info',status='unknown')
  do i=1,Dtot,1
    do k=1,Ntot,1
       do j=1,2*Nsite,1
          write(2,*) phiT(j,k,i)
       end do
    end do
  end do
  close(2)
  if(I_wavefun.eq.2)then
    open(2,file='pairing_matrix.info',status='unknown')
    do j=1,Nsite
      do i=1,Nsite
        write(2,*) FPairing(i,j)
      enddo
    enddo
    close(2)
    if(Nzeta.gt.0)then
      open(2,file='unpaired_orbitals.info',status='unknown')
      do k=1,Nzeta
        do i=1,Nsite
          write(2,*) DUnpaired(i,k)
        enddo
      enddo
      close(2)
    endif
  endif
endif


end subroutine get_phiT

subroutine read_F_matrix()
use param
use phiT_param
use lattice_param
use model_param
use io_module
use mpi_serial_param
implicit none
integer::i,j
dtype='d'
Dtot=1
coe_multi(1)=one
!call allocate_phiT()
open(unit=10,file='pairing.dat',status='old')
do j=1,Nsite,1
     do i=1,Nsite,1
        read(10,*) FPairing(i,j)
     end do
end do
close(10)
end subroutine read_F_matrix




subroutine read_D_unpaired()
use param
use phiT_param
use lattice_param
use model_param
use io_module
use mpi_serial_param
implicit none
integer::i,j
!call allocate_phiT()
open(unit=10,file='unpaired.dat',status='old')
do j=1,Nzeta,1
     do i=1,Nsite,1
        read(10,*) DUnpaired(i,j)
     end do
end do
close(10)
end subroutine read_D_unpaired

!Read PhiT from file. 1.Read the Dtot 2. Read the phiT and coe_multi 
subroutine read_phiT()
use param
use phiT_param
use lattice_param
use model_param
use io_module
use mpi_serial_param
implicit none
integer::i,j,k,j1,j2,npartu,npartd
real(kind=8)::eps
complex(kind=8)::dummy,dummy2
complex(kind=8),dimension(:,:,:),allocatable::phi_read
character(len=300)::filen

eps=1.d-3

!Get Dtot,coe_multi:
open(unit=10,file='phiT_coe.dat',status='old')
  read(10,*) Dtot
  if(I_wavefun.eq.2)then
    Dtot=1
  endif
  call allocate_phiT()
  do i=1,Dtot,1
     read(10,*) coe_multi(i)
  end do
close(10)
call norm_array(Dtot,coe_multi(1))
!Read phiT
allocate(phi_read(2*Nsite,Ntot,Dtot))
open(unit=10,file='phiT.dat',status='old')
do i=1,Dtot,1
  do k=1,Ntot,1
     do j=1,2*Nsite,1
        read(10,*) phi_read(j,k,i)
     end do
  end do
end do
!normalize
do i=1,Dtot,1
 do k=1,Ntot,1
  dummy=zero
  do j=1,2*Nsite,1
    dummy=dummy+conjg(phi_read(j,k,i))*phi_read(j,k,i)
  end do
  do j=1,2*Nsite,1
    phi_read(j,k,i)=phi_read(j,k,i)/dsqrt(dble(dummy))
  end do
 enddo
enddo
close(10)

!write(*,*)'Abbiamo letto phiT.dat '

if(I_wavefun.eq.2)then
  dtype='d'
  if(rank.eq.0)then
    write(*,*)
    write(*,*)'BCS formalism uses decoupled scheme'
    write(*,*)
  endif
endif
if(dtype.eq.'c')then
  phiT=phi_read
elseif(dtype.eq.'d'.or.dtype.eq.'m')then
  do k=1,Dtot,1
    npartu=0
    npartd=0
    do j2=1,Ntot
      dummy=zero
      do j1=1,Nsite,1
        dummy=dummy+conjg(phi_read(j1,j2,k))*phi_read(j1,j2,k)
      enddo
      if(abs(dummy-one).lt.eps)then
        npartu=npartu+1
        do j1=1,Nsite,1
           phiT(j1,npartu,k)=phi_read(j1,j2,k)
        enddo
        do j1=Nsite+1,2*Nsite,1
           phiT(j1,npartu,k)=zero
        enddo
      else
        dummy2=zero
        do j1=Nsite+1,2*Nsite,1
          dummy2=dummy2+conjg(phi_read(j1,j2,k))*phi_read(j1,j2,k)
        enddo
        if(abs(dummy2-one).lt.eps)then
          npartd=npartd+1
          do j1=1,Nsite,1
            phiT(j1,Nspin(1)+npartd,k)=zero
          enddo
          do j1=Nsite+1,2*Nsite,1
            phiT(j1,Nspin(1)+npartd,k)=phi_read(j1,j2,k)
          enddo
        else
          write(*,*)npartu,npartd,dummy2
          write(*,*)'Problems in building phiT'
          stop
        endif
      endif
    enddo     
  enddo
endif

deallocate(phi_read)

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

integer::i,j,k

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
  else
    if(rank.eq.0)write(*,*)'Please use c or d method'
    call mystop
  end if
end if

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call MPI_BCAST(phiT(1,1,1),1,phtype2,0,MPI_COMM_WORLD,IERR)
#endif

end subroutine set_FS_phiT



!Diagonialize h0, give the lowest eigenstate to ph
subroutine input_phiT(nl,nt,h0,ph)
use mpi_serial_param
implicit none
integer,intent(IN)::nl
integer,intent(IN)::nt
integer::i,j
complex(kind=8),intent(IN)::h0(nl,nl)
complex(kind=8),intent(OUT)::ph(nl,nt)
complex(kind=8)::hu(nl,nl)
real(kind=8)::ev(nl)

integer::icall
data icall/0/
save icall

call check_Hermite_c(h0,nl)
call zcopy(nl*nl,h0,1,hu,1)
call eigen(Hu,nl,ev)

if(rank.eq.0)then
  open(2,file='single_particle_basis.info',status='unknown',access='append')
  rewind(2)
  do i=1,nl
    do j=1,nl
      write(2,*)hu(j,i)
    enddo
    write(2,*)
    write(2,*)
  enddo  
  close(2)
  open(2,file='single_particle_levels.info',status='unknown',access='append')
  rewind(2) 
  do i=1,nl
    write(2,*)ev(i)
  enddo
  write(2,*)
  write(2,*)
  close(2)  
endif


call zcopy(nl*nt,hu,1,ph,1)
icall=icall+1
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

phiT=0.d0

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
