
subroutine get_phiT()
use phiT_param
use lattice_param
use model_param
use io_module
use mpi_serial_param

use sc_loop_param
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
  if(I_wavefun.eq.1)then
    if(sc_loop_flag .EQ. 0) then
       if(sc_initial_choose .EQ. 1)then
          !Set from the FS wave function.
          write(*,*)'get phiT from FS'
          call set_FS_phiT()
       elseif(sc_initial_choose .EQ. 2)then
          !Set from the HF wave function.
          write(*,*)'get phiT from UHF'
          call set_HF_phiT()
       elseif(sc_initial_choose .EQ. 3)then
          !Set from the HF wave function.
          write(*,*)'get phiT from GHF'
          call set_GHF_phiT()
       endif 
    endif
  endif   
  
  if(I_wavefun.eq.2)then  
     if(sc_loop_flag .EQ. 1)then
        write(*,*)'get phiT from cc'
        call get_cc_phiT()
        !-------------------------------------------------
        !trans the notation of BCS_sc to Fpairing
        !-------------------------------------------------
        FPairing=BCS_sc
        !-------------------------------------------------
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


!Set the phiT from UHF wave function.
subroutine set_HF_phiT()
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
    write(*,*)'still EF for couple UHF'
    call input_phiT(2*Nsite,Ntot,Hzero,phiT(1,1,1))
  else if(dtype.EQ.'d') then
    call input_phiT_d_HF()
  else
    if(rank.eq.0)write(*,*)'Please use c or d method'
    call mystop
  end if
end if

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call MPI_BCAST(phiT(1,1,1),1,phtype2,0,MPI_COMM_WORLD,IERR)
#endif

end subroutine set_HF_phiT


!Set the phiT from GHF wave function.
subroutine set_GHF_phiT()
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
    call input_phiT_GHF(2*Nsite,Ntot,Hzero,phiT(1,1,1))
  elseif(dtype .EQ. 'c') then
    if(rank.eq.0)write(*,*)'GHF only work in spin_couple situation'
    call mystop
  end if
end if

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call MPI_BCAST(phiT(1,1,1),1,phtype2,0,MPI_COMM_WORLD,IERR)
#endif

end subroutine set_GHF_phiT
    

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


!-------GHF-----
subroutine input_phiT_GHF()
use phiT_param
use lattice_param
use mpi_serial_param
use model_param

use xhf_param
implicit none
complex(kind=8)::ph_save(2*Nsite,Ntot)

call xhf(ph_save)

phiT(1:2*Nsite,1:Ntot,1)=ph_save(1:2*Nsite,1:Ntot)

end subroutine input_phiT_GHF


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
    
    
subroutine input_phiT_d_HF()
use phiT_param
use lattice_param
use mpi_serial_param
use model_param

use HF_param
implicit none
complex(kind=8)::phi_HF(2*Nsite,Ntot)
complex(kind=8)::n_old_in(2*Nsite),n_old_out(2*Nsite),n_new_in(2*Nsite)
complex(kind=8)::a
integer::i

a=HF_a
n_new_in=0
n_old_in=0
n_old_out=0
phi_HF=0
call initial_HF_n(n_new_in)

do i=1,HF_num
   n_old_in=n_new_in 
   call get_phi_from_HF(n_old_in,phi_HF)
   call get_n_from_phi(n_old_out,phi_HF)
   call get_new_n(a,n_old_in,n_old_out,n_new_in)
enddo

phiT=0
phiT(1:Nsite,1:Nspin(1),1)=phi_HF(1:Nsite,1:Nspin(1))
phiT((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,1)=phi_HF(Nsite+1:2*Nsite,Nspin(1)+1:Ntot)

end subroutine input_phiT_d_HF
    
    
!--------------------------------
!get phiT from cc
!--------------------------------
!Set the phiT from free electrons wave function.
subroutine get_cc_phiT()
use phiT_param
use lattice_param
use mpi_serial_param
use model_param

use sc_loop_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

call deallocate_phiT()
call allocate_phiT()
!--------------------
!if want to modefy coe_muti, must consider Bcast it to all rank
coe_multi(1)=1.d0
!--------------------
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
  call MPI_BCAST(BCS_sc,Nsite*Nsite,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
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

phiT=0
phiT(1:Nsite,1:Nspin(1),1:Dtot)=phi_sc(1:Nsite,1:Nspin(1),1:Dtot)
phiT((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,1:Dtot)=phi_sc(1:Nsite,Nspin(1)+1:Ntot,1:Dtot)

end subroutine input_phiT_d_cc
    
!------------------------------
!get phi_sc from |phi_sc><phi_sc|=cicj_sc_global
!------------------------------
subroutine cc_decomposition()
use param
use phiT_param
use lattice_param
use model_param
use sc_loop_param

use io_module
use method_param
implicit none
complex(kind=8)::cicj_sc_global_up(Nsite,Nsite)
complex(kind=8)::cicj_sc_global_dn(Nsite,Nsite)
complex(kind=8)::cicj_sc_global_up1(Nsite,Nsite)
complex(kind=8)::cicj_sc_global_dn1(Nsite,Nsite)

complex(kind=8)::temp_up(Nsite,Nsite)
complex(kind=8)::temp_dn(Nsite,Nsite)
complex(kind=8)::sum_up,sum_dn,sum_ev
real(kind=8)::ev_up(Nsite),ev_dn(Nsite)
complex(kind=8)::ev_m(Nsite,Nsite),det_ev_m(Nsite,Nsite)
complex(kind=8)::temp_BCS_sc(Nsite,Nsite)

integer::i,j
real(kind=8):: a
complex(kind=8)::alpha,delta
complex(kind=8)::b
!FE initial----------------
complex(kind=8)::h0(Nsite,Nsite)
complex(kind=8)::ph(Nsite,Ntot)
!var-------------
integer::var_num
real(kind=8)::var_gap
complex(kind=8)::var_step
!test value-----
complex(kind=8)::a1,b1,c1,sum1


if(sc_ite_flag==1) then
   a=0.75
elseif(sc_ite_flag==0) then
   a=0
endif

!------------------------------------------
!input initial wave function from Green Function cicj_sc_global
!------------------------------------------
if(GM_input_flag .EQ. 1)then
   write(*,*)'read from input cicj_sc_global'
   call read_cicj_sc_global()
   sc_initial_choose=0
endif
!------------------------------------------
!Modify Green Function cicj_sc_global
!Only for 4 by 4 system with PBC
!------------------------------------------
!call GF_modifier()
!write(*,*)'modified_cicj_sc_global:',cicj_sc_global

cicj_sc_global_up1(1:Nsite,1:Nsite)=cicj_sc_global1(1:Nsite,1:Nsite)
if (a .EQ. 0) then
   cicj_sc_global_up(1:Nsite,1:Nsite)=(1-a)*cicj_sc_global(1:Nsite,1:Nsite)
else
   cicj_sc_global_up(1:Nsite,1:Nsite)=(1-a)*cicj_sc_global(1:Nsite,1:Nsite)+a*cicj_sc_global_up1(1:Nsite,1:Nsite)
endif
cicj_sc_global_up=(conjg(transpose(cicj_sc_global_up))+cicj_sc_global_up)/2.0

cicj_sc_global_dn1(1:Nsite,1:Nsite)=cicj_sc_global1(Nsite+1:2*Nsite,Nsite+1:2*Nsite)
if(a .EQ. 0) then
   cicj_sc_global_dn(1:Nsite,1:Nsite)=(1-a)*cicj_sc_global(Nsite+1:2*Nsite,Nsite+1:2*Nsite)
else
   cicj_sc_global_dn(1:Nsite,1:Nsite)=(1-a)*cicj_sc_global(Nsite+1:2*Nsite,Nsite+1:2*Nsite)+a*cicj_sc_global_dn1(1:Nsite,1:Nsite)
endif
cicj_sc_global_dn=(conjg(transpose(cicj_sc_global_dn))+cicj_sc_global_dn)/2.0

cicj_sc_global1(1:Nsite,1:Nsite)=cicj_sc_global_up(1:Nsite,1:Nsite)
cicj_sc_global1(Nsite+1:2*Nsite,Nsite+1:2*Nsite)=cicj_sc_global_dn(1:Nsite,1:Nsite)

call eigen(cicj_sc_global_up,Nsite,ev_up)
call eigen(cicj_sc_global_dn,Nsite,ev_dn)

!----------------------------------------
!check 4*4 exact input degenercy problem
!cicj_sc_global_dn=cicj_sc_global_up
!ev_dn=ev_up
!check cc_decomposition
!----------------------------------------
!--------------------------------
 call openUnit('hubb_gdmc_eigenvalue_up.dat',81,'A')
 call openUnit('hubb_gdmc_eigenvalue_dn.dat',82,'A')
    write(81,*)' '
    write(82,*)' '
    write(81,*)'step: ',sc_step_counter
    write(82,*)'step: ',sc_step_counter
    write(81,*)' '
    write(82,*)' '
 do i=1,Nsite,1
    write(81,*)ev_up(i)
    write(82,*)ev_dn(i)
 enddo

 call openUnit('hubb_gdmc_eigenvector_up.dat',83,'A')
 call openUnit('hubb_gdmc_eigenvector_dn.dat',84,'A')

 write(83,*)' '
 write(84,*)' '
 write(83,*)'step: ',sc_step_counter
 write(84,*)'step: ',sc_step_counter
 write(83,*)' '
 write(84,*)' '
 do i=1,Nsite,1
    write(83,*)'eigenvector: ',i,'real part'
    write(84,*)'eigenvector: ',i,'real part'
    do j=1,Nsite,1
       write(83,*)dble(cicj_sc_global_up(j,i))
       write(84,*)dble(cicj_sc_global_dn(j,i))
    enddo
    write(83,*)'imag part'
    write(84,*)'imag part'
    do j=1,Nsite,1
       write(83,*)aimag(cicj_sc_global_up(j,i))
       write(84,*)aimag(cicj_sc_global_dn(j,i))
    enddo
 enddo 
 
 close(81)
 close(82)
 close(83)
 close(84)

!-----------------------------------------
!-----------------------------------------
!get BCS_SC
!--------
!get ev_m
if (decMethod .EQ. 2)then
  write(*,*)'Analytic'
  ev_m=0
  write(*,*)'ev_up'
  write(*,*)ev_up
  write(*,*)'ev_dn'
  write(*,*)ev_dn
  do i=1,Nsite,1
     b=(ev_up(i)+ev_dn(i))/2
     ev_m(i,i)=sqrt(b/(1-b))
     !write(*,*)'ev_m:',ev_m(i,i)
     write(*,*)real(b)
  enddo
endif
if (decMethod .EQ. 1)then
  write(*,*)'DET'
  ev_m=0
  do i=Nsite-Nspin(1)+1,Nsite,1
     ev_m(i,i)=1
  enddo
endif
!----------------------------------------
!----------------------------------------
!test---------------
!a1=ev_m(Nsite,Nsite)
!b1=ev_m(Nsite-1,Nsite-1)
!c1=ev_m(Nsite-Nspin(1)-1,Nsite-Nspin(1)-1)
!sum1=0
!write(*,*)'a1:',a1
!write(*,*)'b1:',b1
!write(*,*)'c1:',c1
!sum1=sum1+(6*(a1*b1**4*c1**2)**2+6*(a1*c1**4*b1**2)**2+16*(a1*b1**3*c1**3)**2)
!sum1=sum1+(3*(a1*b1**2*c1**4)**2+6*(a1*b1**4*c1**2)**2+3*(b1**3*c1**4)**2&
!&+4*(b1**4*c1**3)**2+12*(a1*b1**3*c1**3)**2)*4
!sum1=sum1+(3*(a1*c1**2*b1**4)**2+6*(a1*c1**4*b1**2)**2+3*(c1**3*b1**4)**2&
!&+4*(c1**4*b1**3)**2+12*(a1*b1**3*c1**3)**2)*4
!write(*,*)'a:',Nspin(1)*(6*(a1*b1**4*c1**2)**2+6*(a1*c1**4*b1**2)**2+16*(a1*b1**3*c1**3)**2)/sum1
!write(*,*)'b:',Nspin(1)*(3*(a1*b1**2*c1**4)**2+6*(a1*b1**4*c1**2)**2&
!&+3*(b1**3*c1**4)**2+4*(b1**4*c1**3)**2+12*(a1*b1**3*c1**3)**2)/sum1
!write(*,*)'c:',Nspin(1)*(3*(a1*c1**2*b1**4)**2+6*(a1*c1**4*b1**2)**2&
!&+3*(c1**3*b1**4)**2+4*(c1**4*b1**3)**2+12*(a1*b1**3*c1**3)**2)/sum1

temp_BCS_sc=zero
BCS_sc=zero
call zgemm('N','N',Nsite,Nsite,Nsite,one,cicj_sc_global_up(1,1),Nsite,ev_m(1,1),Nsite,zero,temp_BCS_sc(1,1),Nsite)
call zgemm('N','T',Nsite,Nsite,Nsite,one,temp_BCS_sc(1,1),Nsite,cicj_sc_global_dn(1,1),Nsite,zero,BCS_sc(1,1),Nsite)

!-------------------------------------------------
do i=1,Dtot,1
   if(decMethod .EQ. 2)then
      h0=0
      ph=0
      phi_sc=0

      h0(1:Nsite,1:Nsite)=Hzero(1:Nsite,1:Nsite)
      call input_phiT(Nsite,Nspin(1),h0,ph(1,1))
      phi_sc(1:Nsite,1:Nspin(1),i)=ph(1:Nsite,1:Nspin(1))

      h0(1:Nsite,1:Nsite)=Hzero((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))
      call input_phiT(Nsite,Nspin(2),h0,ph(1,1))
      phi_sc(1:Nsite,Nspin(1)+1:Ntot,i)=ph(1:Nsite,1:Nspin(2))
   else
      phi_sc=0
      phi_sc(1:Nsite,1:Nspin(1),i)=cicj_sc_global_up(1:Nsite,Nsite-i*Nspin(1)+1:Nsite-(i-1)*Nspin(1))
      phi_sc(1:Nsite,Nspin(1)+1:Ntot,i)=cicj_sc_global_dn(1:Nsite,Nsite-i*Nspin(2)+1:Nsite-(i-1)*Nspin(2)) 
      !coe_multi(i)=sqrt(sum(ev_up(Nsite-i*Nspin(1)+1:Nsite-(i-1)*Nspin(1)))+sum(ev_dn(Nsite-i*Nspin(2)+1:Nsite-(i-1)*Nspin(2))))
   endif
enddo

call norm_array(Dtot,coe_multi(1))


end subroutine cc_decomposition


subroutine BCS_det_sc_ovp(cicj_sc_global_up,cicj_sc_global_dn,ev_m,det_ev_m,alpha)
use param
use phi_param
use phiT_param
use lattice_param
use model_param
use sc_loop_param

use caldet_module

implicit none
complex(kind=8),intent(IN)::cicj_sc_global_up(Nsite,Nsite),cicj_sc_global_dn(Nsite,Nsite)
complex(kind=8),intent(INOUT)::alpha
complex(kind=8),intent(IN)::ev_m(Nsite,Nsite),det_ev_m(Nsite,Nsite)

complex(kind=8)::temp_BCS_sc(Nsite,Nsite),BCS_sc2(Nsite,Nsite)
complex(kind=8)::det_sc(2*Nsite,Ntot)

complex(kind=8)::BCS_ovp(Nspin(1),Nspin(1)),det_ovp(Ntot,Ntot)
complex(kind=8)::M_one

integer::i,j,flag,counter

!---------------
!get BCS state
!---------------
temp_BCS_sc=zero
BCS_sc2=zero
call zgemm('N','N',Nsite,Nsite,Nsite,one,cicj_sc_global_up(1,1),Nsite,ev_m(1,1),Nsite,zero,temp_BCS_sc(1,1),Nsite)
call zgemm('N','T',Nsite,Nsite,Nsite,one,temp_BCS_sc(1,1),Nsite,cicj_sc_global_dn(1,1),Nsite,zero,BCS_sc2(1,1),Nsite)

!-------------
!get det_sc
!-------------
det_sc=zero
counter=0
do i=1,Nsite
   if(det_ev_m(i,i) .EQ. 1)then
      counter=counter+1
      det_sc(1:Nsite,counter)=cicj_sc_global_up(1:Nsite,i)
      det_sc(Nsite+1:2*Nsite,Nspin(1)+counter)=cicj_sc_global_dn(1:Nsite,i)
   endif
enddo
if(counter .NE. Nspin(1))then
   write(*,*)'some problems in det_ev_m'
   call mystop
endif

!--------------------------------------------
!calculate alpha=<BCS|det_sc>/<det_sc|det_sc>
!--------------------------------------------
alpha=zero
BCS_ovp=zero
call bcs_over_lap_dc(BCS_sc2,det_sc,BCS_ovp)
call caldet(Nspin(1),BCS_ovp(1:Nspin(1),1:Nspin(1)),alpha)

M_one=zero
det_ovp=zero
call over_lap_dc(det_sc,det_sc,det_ovp)
call caldet(Ntot,det_ovp(1:Ntot,1:Ntot),M_one)

alpha=alpha/M_one
write(*,*)'M_one:',M_one

end subroutine BCS_det_sc_ovp


subroutine read_cicj_sc_global()
use sc_loop_param
use lattice_param

implicit none
integer::i,j

open(unit=10,file='cicj_sc_global_input.inputdat',status='old')
do i=1,2*Nsite,1
  do j=1,2*Nsite,1
        read(10,*) cicj_sc_global(i,j)
  end do
end do
close(10)

end subroutine read_cicj_sc_global


subroutine GF_Modifier()
use sc_loop_param
use lattice_param
use model_param

implicit none
complex(kind=8),allocatable::rho(:)
complex(kind=8),allocatable::rho_counter(:)
complex(kind=8)::rho01,rho02
real(kind=8)::distance
integer::sitex,sitey,Nsitex,Nsitey
integer::i,j

!-------------
!collect
!-------------
allocate(rho(5))
allocate(rho_counter(5))

rho=0
rho_counter=0
Nsitex=Nl(1)
Nsitey=Nl(2)

do sitex=1,Nsitex
do sitey=1,Nsitey
   do i=1,Nsitex
      do j=1,Nsitey
         distance=(sitex-i)**2+(sitey-j)**2
         !spin up
         if(distance .EQ. 1)rho(1)=rho(1)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)   
         if(distance .EQ. 9)rho(1)=rho(1)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)
         if(distance .EQ. 2)rho(2)=rho(2)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)
         if(distance .EQ. 10)rho(2)=rho(2)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)
         if(distance .EQ. 18)rho(2)=rho(2)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)
         if(distance .EQ. 4)rho(3)=rho(3)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)
         if(distance .EQ. 5)rho(4)=rho(4)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)
         if(distance .EQ. 13)rho(4)=rho(4)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)
         if(distance .EQ. 8)rho(5)=rho(5)+cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)
         !spin dn
         if(distance .EQ. 1)rho(1)=rho(1)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex) 
         if(distance .EQ. 9)rho(1)=rho(1)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex) 
         if(distance .EQ. 2)rho(2)=rho(2)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex) 
         if(distance .EQ. 10)rho(2)=rho(2)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)
         if(distance .EQ. 18)rho(2)=rho(2)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex) 
         if(distance .EQ. 4)rho(3)=rho(3)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)
         if(distance .EQ. 5)rho(4)=rho(4)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)
         if(distance .EQ. 13)rho(4)=rho(4)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)
         if(distance .EQ. 8)rho(5)=rho(5)+cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)
      enddo
   enddo
enddo
enddo

rho_counter(1)=8*Nsite
rho_counter(2)=8*Nsite
rho_counter(3)=4*Nsite
rho_counter(4)=8*Nsite
rho_counter(5)=2*Nsite

rho=rho/rho_counter
!rho01=(real(Nspin(1))/real(Nsite))*1.1
!rho02=(real(Nspin(1))/real(Nsite))*0.9
rho01=0.625
rho02=0.25

write(*,*)'hhhhhhhhhhhhhhh',rho
!------------------
!set 
!------------------
do sitex=1,Nsitex
do sitey=1,Nsitey
   if( mod((sitex+sitey),2) .EQ. 1 )cicj_sc_global(sitex+(sitey-1)*Nsitex,sitex+(sitey-1)*Nsitex)=rho02
   if( mod((sitex+sitey),2) .EQ. 0 )cicj_sc_global(sitex+(sitey-1)*Nsitex,sitex+(sitey-1)*Nsitex)=rho01
   if( mod((sitex+sitey),2) .EQ. 1 )cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+sitex+(sitey-1)*Nsitex)=rho01
   if( mod((sitex+sitey),2) .EQ. 0 )cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+sitex+(sitey-1)*Nsitex)=rho02
   do i=1,Nsitex
      do j=1,Nsitey
         distance=(sitex-i)**2+(sitey-j)**2
         !spin up
         if(distance .EQ. 1)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(1)   
         if(distance .EQ. 9)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(1)  
         if(distance .EQ. 2)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(2) 
         if(distance .EQ. 10)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(2)
         if(distance .EQ. 18)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(2) 
         if(distance .EQ. 4)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(3) 
         if(distance .EQ. 5)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(4)
         if(distance .EQ. 13)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(4) 
         if(distance .EQ. 8)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho(5)
         !spin dn
         if(distance .EQ. 1)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(1)   
         if(distance .EQ. 9)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(1)  
         if(distance .EQ. 2)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(2) 
         if(distance .EQ. 10)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(2)
         if(distance .EQ. 18)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(2) 
         if(distance .EQ. 4)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(3) 
         if(distance .EQ. 5)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(4)
         if(distance .EQ. 13)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(4) 
         if(distance .EQ. 8)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho(5)
      enddo
   enddo
enddo
enddo



end subroutine GF_Modifier









