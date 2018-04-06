
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
       !Set from the FS wave function.
       write(*,*)'get phiT from FS'
       call set_FS_phiT()
    elseif(sc_loop_flag .EQ. 1)then
       write(*,*)'get phiT from cc'
       call get_cc_phiT()
    endif
  endif   
  
  if(I_wavefun.eq.2)then
    !call read_F_matrix()
    !if(Nzeta.gt.0)then
    !  call read_D_unpaired()
    !endif
    !-----------------------------------
    !-----------------------------------
    if(sc_loop_flag .EQ. 0) then
       !Set from the FS wave function.
       write(*,*)'get phiT from FS'
       call set_FS_phiT()
    elseif(sc_loop_flag .EQ. 1)then
       write(*,*)'get phiT from cc'
       call get_cc_phiT()
!-------------------------------------------------
!trans the notation of BCS_sc to Fpairing
!-------------------------------------------------
FPairing=BCS_sc
!-------------------------------------------------
    endif
    !-----------------------------------
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

write(*,*)h0

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

integer::i
real(kind=8):: a
complex(kind=8)::alpha
complex(kind=8)::b,delta
!var-------------
integer::var_num
real(kind=8)::var_gap
complex(kind=8)::var_step
!----------------

if(sc_ite_flag==1) then
   a=0.75
elseif(sc_ite_flag==0) then
   a=0
endif

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
!check cc_decomposition
!----------------------------------------
!--------------------------------
 call openUnit('EV_UP.dat',81,'R')
 call openUnit('EV_DN.dat',82,'R')
 do i=1,Nsite,1
    write(81,*)'This is ev_up from cicj_sc_global_up, the order is',i,'the value is',ev_up(i)
    write(82,*)'This is ev_dn from cicj_sc_global_dn, the order is',i,'the value is',ev_dn(i)
 enddo

 write(81,*)'coe_muti:',coe_multi,'sum_ev',sum_ev
 write(82,*)'coe_muti:',coe_multi,'sum-ev',sum_ev

 temp_up=0
 temp_dn=0
 do i=1,Dtot,1
 temp_up=temp_up+phi_sc(1:Nsite,1:Nspin(1),i)*transpose(conjg(phi_sc(1:Nsite,1:Nspin(1),i)))*coe_multi(i)**2
 temp_dn=temp_dn+phi_sc(1:Nsite,Nspin(1)+1:Ntot,i)*transpose(conjg(phi_sc(1:Nsite,Nspin(1)+1:Ntot,i)))*coe_multi(i)**2
 enddo
 temp_up=abs(cicj_sc_global_up-temp_up)
 temp_dn=abs(cicj_sc_global_dn-temp_dn)
 sum_up=sum(temp_up(1:Nsite,1:Nsite))
 sum_dn=sum(temp_dn(1:Nsite,1:Nsite))
 write(81,*)'cc_decomp_error:',sum_up
 write(82,*)'cc_decomp_error:',sum_dn

 close(81)
 close(82)
!-----------------------------------------
!-----------------------------------------
!get BCS_SC
!--------
!get ev_m
ev_m=0
sum_ev=0
do i=Nsite-Nspin(1)+1,Nsite
    ev_m(i,i)=1
enddo
!--------------------
!for 1 x
!--------------------
!b=ev_up(Nsite-Nspin(1))
!do i=Nsite-Nspin(1),Nsite-Nspin(1)
!    ev_m(i,i)=sqrt(b/(Nspin(1)-Nspin(1)*b))
!enddo
!--------------------
!for 2 x
!--------------------
b=(ev_up(Nsite-Nspin(1)-1)+ev_dn(Nsite-Nspin(1)-1))/2
delta=(Nspin(1)*(1-2*b))**2+2*b*(Nspin(1)-1)*(Nspin(1)-Nspin(1)*b)
do i=Nsite-Nspin(1)-1,Nsite-Nspin(1)
    ev_m(i,i)=sqrt( ((2*b-1)*Nspin(1)+sqrt(delta))/((Nspin(1)-1)*(Nspin(1)-Nspin(1)*b)) )
enddo
write(*,*)'ev_m(Nsite-Nspin(1),Nsite-Nspin(1)):',ev_m(Nsite-Nspin(1),Nsite-Nspin(1))
!----------------------------------------
!get det_ev_m and check the decomposition
!----------------------------------------
!det_ev_m=0
!do i=Nsite-Nspin(1)+1,Nsite
!   det_ev_m(i,i)=1
!enddo
!do i=1,Nsite-Nspin(1),1
!   !det_ev_m(i,i)=1
!enddo
!------------------------------------------
!call BCS_det_sc_ovp(cicj_sc_global_up,cicj_sc_global_dn,ev_m,det_ev_m,alpha)
!write(*,*)'alpha=<BCS|det_sc>/<det_sc|det_sc>:',alpha
!---------------------------------
!ev_m=ev_m/(sum_ev/real(Nspin(1)))
!BCS_sc=matmul(matmul(cicj_sc_global_up,sqrt(ev_m)),transpose(cicj_sc_global_dn))
temp_BCS_sc=zero
BCS_sc=zero
call zgemm('N','N',Nsite,Nsite,Nsite,one,cicj_sc_global_up(1,1),Nsite,ev_m(1,1),Nsite,zero,temp_BCS_sc(1,1),Nsite)
call zgemm('N','T',Nsite,Nsite,Nsite,one,temp_BCS_sc(1,1),Nsite,cicj_sc_global_dn(1,1),Nsite,zero,BCS_sc(1,1),Nsite)

!-------------------------------------------------
!check the phi_initial
!-------------------------------------------------
!cicj_sc_global_up(1:Nsite,1:Nsite)=cicj_sc_global_up1(1:Nsite,1:Nsite)
!cicj_sc_global_dn(1:Nsite,1:Nsite)=cicj_sc_global_dn1(1:Nsite,1:Nsite)
!call eigen(cicj_sc_global_up,Nsite,ev_up)
!call eigen(cicj_sc_global_dn,Nsite,ev_dn)
!-------------------------------------------------
!-------------------------------------------------
do i=1,Dtot,1
   phi_sc(1:Nsite,1:Nspin(1),i)=cicj_sc_global_up(1:Nsite,Nsite-i*Nspin(1)+1:Nsite-(i-1)*Nspin(1))
   phi_sc(1:Nsite,Nspin(1)+1:Ntot,i)=cicj_sc_global_dn(1:Nsite,Nsite-i*Nspin(2)+1:Nsite-(i-1)*Nspin(2)) 
   !coe_multi(i)=sqrt(sum(ev_up(Nsite-i*Nspin(1)+1:Nsite-(i-1)*Nspin(1)))+sum(ev_dn(Nsite-i*Nspin(2)+1:Nsite-(i-1)*Nspin(2))))
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
