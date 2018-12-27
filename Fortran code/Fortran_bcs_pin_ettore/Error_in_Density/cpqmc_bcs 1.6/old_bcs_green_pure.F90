subroutine bcs_green_pure(i,phR,Green)
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param
use method_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer,intent(IN)::i
complex(kind=8),intent(IN)::phR(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::Green(2*Nsite,2*Nsite)
integer::ib,jb,istep,Nsteps,l,j
complex(kind=8)::ph(2*Nsite,Ntot)
complex(kind=8)::Bforward_Up(Nsite,Nsite),Bforward_Dn(Nsite,Nsite)
complex(kind=8)::Bbackward_Up(Nsite,Nsite),Bbackward_Dn(Nsite,Nsite)
complex(kind=8)::DotB_up(Nsite,Nsite),DotB_dn(Nsite,Nsite),f_up(Nsite),f_dn(Nsite)
complex(kind=8)::AA(Nspin(1),Nspin(1))
complex(kind=8)::FstarPhidn(Nsite,Nspin(1)),NewFPhd(Nsite,Nspin(1)),NewPhu(Nsite,Nspin(1))
complex(kind=8)::Am1Phut(Nspin(1),Nsite)
complex(kind=8)::imp
integer::i_deb,j_deb

ph(:,:)=phR(:,:)
Bforward_Up=zero;Bforward_Dn=zero;Bbackward_Up=zero;Bbackward_Dn=zero
DotB_up=zero;DotB_dn=zero
f_up=one;f_dn=one

do ib=1,Nsite
  Bforward_Up(ib,ib)=one
  Bforward_Dn(ib,ib)=one
  Bbackward_Up(ib,ib)=one
  Bbackward_Dn(ib,ib)=one
  DotB_up(ib,ib)=one
  DotB_dn(ib,ib)=one
enddo

if(back_pro)then
  Nsteps=Nstps_fwd
else
  Nsteps=0
endif

!DEBUG
write(*,*)
write(*,*)'Pure subroutine '
write(*,*)
write(*,*)'Nsteps = ',Nsteps
write(*,*)

do istep=1,Nsteps
  call move(i,istep,ph,Bforward_Up,Bforward_Dn,Bbackward_Up,Bbackward_Dn)
  if(mod(istep,StepforGram).eq.0)then
!DEBUG
    write(*,*)
    write(*,*)'Call stabilization  '
    write(*,*)
    call stabilize(ph,Bforward_Up,Bforward_Dn,Bbackward_Up,Bbackward_Dn,DotB_up,DotB_dn,f_up,f_dn)
  endif
enddo

call bcs_over_lap_dc(Fpairing,ph,AA(1:Nspin(1),1:Nspin(1)))
call inverse_d(AA(1:Nspin(1),1:Nspin(1)),Nspin(1),imp)
call zgemm('N','N',Nsite,   Nspin(1),Nsite,   one,conjg(Fpairing(:,:)),Nsite      &
     &     ,ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,FstarPhidn,Nsite)
call zgemm('T','N',Nsite,Nspin(1),Nsite,   one,conjg(Bforward_Up(:,:)),Nsite      &
     &     ,FstarPhidn,Nsite,zero,NewFPhd,Nsite)
call zgemm('T','N',Nsite,Nspin(1),Nsite,   one,Bbackward_Up,Nsite                 &
     &     ,ph(1:Nsite,1:Nspin(1)),Nsite,zero,NewPhu,Nsite)
call zgemm('N','T',Nspin(1),Nsite,Nspin(1),one,AA,Nspin(1)                        &
     &      ,NewPhu,Nsite,zero,Am1Phut,Nspin(1))
call zgemm('N','N',Nsite   ,Nsite,   Nspin(1),one,NewFPhd,Nsite                   &
     &      ,Am1Phut,Nspin(1),zero,Green(1:Nsite,1:Nsite),Nsite)


do l=1,Nsite
  do j=1,Nsite
    Green(j,l)=Green(j,l)-DotB_up(l,j)
  enddo
  Green(l,l)=Green(l,l)+one
enddo


!DEBUG
write(*,*)
write(*,*)'Pure subroutine '
write(*,*)
write(*,*)'Nsteps = ',Nsteps
write(*,*)
   write(*,*)'walker ',i
   write(*,*)
   write(*,*)
   write(*,*)'Fpairing '
   do i_deb=1,Nsite
     write(*,*)(Fpairing(i_deb,j_deb),j_deb=1,Nsite)
   enddo
   write(*,*)
   write(*,*)'phi up '
   do i_deb=1,Nsite
     write(*,*)(ph(i_deb,j_deb),j_deb=1,Nspin(1))
   enddo 
   write(*,*)
   write(*,*)'phi dn '
   do i_deb=Nsite+1,2*Nsite
     write(*,*)(ph(i_deb,j_deb),j_deb=Nspin(1)+1,Ntot)
   enddo
   write(*,*)
write(*,*)
write(*,*)'<BCS|ph> = ',imp
write(*,*)'Bforward_Up '
   do i_deb=1,Nsite
     write(*,*)(Bforward_Up(i_deb,j_deb),j_deb=1,Nsite)
   enddo
write(*,*)
write(*,*)'Bbackward_Up '
   do i_deb=1,Nsite
     write(*,*)(Bbackward_Up(i_deb,j_deb),j_deb=1,Nsite)
   enddo
write(*,*)
write(*,*)'FstarPhidn '
   do i_deb=1,Nsite
     write(*,*)(FstarPhidn(i_deb,j_deb),j_deb=1,Nspin(1))
   enddo
write(*,*)
write(*,*)'NewFPhd '
   do i_deb=1,Nsite
     write(*,*)(NewFPhd(i_deb,j_deb),j_deb=1,Nspin(1))
   enddo
write(*,*)

end subroutine bcs_green_pure



subroutine move(i,i_b,ph,Bf_up,Bf_dn,Bb_up,Bb_dn)
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param
use method_param
use project_param
implicit none
integer,intent(IN)::i,i_b
complex(kind=8),intent(INOUT)::ph(2*Nsite,Ntot)
complex(kind=8),intent(INOUT)::Bf_up(Nsite,Nsite),Bf_dn(Nsite,Nsite),Bb_up(Nsite,Nsite),Bb_dn(Nsite,Nsite)
complex(kind=8)::temp(2*Nsite,Ntot),Btmp(Nsite,Nsite)
integer::j,l
complex(kind=8)::explr_up,explr_dn
real(kind=8)::x
integer::aux

!T/2                                  
call copy_wf_dc(ph(1,1),temp(1,1))
    !complexity L^2 * N
call k_to_ph_dc(exp_halfK,temp(1,1),ph(1,1))   
    !complexity L^3
call zcopy(Nsite*Nsite,Bf_up,1,Btmp,1)
call zgemm('N','N',Nsite,Nsite,Nsite,one,exp_halfK(1,1),2*Nsite,Btmp,Nsite,zero,Bf_up,Nsite)
call zcopy(Nsite*Nsite,Bf_dn,1,Btmp,1)
call zgemm('N','N',Nsite,Nsite,Nsite,one,exp_halfK(Nsite+1,Nsite+1),2*Nsite,Btmp,Nsite,zero,Bf_dn,Nsite)
     !kinetic energy is hermitian
call zcopy(Nsite*Nsite,Bb_up,1,Btmp,1)
call zgemm('N','N',Nsite,Nsite,Nsite,one,exp_mhalfK(1,1),2*Nsite,Btmp,Nsite,zero,Bb_up,Nsite)
call zcopy(Nsite*Nsite,Bb_dn,1,Btmp,1)
call zgemm('N','N',Nsite,Nsite,Nsite,one,exp_mhalfK(Nsite+1,Nsite+1),2*Nsite,Btmp,Nsite,zero,Bb_dn,Nsite)


!V
do j=1,Nsite,1
  x=back_store(j,i_b,i)
  aux=NINT(x)
  explr_up=expln_up(aux)
  explr_dn=expln_dn(aux)
  do l=1,Nspin(1),1          !complexity L*N
    ph(j,l)=ph(j,l)*(explr_up+one)
  end do
  do l=Nspin(1)+1,Ntot,1
    ph(j+Nsite,l)=ph(j+Nsite,l)*(explr_dn+one)
  enddo
  do l=1,Nsite
    Bf_up(j,l)=Bf_up(j,l)*(explr_up+one)
    Bf_dn(j,l)=Bf_dn(j,l)*(explr_dn+one)
    Bb_up(j,l)=Bb_up(j,l)*conjg(one/(explr_up+one))
    Bb_dn(j,l)=Bb_dn(j,l)*conjg(one/(explr_dn+one))
  enddo
enddo

!T/2
call copy_wf_dc(ph(1,1),temp(1,1))
call k_to_ph_dc(exp_halfK,temp(1,1),ph(1,1))  
call zcopy(Nsite*Nsite,Bf_up,1,Btmp,1)
call zgemm('N','N',Nsite,Nsite,Nsite,one,exp_halfK(1,1),2*Nsite,Btmp,Nsite,zero,Bf_up,Nsite)
call zcopy(Nsite*Nsite,Bf_dn,1,Btmp,1)
call zgemm('N','N',Nsite,Nsite,Nsite,one,exp_halfK(Nsite+1,Nsite+1),2*Nsite,Btmp,Nsite,zero,Bf_dn,Nsite)
call zcopy(Nsite*Nsite,Bb_up,1,Btmp,1)
call zgemm('N','N',Nsite,Nsite,Nsite,one,exp_mhalfK(1,1),2*Nsite,Btmp,Nsite,zero,Bb_up,Nsite)
call zcopy(Nsite*Nsite,Bb_dn,1,Btmp,1)
call zgemm('N','N',Nsite,Nsite,Nsite,one,exp_mhalfK(Nsite+1,Nsite+1),2*Nsite,Btmp,Nsite,zero,Bb_dn,Nsite)


end subroutine move


subroutine stabilize(ph,Bf_up,Bf_dn,Bb_up,Bb_dn,D_up,D_dn,f_up,f_dn)
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param
use method_param
implicit none
complex(kind=8),intent(INOUT)::ph(2*Nsite,Ntot)
complex(kind=8),intent(INOUT)::Bf_up(Nsite,Nsite),Bf_dn(Nsite,Nsite),Bb_up(Nsite,Nsite),Bb_dn(Nsite,Nsite)
complex(kind=8),intent(INOUT)::D_up(Nsite,Nsite),D_dn(Nsite,Nsite),f_up(Nsite),f_dn(Nsite)
complex(kind=8)::Rmat(Ntot,Ntot)
complex(kind=8)::mup(Nsite),ui(Nsite),mup_tilde(Nsite)
real(kind=8)::norm(Nsite)
real(kind=8)::anm1,anm2
integer::j,l,ip
complex(kind=8)::proj(Ntot,Nsite),projb(Ntot,Nsite)
complex(kind=8),external::zdotc
complex(kind=8)::dummy
integer::i_deb,j_deb


!DEBUG
write(*,*)'Call Gram Schmidt'
write(*,*)

!make the orbitals in the wave function orthonormal
call modGS(ph(1:Nsite,1:Nspin(1)),Nsite,Nspin(1),anm1,Rmat(1:Nspin(1),1:Nspin(1)))
call modGS(ph((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot),Nsite,Nspin(2),anm2,Rmat(1:Nspin(2),1:Nspin(2)))

!DEBUG
write(*,*)'Exit from Gram Schmidt'
write(*,*)
write(*,*)
   write(*,*)'phi up '
   do i_deb=1,Nsite
     write(*,*)(ph(i_deb,j_deb),j_deb=1,Nspin(1))
   enddo
   write(*,*)
   write(*,*)'phi dn '
   do i_deb=Nsite+1,2*Nsite
     write(*,*)(ph(i_deb,j_deb),j_deb=Nspin(1)+1,Ntot)
   enddo
   write(*,*)
write(*,*)'Bforward_Up '
   do i_deb=1,Nsite
     write(*,*)(Bf_Up(i_deb,j_deb),j_deb=1,Nsite)
   enddo
write(*,*)
write(*,*)'Bbackward_Up '
   do i_deb=1,Nsite
     write(*,*)(Bb_Up(i_deb,j_deb),j_deb=1,Nsite)
   enddo
write(*,*)

!Remember we are buildind  
!  c_{B^{-1,T} alpha} c+_{B beta} phi+_1 .... phi+_N |0>
! We have applied GS building an orthonormal u+_1 ..... u+_N
! We want to change
! |B beta> ---> |beta'> = B beta> - sum_i=1^N  < u_i | B beta > |u _i >
! normalized
! this does not change the ket

do l=1,Nsite
  mup(:)=Bf_up(:,l)

!DEBUG
!write(*,*)
!write(*,*)
!write(*,*)'l = ',l
!write(*,*)'|mu l> = '
!do i_deb=1,Nsite
!  write(*,*)mup(i_deb)
!enddo
!write(*,*)
!write(*,*)
  

  do ip=1,Nspin(1)
    ui(1:Nsite)=ph(1:Nsite,ip)
    proj(ip,l)=zdotc(Nsite,ui,1,mup,1)


!DEBUG
!write(*,*)
!write(*,*)'ip = ',ip
!write(*,*)'| u ip > = '
!do i_deb=1,Nsite
!  write(*,*)ui(i_deb)
!enddo
!write(*,*)
!write(*,*)' < u ip | mu l > = ',proj(ip,l)
!write(*,*)
!write(*,*)


  enddo
enddo

!DEBUG
write(*,*)'Computed projections < u_ip | B l >, l=beta'
write(*,*)

! proj(ip,l) = < u_ip | B l >, l=beta

do l=1,Nsite
  do j=1,Nsite
    mup_tilde(j)=Bf_up(j,l)
    do ip=1,Nspin(1)
      mup_tilde(j)=mup_tilde(j)-proj(ip,l)*ph(j,ip)
    enddo
  enddo
  norm(l)=1.d0 !dsqrt(dble(zdotc(Nsite,mup_tilde,1,mup_tilde,1)))
  f_up(l)=f_up(l)*norm(l)
  if(norm(l).gt.1.d-10)then
    do j=1,Nsite
      Bf_up(j,l)=mup_tilde(j)/norm(l)
    enddo
  else
    do j=1,Nsite
      Bf_up(j,l)=zero
    enddo
  endif


enddo

!DEBUG
write(*,*)'Updated Bf_up'
write(*,*)
write(*,*)'Bforward_Up '
   do i_deb=1,Nsite
     write(*,*)(Bf_Up(i_deb,j_deb),j_deb=1,Nsite)
   enddo
write(*,*)
write(*,*)

!replaced  |B beta> ---> |beta'> = B beta> - sum_i=1^N  < u_i | B beta > |u _i >
!we now want to update the scalar product matrix
! < B^{-1,T} alpha | beta'> = < B^{-1,T} alpha | B beta> 
! - sum_i=1^N  < u_i | B beta >  < B^{-1,T} alpha | u _i >


do l=1,Nsite
  mup(:)=Bb_up(:,l)
  do ip=1,Nspin(1)
    ui(1:Nsite)=ph(1:Nsite,ip)
    projb(ip,l)=zdotc(Nsite,ui,1,mup,1)
  enddo
enddo

!DEBUG
write(*,*)'Computed projections < u_ip | B^{-1,T} l >, l=beta'
write(*,*)
do i_deb=1,Nsite
  write(*,*)(projb(j_deb,i_deb),j_deb=1,Nspin(1))
enddo


! projb(ip,l) = < u_ip | B^{-1,T} l >

!DEBUG
write(*,*)' D_up'
write(*,*)
write(*,*)
write(*,*)'D_up old '
   do i_deb=1,Nsite
     write(*,*)(D_up(i_deb,j_deb),j_deb=1,Nsite)
   enddo
write(*,*)
write(*,*)

do l=1,Nsite

!DEBUG
  write(*,*)'norm ',l,norm(l)

  do j=1,Nsite
    dummy=D_up(j,l)
    if(norm(l).gt.1.d-10)then
      D_up(j,l)=dummy
      do ip=1,Nspin(1)
        D_up(j,l)=D_up(j,l)-conjg(projb(ip,j))*proj(ip,l)
      enddo
      D_up(j,l)=D_up(j,l)/norm(l)
    else
      D_up(j,l)=zero
    endif
  enddo
enddo

!DEBUG
write(*,*)'Updated D_up'
write(*,*)
write(*,*)
write(*,*)'D_up new '
   do i_deb=1,Nsite
     write(*,*)(D_up(i_deb,j_deb),j_deb=1,Nsite)
   enddo
write(*,*)
write(*,*)
!TUTTO OK FINO QUI

!  D_up(j,l) = < B^{-1,T} alpha | beta' >

!Now we want to replace
! | B^{-1,T} alpha > --> | alpha' > = sum_i=1^N < u_ip | B^{-1,T} alpha > u_i
! + < beta' | B^{-1,T} alpha > | beta '>

do l=1,Nsite
  do j=1,Nsite
   ! mup_tilde(j)=Bf_up(j,l)*conjg(D_up(l,j))
    mup_tilde(j)=zero !Bb_up(j,l)
    do ip=1,Nspin(1)
      mup_tilde(j)=mup_tilde(j)+projb(ip,l)*ph(j,ip)
    enddo
  enddo
  do j=1,Nsite
    Bb_up(j,l)=mup_tilde(j)
  enddo
enddo

!DEBUG
write(*,*)'Updated Bb_up'
write(*,*)
write(*,*)'Bbackword_Up '
   do i_deb=1,Nsite
     write(*,*)(Bb_Up(i_deb,j_deb),j_deb=1,Nsite)
   enddo
write(*,*)
write(*,*)

!last point update D = <alpha' | beta'> but this already coincides with
!D_up(j,l)


!---- SPIN DOWN SAME PROCEDURE

do l=1,Nsite
  mup(:)=Bf_dn(:,l)
  do ip=Nspin(1)+1,Ntot
    ui(1:Nsite)=ph(Nsite+1:2*Nsite,ip)
    proj(ip,l)=zdotc(Nsite,ui,1,mup,1)
  enddo
enddo

!DEBUG
write(*,*)'DOWN Computed projections < u_ip | B l >, l=beta'
write(*,*)

do l=1,Nsite
  do j=1,Nsite
    mup_tilde(j)=mup(j)
    do ip=Nspin(1)+1,Ntot
      mup_tilde(j)=mup_tilde(j)-proj(ip,l)*ph(Nsite+j,ip)
    enddo
  enddo
  norm(l)=dsqrt(dble(zdotc(Nsite,mup_tilde,1,mup_tilde,1)))
  f_dn(l)=f_dn(l)*norm(l)
  if(norm(l).gt.1.d-10)then
    do j=1,Nsite
      Bf_dn(j,l)=mup_tilde(j)/norm(l)
    enddo
  else
    do j=1,Nsite
      Bf_dn(j,l)=zero
    enddo
  endif
enddo

!DEBUG
write(*,*)'Updated Bf_dn'
write(*,*)

do l=1,Nsite
  mup(:)=Bb_dn(:,l)
  do ip=Nspin(1)+1,Ntot
    ui(1:Nsite)=ph(Nsite+1:2*Nsite,ip)
    projb(ip,l)=zdotc(Nsite,ui,1,mup,1)
  enddo
enddo

do l=1,Nsite
  do j=1,Nsite
    if(norm(l).gt.1.d-10)then
      do ip=1,Nspin(1)
        D_dn(j,l)=D_dn(j,l)-conjg(projb(ip,j))*proj(ip,l)
      enddo
    else
      D_dn(j,l)=zero
    endif
  enddo
enddo

!DEBUG
write(*,*)'DOWN Computed projections < u_ip | B^{-1,T} l >, l=beta'
write(*,*)


do l=1,Nsite
  mup(:)=Bf_dn(:,l)
  do j=1,Nsite
    mup_tilde(j)=mup(j)*conjg(D_dn(j,l))
    do ip=Nspin(1)+1,Ntot
      mup_tilde(j)=mup_tilde(j)+projb(ip,l)*ph(j+Nsite,ip)
    enddo
  enddo
  do j=1,Nsite
    Bb_dn(j,l)=mup_tilde(j)
  enddo
enddo


!DEBUG
write(*,*)'Updated Bb_dn'
write(*,*)


end subroutine stabilize
