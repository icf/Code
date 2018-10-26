subroutine bcs_green_pure(i,phR,Green,when)
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
logical, intent(IN)::when
integer::ib,jb,istep,Nsteps,l,j
complex(kind=8)::ph(2*Nsite,Ntot)
complex(kind=8)::Bforward_Up(Nsite,Nsite),Bforward_Dn(Nsite,Nsite)
complex(kind=8)::Bbackward_Up(Nsite,Nsite),Bbackward_Dn(Nsite,Nsite)
complex(kind=8)::DotB_up(Nsite,Nsite),DotB_dn(Nsite,Nsite)

!---------------------------------------
!parameters for coupled form
!---------------------------------------
complex(kind=8)::Bforward(2*Nsite,2*Nsite)
complex(kind=8)::Bbackward(2*Nsite,2*Nsite)
complex(kind=8)::DotB(2*Nsite,2*Nsite)
complex(kind=8)::dummy1(2*Nsite,Ntot),dummy2(Ntot,2*Nsite),dummy3(2*Nsite,2*Nsite),dummy4(2*Nsite,2*Nsite)
real(kind=8)::anm

complex(kind=8)::AA(Nspin(1),Nspin(1)),AAsmall(Nspin(2),Nspin(1))
complex(kind=8)::cAA(Ntot,Ntot),cAA_temp(Ntot,Ntot)
complex(kind=8)::phiu(2*Nsite,Nspin(1)),phid(2*Nsite,Nspin(2))
complex(kind=8)::TrueFstarPhidn(2*Nsite,Nspin(2))
complex(kind=8)::FstarPhidn(2*Nsite,Nspin(1)),FDPhiup(2*Nsite,Nspin(1))
complex(kind=8)::NewPhd(2*Nsite,Nspin(2)),NewPhu(2*Nsite,Nspin(1))
complex(kind=8)::NewFDPhu(2*Nsite,Nspin(1)),NewFPhd(2*Nsite,Nspin(1))
complex(kind=8)::Am1Phut(Nspin(1),2*Nsite),Am1tPhdt(Nspin(1),2*Nsite)
complex(kind=8)::Anomalous(2*Nsite,2*Nsite),FBd(2*Nsite,2*Nsite),NewF(2*Nsite,2*Nsite)
complex(kind=8)::Rmat(Ntot,Ntot)
real(kind=8)::anm1,anm2,sommau,sommad
complex(kind=8)::imp
integer::i_deb,j_deb
logical::i_debug

i_debug=.false.

ph(:,:)=phR(:,:)

!----------------------------------------------------
!for dtype = 'd'
!----------------------------------------------------
if(dtype .EQ. 'd')then

Bforward_Up=zero;Bforward_Dn=zero;Bbackward_Up=zero;Bbackward_Dn=zero
DotB_up=zero;DotB_dn=zero

do ib=1,Nsite
  Bforward_Up(ib,ib)=one
  Bforward_Dn(ib,ib)=one
  Bbackward_Up(ib,ib)=one
  Bbackward_Dn(ib,ib)=one
!  DotB_up(ib,ib)=one
!  DotB_dn(ib,ib)=one
enddo

if(back_pro)then
  Nsteps=Nstps_fwd
else
  Nsteps=0
endif

if(when)then
  do istep=1,Nsteps
    call move(i,istep,ph,Bforward_Up,Bforward_Dn,Bbackward_Up,Bbackward_Dn)
    if(mod(istep,StepforGram).eq.0)then
      call modGS(ph(1:Nsite,1:Nspin(1)),Nsite,Nspin(1),anm1,Rmat(1:Nspin(1),1:Nspin(1)))
      call modGS(ph((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot),Nsite,Nspin(2),anm2,Rmat(1:Nspin(2),1:Nspin(2)))
      if(mod(istep,StepforPure).eq.0)then
        call stabilize(ph,Bforward_Up,Bforward_Dn,Bbackward_Up,Bbackward_Dn,DotB_up,DotB_dn)
      endif
    endif
  enddo
endif

if(Nzeta.eq.0)then
  call bcs_over_lap_d(Fpairing,ph,AA(1:Nspin(1),1:Nspin(1)))
else
  call unp_bcs_over_lap_dc(Fpairing,Dunpaired,ph,AA(1:Nspin(1),1:Nspin(1)))
endif
call inverse_d(AA(1:Nspin(1),1:Nspin(1)),Nspin(1),imp)
AAsmall(1:Nspin(2),1:Nspin(1))=AA(Nzeta+1:Nspin(1),1:Nspin(1))

phiu(1:Nsite,1:Nspin(1))=ph(1:Nsite,1:Nspin(1))
phid(1:Nsite,1:Nspin(2))=ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot)

!Green up up

call zgemm('N','N',Nsite,   Nspin(2),Nsite,   one,conjg(Fpairing(:,:)),Nsite      &
     &     ,phid,Nsite,zero,TrueFstarPhidn,Nsite)

if(Nzeta.eq.0)then
  FstarPhidn(:,:)=TrueFstarPhidn(:,:)
else
  FstarPhidn(1:Nsite,1:Nzeta)=conjg(Dunpaired(1:Nsite,1:Nzeta))
  FstarPhidn(1:Nsite,Nzeta+1:Nspin(1))=TrueFstarPhidn(1:Nsite,1:Nspin(2))
endif

call zgemm('T','N',Nsite,Nspin(1),Nsite,   one,conjg(Bbackward_Up(:,:)),Nsite     &
     &     ,phiu,Nsite,zero,NewPhu,Nsite)

call zgemm('T','N',Nsite,Nspin(1),Nsite,   one,Bforward_Up,Nsite                  &
     &     ,FstarPhidn,Nsite,zero,NewFPhd,Nsite)

call zgemm('N','T',Nspin(1),Nsite,Nspin(1),one,AA,Nspin(1)                        &
     &      ,NewPhu,Nsite,zero,Am1Phut,Nspin(1))

call zgemm('N','N',Nsite   ,Nsite,   Nspin(1),one,NewFPhd,Nsite                   &
     &      ,Am1Phut,Nspin(1),zero,Green(1:Nsite,1:Nsite),Nsite)

do l=1,Nsite
  do j=1,Nsite
    Green(j,l)=Green(j,l)+DotB_up(l,j)
  enddo
enddo

sommau=0.d0
do j=1,Nsite
!  write(*,*)'nup(j) ',j,Green(j,j)
  sommau=sommau+Green(j,j)
enddo
!write(*,*)'Nup = ',sommau

!Green dn dn


call zgemm('T','N',Nsite,   Nspin(1),Nsite,   one,conjg(Fpairing(:,:)),Nsite      &
     &     ,phiu,Nsite,zero,FDPhiup,Nsite)

call zgemm('T','N',Nsite,Nspin(1),Nsite,   one,Bforward_Dn,Nsite                  &
     &     ,FDPhiup,Nsite,zero,NewFDPhu,Nsite)

call zgemm('T','N',Nsite,Nspin(2),Nsite,   one,conjg(Bbackward_Dn(:,:)),Nsite     &
     &     ,phid,Nsite,zero,NewPhd,Nsite)

call zgemm('T','T',Nspin(1),Nsite,Nspin(2),one,AAsmall,Nspin(2)                   &
     &      ,NewPhd,Nsite,zero,Am1tPhdt,Nspin(1))

call zgemm('N','N',Nsite   ,Nsite,   Nspin(1),one,NewFDPhu,Nsite                   &
     &      ,Am1tPhdt,Nspin(1),zero,Green(Nsite+1:2*Nsite,Nsite+1:2*Nsite),Nsite)


do l=1,Nsite
  do j=1,Nsite
    Green(j+Nsite,l+Nsite)=Green(j+Nsite,l+Nsite)+DotB_dn(l,j)
  enddo
enddo

sommau=0.d0
do j=1,Nsite
!  write(*,*)'ndn(j) ',j,Green(j+Nsite,j+Nsite)
  sommau=sommau+Green(j+Nsite,j+Nsite)
enddo
!write(*,*)'Ndn = ',sommau
!stop'deb'

!Anomalous correlators

call zgemm('T','N',Nsite,Nspin(1),Nsite,   one,Bforward_Dn,Nsite                     &
     &     ,FDPhiup,Nsite,zero,NewFDPhu,Nsite)

call zgemm('N','T',Nspin(1),Nsite,Nspin(1),one,AA,Nspin(1)                            &
     &     ,NewFDPhu,Nsite,zero,Am1Phut,Nspin(1))

call zgemm('N','N',Nsite   ,Nsite,   Nspin(1),one,NewFPhd,Nsite                       &
     &      ,Am1Phut,Nspin(1),zero,Anomalous,Nsite)

call zgemm('N','N',Nsite   ,Nsite,Nsite,one,conjg(Fpairing(:,:)),Nsite                &
     &      ,Bforward_Dn,Nsite,zero,FBd,Nsite)

call zgemm('T','N',Nsite,Nsite,Nsite,   one,Bforward_Up,Nsite                         &
     &     ,FBd,Nsite,zero,NewF,Nsite)

do l=1,Nsite
  do j=1,Nsite
    Green(j,l+Nsite)=NewF(j,l)-Anomalous(j,l)     !<non projected BCS|c+_{j,up} c+_{l,dn}|W>
  enddo
enddo



call zgemm('N','N',Nsite   ,Nsite,   Nspin(1),one,NewPhu,Nsite                        &
     &      ,Am1tPhdt,Nspin(1),zero,Anomalous,Nsite)

do l=1,Nsite
  do j=1,Nsite
    Green(j+Nsite,l)=-Anomalous(j,l)       !   <non projected BCS|c_{j,up} c_{l,dn}|W>
  enddo
enddo

!----------------------------------------------------
!for dtype = 'c'
!----------------------------------------------------
elseif(dtype .EQ. 'c')then

   Bforward=zero;Bbackward=zero;DotB=zero

   do ib=1,2*Nsite
     Bforward(ib,ib)=one
     Bbackward(ib,ib)=one
   enddo

   if(back_pro)then
     Nsteps=Nstps_fwd
   else
     Nsteps=0
   endif

   if(when)then
     do istep=1,Nsteps
       call move_c(i,istep,ph,Bforward,Bbackward)
       if(mod(istep,StepforGram).eq.0)then
         call modGS(ph(1:2*Nsite,1:Ntot),2*Nsite,Ntot,anm,Rmat)
         if(mod(istep,StepforPure).eq.0)then
           call stabilize_c(ph,Bforward,Bbackward,DotB)
         endif
       endif
     enddo
   endif

   if(Nzeta.eq.0)then
     call bcs_over_lap_c(transpose(Fpairing),ph,cAA(1:Ntot,1:Ntot))
   else
     !call unp_bcs_over_lap_dc(Fpairing,Dunpaired,ph,AA(1:Nspin(1),1:Nspin(1)))
   endif

   cAA_temp=cAA
   call inverse_pf(cAA(1:Ntot,1:Ntot),Ntot,imp)

!Green=matmul( matmul( Transpose(Fpairing) ,ph ), matmul( cAA ,Transpose(ph) ) )
call zgemm('T','N',2*Nsite,Ntot,2*Nsite,one,conjg(Fpairing),2*Nsite,ph,2*Nsite,zero,dummy1,2*Nsite)
call zgemm('N','T',Ntot,2*Nsite,Ntot,one,cAA,Ntot,ph,2*Nsite,zero,dummy2,Ntot)
call zgemm('N','N',2*Nsite,2*Nsite,Ntot,one,dummy1,2*Nsite,dummy2,Ntot,zero,dummy3,2*Nsite)

!Green=matmul( matmul( transpose(Bforward) ,Green), Bbackward )
call zgemm('T','N',2*Nsite,2*Nsite,2*Nsite,one,Bforward,2*Nsite,dummy3,2*Nsite,zero,dummy4,2*Nsite)
call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,dummy4,2*Nsite,Bbackward,2*Nsite,zero,Green,2*Nsite)


   Green=Green+transpose(DotB)
   Green=transpose(Green)

   sommau=0.d0
   do j=1,Nsite
      !write(*,*)'nup(j) ',j,Green(j,j)
      sommau=sommau+Green(j,j)
   enddo
   !write(*,*)'Nup = ',sommau

   sommau=0.d0
   do j=Nsite+1,2*Nsite
     !write(*,*)'ndn(j) ',j,Green(j,j)
     sommau=sommau+Green(j,j)
   enddo
   !write(*,*)'Ndn = ',sommau


endif



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


subroutine move_c(i,i_b,ph,Bf,Bb)
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
complex(kind=8),intent(INOUT)::Bf(2*Nsite,2*Nsite),Bb(2*Nsite,2*Nsite)
complex(kind=8)::temp(2*Nsite,Ntot),Btmp(2*Nsite,2*Nsite)
integer::j,l
complex(kind=8)::explr_up,explr_dn
real(kind=8)::x
integer::aux

!T/2                                  
call copy_wf_dc(ph(1,1),temp(1,1))
    !complexity L^2 * N
call k_to_ph_dc(exp_halfK,temp(1,1),ph(1,1))   
    !complexity L^3
call zcopy(2*Nsite*2*Nsite,Bf,1,Btmp,1)
call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,exp_halfK,2*Nsite,Btmp,2*Nsite,zero,Bf,2*Nsite)
     !kinetic energy is hermitian
call zcopy(2*Nsite*2*Nsite,Bb,1,Btmp,1)
call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,exp_mhalfK,2*Nsite,Btmp,2*Nsite,zero,Bb,2*Nsite)
!V
do j=1,Nsite,1
  x=back_store(j,i_b,i)
  aux=NINT(x)
  explr_up=expln_up(aux)
  explr_dn=expln_dn(aux)
  do l=1,Ntot,1          !complexity L*N
    ph(j,l)=ph(j,l)*(explr_up+one)
  end do
  do l=1,Ntot,1
    ph(j+Nsite,l)=ph(j+Nsite,l)*(explr_dn+one)
  enddo
  do l=1,2*Nsite
    Bf(j,l)=Bf(j,l)*(explr_up+one)
    Bf(j+Nsite,l)=Bf(j+Nsite,l)*(explr_dn+one)
    Bb(j,l)=Bb(j,l)*conjg(one/(explr_up+one))
    Bb(j+Nsite,l)=Bb(j+Nsite,l)*conjg(one/(explr_dn+one))
  enddo
enddo

!T/2
call copy_wf_dc(ph(1,1),temp(1,1))
call k_to_ph_dc(exp_halfK,temp(1,1),ph(1,1))  
call zcopy(2*Nsite*2*Nsite,Bf,1,Btmp,1)
call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,exp_halfK,2*Nsite,Btmp,2*Nsite,zero,Bf,2*Nsite)
call zcopy(2*Nsite*2*Nsite,Bb,1,Btmp,1)
call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,exp_mhalfK,2*Nsite,Btmp,2*Nsite,zero,Bb,2*Nsite)

end subroutine move_c


subroutine stabilize(ph,Bf_up,Bf_dn,Bb_up,Bb_dn,D_up,D_dn)
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param
use method_param
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(INOUT)::Bf_up(Nsite,Nsite),Bf_dn(Nsite,Nsite),Bb_up(Nsite,Nsite),Bb_dn(Nsite,Nsite)
complex(kind=8),intent(INOUT)::D_up(Nsite,Nsite),D_dn(Nsite,Nsite)
complex(kind=8)::Rmat(Ntot,Ntot)
complex(kind=8)::mup(Nsite),ui(Nsite),mup_tilde(Nsite)
real(kind=8)::norm(Nsite)
real(kind=8)::anm1,anm2
integer::j,l,ip
complex(kind=8)::proj(Ntot,Nsite),projb(Ntot,Nsite)
complex(kind=8),external::zdotc
complex(kind=8)::dummy
integer::i_deb,j_deb
logical::i_debug

i_debug=.false.

!make the orbitals in the wave function orthonormal
!call modGS(ph(1:Nsite,1:Nspin(1)),Nsite,Nspin(1),anm1,Rmat(1:Nspin(1),1:Nspin(1)))
!call modGS(ph((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot),Nsite,Nspin(2),anm2,Rmat(1:Nspin(2),1:Nspin(2)))

! |B^{-1,dagger} alpha> ---> |alpha'> =  sum_i=1^N  < u_i | B^{-1,dagger} alpha > |u _i >
! this does not change the ket

proj=zero
projb=zero
do l=1,Nsite
  mup(:)=Bb_up(:,l)
  do ip=1,Nspin(1)
    ui(1:Nsite)=ph(1:Nsite,ip)
    projb(ip,l)=zdotc(Nsite,ui,1,mup,1)    !< u_ip | B^{-1,dagger} alpha >
  enddo
enddo
do l=1,Nsite
  mup(:)=Bf_up(:,l)
  do ip=1,Nspin(1)
    ui(1:Nsite)=ph(1:Nsite,ip)
    proj(ip,l)=zdotc(Nsite,ui,1,mup,1)  !< u_ip | B beta >
  enddo
enddo



!DEBUG
if(i_debug)then
write(*,*)
write(*,*)'phiup 1 = ',dble(ph(1:Nsite,1))
write(*,*)
write(*,*)'phiup 2 = ',dble(ph(1:Nsite,2))
write(*,*)
write(*,*)'phidn 1 = ',dble(ph(Nsite+1:2*Nsite,3))
write(*,*)
write(*,*)'phidn 2 = ',dble(ph(Nsite+1:2*Nsite,4))
write(*,*)
write(*,*)'|B-1T 1> =  ',dble(Bb_up(:,1))
write(*,*)
write(*,*)'|B 1 > =  ',dble(Bf_up(:,1))
write(*,*)
write(*,*)'< u ip | |B-1T 1> = ',dble(projb(:,1))  
write(*,*)
write(*,*)'< u ip | |B  1> = ',dble(proj(:,1))
endif

!Update Bb_up
! |B^{-1,dagger} alpha> ---> |alpha'> =  sum_i=1^N  < u_i | B^{-1,dagger} alpha >   |u _i >
do l=1,Nsite
  do j=1,Nsite
    mup_tilde(j)=zero 
    do ip=1,Nspin(1)
      mup_tilde(j)=mup_tilde(j)+projb(ip,l)*ph(j,ip)
    enddo
  enddo
  do j=1,Nsite
    Bb_up(j,l)=mup_tilde(j)
  enddo
enddo
!if(.false.)then
!compute <alpha' | B beta> and update overlaps

do l=1,Nsite
  do j=1,Nsite
    dummy=zero
!    D_up(j,l)=dummy
    do ip=1,Nspin(1)
 !     D_up(j,l)=D_up(j,l)+conjg(projb(ip,j))*proj(ip,l)
      dummy=dummy+conjg(projb(ip,j))*proj(ip,l)
    enddo
    D_up(j,l)=D_up(j,l)+dummy
  enddo
enddo

!DEBUG
if(i_debug)then
write(*,*)
write(*,*)'|1 prime (alpha) > =  ',dble(Bb_up(:,1))
write(*,*)
write(*,*)'|B 1 > =  ',dble(Bf_up(:,1))
write(*,*)
write(*,*)'<1 prime (alpha)  |B 1 > = ',dble(D_up(1,1))
write(*,*)
endif

!update Bf_up
! | B beta> ---> |beta'> = | B beta> - sum_i=1^N  < u_i | B beta >  |u _i >
!if(.false.)then
do l=1,Nsite
  do j=1,Nsite
!    mup_tilde(j)=Bf_up(j,l)
    mup_tilde(j)=zero
    do ip=1,Nspin(1)
   !   mup_tilde(j)=mup_tilde(j)-proj(ip,l)*ph(j,ip)
      mup_tilde(j)=mup_tilde(j)+proj(ip,l)*ph(j,ip)
    enddo
  enddo
  do j=1,Nsite
    !Bf_up(j,l)=mup_tilde(j)
     Bf_up(j,l)=Bf_up(j,l)-mup_tilde(j)
  enddo
enddo
!endif

!DEBUG
if(i_debug)then
write(*,*)
write(*,*)'|1 prime (beta) > =  ',dble(Bf_up(:,1))
write(*,*)
endif



!---- SPIN DOWN SAME PROCEDURE
do l=1,Nsite
  mup(:)=Bb_dn(:,l)
  do ip=Nspin(1)+1,Ntot
    ui(1:Nsite)=ph(Nsite+1:2*Nsite,ip)
    projb(ip,l)=zdotc(Nsite,ui,1,mup,1)    !< u_ip | B^{-1,dagger} alpha >
  enddo
enddo
do l=1,Nsite
  mup(:)=Bf_dn(:,l)
  do ip=Nspin(1)+1,Ntot
    ui(1:Nsite)=ph(Nsite+1:2*Nsite,ip)
    proj(ip,l)=zdotc(Nsite,ui,1,mup,1)  !< u_ip | B beta >
  enddo
enddo
!Update Bb_dn
! |B^{-1,dagger} alpha> ---> |alpha'> =  sum_i=1^N  < u_i | B^{-1,dagger} alpha
! >   |u _i >
do l=1,Nsite
  do j=1,Nsite
    mup_tilde(j)=zero
    do ip=Nspin(1)+1,Ntot
      mup_tilde(j)=mup_tilde(j)+projb(ip,l)*ph(Nsite+j,ip)
    enddo
  enddo
  do j=1,Nsite
    Bb_dn(j,l)=mup_tilde(j)
  enddo
enddo
!compute <alpha' | B beta> and update overlaps
do l=1,Nsite
  do j=1,Nsite
    dummy=zero
!    D_dn(j,l)=dummy
    do ip=Nspin(1)+1,Ntot
    !  D_dn(j,l)=D_dn(j,l)+conjg(projb(ip,j))*proj(ip,l)
      dummy=dummy+conjg(projb(ip,j))*proj(ip,l)
    enddo
    D_dn(j,l)=D_dn(j,l)+dummy
  enddo
enddo

!update Bf_dn
! | B beta> ---> |beta'> = | B beta> - sum_i=1^N  < u_i | B beta >  |u _i >
do l=1,Nsite
  do j=1,Nsite
!    mup_tilde(j)=Bf_dn(j,l)
    mup_tilde(j)=zero
    do ip=Nspin(1)+1,Ntot
      mup_tilde(j)=mup_tilde(j)+proj(ip,l)*ph(Nsite+j,ip)
    enddo
  enddo
  do j=1,Nsite
    Bf_dn(j,l)=Bf_dn(j,l)-mup_tilde(j)
  enddo
enddo

end subroutine stabilize


subroutine stabilize_c(ph,Bf,Bb,D)
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param
use method_param
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(INOUT)::Bf(2*Nsite,2*Nsite),Bb(2*Nsite,2*Nsite)
complex(kind=8),intent(INOUT)::D(2*Nsite,2*Nsite)
complex(kind=8)::Rmat(Ntot,Ntot)
complex(kind=8)::mup(2*Nsite),ui(2*Nsite),mup_tilde(2*Nsite)
integer::j,l,ip
complex(kind=8)::proj(Ntot,2*Nsite),projb(Ntot,2*Nsite)
complex(kind=8),external::zdotc
complex(kind=8)::dummy
integer::i_deb,j_deb
logical::i_debug

i_debug=.false.

! |B^{-1,dagger} alpha> ---> |alpha'> =  sum_i=1^N  < u_i | B^{-1,dagger} alpha > |u _i >
! this does not change the ket

proj=zero
projb=zero
do l=1,2*Nsite
  mup(:)=Bb(:,l)
  do ip=1,Ntot
    ui(1:2*Nsite)=ph(1:2*Nsite,ip)
    projb(ip,l)=zdotc(2*Nsite,ui,1,mup,1)    !< u_ip | B^{-1,dagger} alpha >
  enddo
enddo
do l=1,2*Nsite
  mup(:)=Bf(:,l)
  do ip=1,Ntot
    ui(1:2*Nsite)=ph(1:2*Nsite,ip)
    proj(ip,l)=zdotc(2*Nsite,ui,1,mup,1)  !< u_ip | B beta >
  enddo
enddo


!Update Bb_up
! |B^{-1,dagger} alpha> ---> |alpha'> =  sum_i=1^N  < u_i | B^{-1,dagger} alpha >   |u _i >
do l=1,2*Nsite
  do j=1,2*Nsite
    mup_tilde(j)=zero 
    do ip=1,Ntot
      mup_tilde(j)=mup_tilde(j)+projb(ip,l)*ph(j,ip)
    enddo
  enddo
  do j=1,2*Nsite
    Bb(j,l)=mup_tilde(j)
  enddo
enddo
!if(.false.)then
!compute <alpha' | B beta> and update overlaps

do l=1,2*Nsite
  do j=1,2*Nsite
    dummy=zero
    do ip=1,Ntot
      dummy=dummy+conjg(projb(ip,j))*proj(ip,l)
    enddo
    D(j,l)=D(j,l)+dummy
  enddo
enddo

!update Bf_up
! | B beta> ---> |beta'> = | B beta> - sum_i=1^N  < u_i | B beta >  |u _i >
!if(.false.)then
do l=1,2*Nsite
  do j=1,2*Nsite
    mup_tilde(j)=zero
    do ip=1,Ntot
      mup_tilde(j)=mup_tilde(j)+proj(ip,l)*ph(j,ip)
    enddo
  enddo
  do j=1,2*Nsite
     Bf(j,l)=Bf(j,l)-mup_tilde(j)
  enddo
enddo
!endif

end subroutine stabilize_c
