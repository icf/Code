!----------------------------------------------------------------
!The subroutine is the forward propagation before progapagte back
!----------------------------------------------------------------
subroutine back_propag(sample,i_pop,i_GS,i_l)
use phi_param
use phiT_param
use lattice_param
use mc_loop_param
use model_param
use method_param
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer,intent(IN)::sample
integer,intent(INOUT)::i_pop,i_GS
integer,intent(IN)::i_l
integer::i,j,k
real(kind=8)::t0,t2,t3


!Get the backup of phi0
!call zcopy(2*Nsite*Ntot*Nwalkers,phi(1,1,1),1,phi0(1,1,1),1)
call copy_wf_W_dc(phi(1,1,1),phi0(1,1,1))

!prepare for free projection
if(crn.LT.0.d0) then
  call Modified_GS();i_GS=0
  crn=1.d0
  do i=1,Nwalkers,1
     rx(i)=rx(i)/tot_imp(i)
  end do
end if

do i_back=1,Nbk(1),1
   call one_propagation(i_pop,i_GS)
end do


!prepare for constraint path
if(crn.GT.0.d0) then
  call Modified_GS();i_GS=0
  crn=-1.d0
  call cal_imp_ovlap()
  do i=1,Nwalkers,1
     rx(i)=rx(i)*tot_imp(i)
  end do
end if

do i_back=Nbk(1)+1,Nstps_fwd,1
   call one_propagation(i_pop,i_GS)
end do

!project back and measure
if(max_crn.GE.0) then
  call add_measure(sample,i_l)
else
  call freep_measure(sample,i_l)
end if


end subroutine back_propag



!----------------------------------------
!This subroutine propagate the phT to phL for BCS
!----------------------------------------
subroutine bcs_back_progF(i,F)
use param
use lattice_param
use model_param
use method_param
use phiT_param
use mc_loop_param
use phi_param
use project_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(OUT)::F(Nsite,Nsite)
real(kind=8)::norm
integer::j,k,i_b,i_deb,j_deb

call zcopy(Nsite*Nsite,Fpairing(1,1),1,F(1,1),1)
!do i_b=Nstps_fwd,1,-1
!DEBUG
!  write(*,*)
!  write(*,*)'bcs_back_progF'
!  write(*,*)
!  write(*,*)'i_b = ',i_b
!  write(*,*)
!  write(*,*)'Matrix F = '
!  do i_deb=1,Nsite
!    write(*,*)(F(i_deb,j_deb),j_deb=1,Nsite)
!  enddo
!  write(*,*)

!  call bcs_back_one(i_b,i,F)

!DEBUG
!  write(*,*)'Matrix F after back one = '
!  do i_deb=1,Nsite
!    write(*,*)(F(i_deb,j_deb),j_deb=1,Nsite)
!  enddo
!  write(*,*)
!  write(*,*)
!  stop 'deb'

!  norm=0.d0
!  do j=1,Nsite
!    do k=1,Nsite
!      norm=norm+abs(F(k,j))**2
!    enddo
!  enddo
!  norm=abs(F(1,1))
!  do j=1,Nsite
!    do k=1,Nsite
!      F(k,j)=F(k,j)/norm
!    enddo
!  enddo

!enddo

!DEBUG
!  write(*,*)'Matrix F after back one = '
!  do i_deb=1,Nsite
!    write(*,*)(F(i_deb,j_deb),j_deb=1,Nsite)
!  enddo
!  write(*,*)
!  write(*,*)

end subroutine bcs_back_progF


!----------------------------------------
!This subroutine propagate the phi to phL
!---------------------------------------- 
subroutine back_prog_phL(i,phL,coe)
use param
use lattice_param
use model_param
use method_param
use phiT_param
use mc_loop_param
use phi_param
use project_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(OUT)::phL(2*Nsite,Ntot,Dtot)
complex(kind=8),intent(OUT)::coe(Dtot)
complex(kind=8)::phtmp(2*Nsite,Ntot)
integer::i_gs
integer::j,k,i_b

do k=1,Dtot,1
   !call zgemm('C','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phiT(1,1,k),2*Nsite,zero,phL(1,1,k),2*Nsite)
   call dk_to_ph_dc(exp_halfK,phiT(1,1,k),phL(1,1,k))
end do
call zcopy(Dtot,coe_multi(1),1,coe(1),1)


i_gs=0
do i_b=Nstps_fwd,1,-1
   call back_one(i_b,i,phL)
   i_gs=i_gs+1
   if(i_gs.eq.StepforGram) then
     i_gs=0
     call gs_phL(phL,coe)
   end if
end do

do k=1,Dtot,1
   call copy_wf_dc(phL(1,1,k),phtmp(1,1))
   call dk_to_ph_dc(exp_mhalfK,phtmp(1,1),phL(1,1,k))
end do
end subroutine back_prog_phL


!------------------------------------------------------------
!This subroutine do one back propagation (expK. expV)^+ |BCS>
!------------------------------------------------------------

subroutine bcs_back_one(i_b,i,F)
use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
implicit none
integer,intent(IN)::i_b,i
complex(kind=8),intent(INOUT)::F(Nsite,Nsite)
complex(kind=8)::Ftemp(Nsite,Nsite),Buptemp(Nsite,Nsite),Bdntemp(Nsite,Nsite)
complex(kind=8)::Bup(Nsite,Nsite),Bdn(Nsite,Nsite),Kappa(Nsite,Nsite)
complex(kind=8)::explr_up,explr_dn,lr_up,lr_dn
!complex(kind=8)::Prova(Nsite,Nsite)
real(kind=8)::x
integer::aux
integer::j,k,l,i_deb,j_deb

!DEBUG
!write(*,*)
!write(*,*)'In bcs_back_one '
!write(*,*)

do j=1,Nsite,1
  x=back_store(j,i_b,i)
  aux=NINT(x)
  explr_up=expln_up(aux)
  explr_dn=expln_dn(aux)
  do k=1,Nsite,1
    Bup(k,j)=exp_halfK(k,j)*(explr_up+one)
    Bdn(k,j)=exp_halfK(k,j)*(explr_dn+one)
  enddo
enddo
Buptemp=Bup
Bdntemp=Bdn

!DEBUG
!write(*,*)'K/2 '
!do i_deb=1,Nsite
!  write(*,*)(Buptemp(i_deb,j_deb),j_deb=1,Nsite)
!enddo
!write(*,*)'exp_halfK '
!do i_deb=1,Nsite
!  write(*,*)(exp_halfK(i_deb,j_deb),j_deb=1,Nsite)
!enddo
!write(*,*)
!write(*,*)'one , zero ',one,zero
!write(*,*)

Kappa(:,:)=exp_halfK(1:Nsite,1:Nsite)
call zgemm('N','N',Nsite,Nsite,Nsite,one,Buptemp,Nsite,Kappa,Nsite,zero,Bup,Nsite)
call zgemm('N','N',Nsite,Nsite,Nsite,one,Bdntemp,Nsite,Kappa,Nsite,zero,Bdn,Nsite)


!DEBUG
!write(*,*)'Bup '
!do i_deb=1,Nsite
!  write(*,*)(Bup(i_deb,j_deb),j_deb=1,Nsite)
!enddo
!write(*,*)
!write(*,*)'Bdn '
!do i_deb=1,Nsite
!  write(*,*)(Bdn(i_deb,j_deb),j_deb=1,Nsite)
!enddo
!write(*,*)


Bup=conjg(Bup)
Bdn=conjg(Bdn)
call zgemm('N','N',Nsite,Nsite,Nsite,one,F,Nsite,Bdn,Nsite,zero,Ftemp,Nsite)
call zgemm('T','N',Nsite,Nsite,Nsite,one,Bup,Nsite,Ftemp,Nsite,zero,F,Nsite)

end subroutine bcs_back_one

!------------------------------------------------------------
!This subroutine do one back propagation (expK. expV)^+ |phL>
!------------------------------------------------------------
subroutine back_one(i_b,i,phL)
use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
implicit none
integer,intent(IN)::i_b,i
complex(kind=8),intent(INOUT)::phL(2*Nsite,Ntot,Dtot)
complex(kind=8)::phtmp(2*Nsite,Ntot)
complex(kind=8)::explr_up,explr_dn,lr_up,lr_dn
real(kind=8)::x
integer::aux
integer::j,k,l


if(Nbands.eq.3)then

 do j=1,Nsite,1
 
   x=back_store(j,i_b,i)
   aux=NINT(x)

   if(j.le.Nbravais)then
     explr_up=explnd_up(aux)
     explr_dn=explnd_dn(aux)
   elseif(j.le.2*Nbravais)then
     explr_up=explnx_up(aux)
     explr_dn=explnx_dn(aux)
   else
     explr_up=explny_up(aux)
     explr_dn=explny_dn(aux)
   endif

   if(dtype.EQ.'c') then
     do k=1,Dtot,1
       do l=1,Ntot,1
          phL(j,l,k)=phL(j,l,k)*conjg(explr_up+one)
          phL(j+Nsite,l,k)=phL(j+Nsite,l,k)*conjg(explr_dn+one)
       end do
     end do
   else if(dtype.EQ.'d') then
     do k=1,Dtot,1
       do l=1,Nspin(1),1
          phL(j,l,k)=phL(j,l,k)*conjg(explr_up+one)
       end do
       do l=Nspin(1)+1,Ntot,1
          phL(j+Nsite,l,k)=phL(j+Nsite,l,k)*conjg(explr_dn+one)
       end do
     end do
   end if

 enddo


elseif(Nbands.eq.1)then

 do j=1,Nsite,1

   x=back_store(j,i_b,i)
   if(i_b.LE.Nbk(1).AND.i_b.GE.1) then
     if(kcrn.eq.1) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.2) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.3) then
       lr_up=x*sqrt(dcmplx(dt*onsitU))
       lr_up=lr_up-(dt*onsitU/2.d0)

       lr_dn=-x*sqrt(dcmplx(dt*onsitU))
       lr_dn=lr_dn-(dt*onsitU/2.d0)

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.4) then
       lr_up=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_up=lr_up+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       lr_dn=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_dn=lr_dn+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else
       write(*,*) "Something is wrong with kcrn input:",kcrn
       call mystop
     end if
   else if(i_b.GT.Nbk(1).AND.i_b.LE.Nstps_fwd) then
     aux=NINT(x)
     explr_up=expln_up(aux)
     explr_dn=expln_dn(aux)
   else
     write(*,*) "Something is wrong of i_b",i_b
     call mystop
   end if


   if(dtype.EQ.'c') then
     do k=1,Dtot,1
       do l=1,Ntot,1
          phL(j,l,k)=phL(j,l,k)*conjg(explr_up+one)
          phL(j+Nsite,l,k)=phL(j+Nsite,l,k)*conjg(explr_dn+one)
       end do
     end do
   else if(dtype.EQ.'d') then
     do k=1,Dtot,1
       do l=1,Nspin(1),1
          phL(j,l,k)=phL(j,l,k)*conjg(explr_up+one)
       end do
       do l=Nspin(1)+1,Ntot,1
          phL(j+Nsite,l,k)=phL(j+Nsite,l,k)*conjg(explr_dn+one)
       end do
     end do
   end if

 end do 

endif


do k=1,Dtot,1
   !call zcopy(2*Nsite*Ntot,phL(1,1,k),1,phtmp(1,1),1)
   call copy_wf_dc(phL(1,1,k),phtmp(1,1))
   !call
   !zgemm('C','N',2*Nsite,Ntot,2*Nsite,one,exp_K,2*Nsite,phtmp(1,1),2*Nsite,zero,phL(1,1,k),2*Nsite)
   call dk_to_ph_dc(exp_K,phtmp(1,1),phL(1,1,k))
end do
end subroutine back_one


!--------------------------------------------------
!Do the modified gs for the phl, and normalized coe
!--------------------------------------------------
subroutine gs_phL(phL,coe)
use lattice_param
use model_param
use phiT_param
implicit none
complex(kind=8),intent(INOUT)::phL(2*Nsite,Ntot,Dtot)
complex(kind=8),intent(INOUT)::coe(Dtot)
real(kind=8)::anm,anm1,anm2
complex(kind=8)::Rmat(Ntot,Ntot)
integer::k

do k=1,Dtot,1
   if(dtype.EQ.'c') then
     call modGS(phL(1,1,k),2*Nsite,Ntot,anm,Rmat)
   else if(dtype.EQ.'d') then
     call modGS(phL(1:Nsite,1:Nspin(1),k),Nsite,Nspin(1),anm1,Rmat(1:Nspin(1),1:Nspin(1)))
     call modGS(phL((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,k),Nsite,Nspin(2),anm2,Rmat(1:Nspin(2),1:Nspin(2)))
     anm=anm1*anm2
   end if
   coe(k)=coe(k)/anm
end do

call norm_array(Dtot,coe)
end subroutine gs_phL


!-------------------------
!normlize the array of coe
!-------------------------
subroutine norm_array(Dtot,coe)
implicit none
integer,intent(IN)::Dtot
complex(kind=8),intent(INOUT)::coe(Dtot)
complex(kind=8)::norm
real(kind=8),external::dznrm2
complex(kind=8)::one=dcmplx(1d0,0d0)
norm=one/dcmplx(dznrm2(Dtot,coe(1),1))
call zscal(Dtot,norm,coe(1),1)
end subroutine norm_array
