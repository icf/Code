!One propagation in update, include GS,pop control
subroutine one_propagation(i_pop,i_GS)
use mc_loop_param
use method_param
!DEBUG
!use phi_param
implicit none
integer,intent(INOUT)::i_pop,i_GS

 i_pop=i_pop+1;i_GS=i_GS+1


!DEBUG
!  write(*,*)'weight before update ',weight

 !propagate one step in MC
 call K_V_update()

!  !DEBUG
!  write(*,*)'weight after update ',weight


 !Do periodically Modified GS
 !While free-projection it is done
 !in the population control step
 if(crn.LT.0.d0) then
   if(i_GS.EQ.StepforGram) then
     call Modified_GS()
     i_GS=0
   end if
 end if


 !Do periodically population control and adjust ET
 !If crn.GT.0.d0  do modified gs
 if(i_pop.EQ.PopContrlstep) then
   if(crn.GT.0.d0) then
     call Modified_GS()
     i_GS=0
   end if
   call PopControl()
   i_pop=0
 end if

end subroutine one_propagation


!------------------------------------------
!Exp_U.Exp_K.phi,it is a first order update
!------------------------------------------
subroutine K_V_update()
use phi_param
use mc_loop_param
use project_param
use adET_param
use model_param
use phiT_param
implicit none
complex(kind=8)::ovlp(Ntot,Ntot,Dtot)
logical::inv_flag=.true.
integer::i,j


!DEBUG
!write(*,*)'weight before update ',weight

do i=1,Nwalkers,1
  !first do K then do U, use the information of inverse
  call K_phi(2,i,ovlp,inv_flag)
  !write(*,*) "imptant function is:",tot_imp(i)
  !write(*,*) "K done",i;pause
  if(I_wavefun.eq.1)then 
    call U_phi(i,ovlp)
  elseif(I_wavefun.eq.2)then
    call BCS_U_phi(i,ovlp)
  endif
  !write(*,*) "imptant function is:",tot_imp(i)
  !write(*,*) "U done",i
  if(weight(i).le.0.d0) cycle

  dlogw(i)=dlogw(i)+dt*ET
  weight(i)=exp(dlogw(i))
end do


!DEBUG
!write(*,*)'weight after update ',weight


end subroutine K_V_update



!---------------------------------------------------------------------------------
!update exp_mhalfK or exp_K to the wave function, whole walkers,used for beginning
!---------------------------------------------------------------------------------
subroutine K_phi_W(input)
use mc_loop_param
use model_param
use phiT_param
implicit none
integer,intent(IN)::input !input=1 exp_half_K; input=2 exp_K
complex(kind=8)::ovlp(Ntot,Ntot,Dtot)
logical::inv_flag=.false.
integer::i
do i=1,Nwalkers,1
   call K_phi(input,i,ovlp,inv_flag)
end do
end subroutine K_phi_W


!---------------------------------------------
!update exp_mhalfK or exp_K to the wave function
!---------------------------------------------
subroutine K_phi(input,i,ovlp,inv_flag)
use param
use mc_loop_param
use phi_param
use lattice_param
use model_param
use project_param
use method_param
use phiT_param
implicit none
integer,intent(IN)::input !input=1 exp_half_K; input=2 exp_K
integer,intent(IN)::i !the label of walkers
complex(kind=8),intent(OUT)::ovlp(Ntot,Ntot,Dtot)
logical,intent(IN)::inv_flag
complex(kind=8)::temp(2*Nsite,Ntot)
complex(kind=8)::tot_tmp

if(weight(i).le.0.d0) return

!call zcopy(2*Nsite*Ntot,phi(1,1,i),1,temp(1,1),1)
call copy_wf_dc(phi(1,1,i),temp(1,1))
if(input.eq.1) then
   !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_mhalfK,2*Nsite,temp,2*Nsite,zero,phi(1,1,i),2*Nsite)
   call k_to_ph_dc(exp_mhalfK,temp(1,1),phi(1,1,i))
else if(input.eq.2) then
   !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_K,2*Nsite,temp,2*Nsite,zero,phi(1,1,i),2*Nsite)
   call k_to_ph_dc(exp_K,temp(1,1),phi(1,1,i))
else
   write(*,*) "Wrong in put in subroutine K_phi:",input
   write(*,*) "input=1 exp_half_K; input=2 exp_K"
   call mystop
end if

!If CPMC , then calculate the impfunction
if(crn.LT.0.d0) then

   tot_tmp=tot_imp(i)
   if(inv_flag) then
     call get_imp_inv(i,ovlp)
   else
     call get_imp(i)
   end if
   if(dble(tot_imp(i)/tot_tmp).LE.0.d0) then
      weight(i)=0.d0
      dlogw(i)=-1d100
      return
   end if

   weight(i)=weight(i)*dble(tot_imp(i)/tot_tmp)
   dlogw(i)=dlogw(i)+dlog(dble(tot_imp(i)/tot_tmp))

   if(weight(i).le.0.d0)  then
      write(*,*) "something is wrong in K_phi."
      call mystop
   end if

end if
end subroutine K_phi




!---------------------------------
!update exp_U to the wave function
!---------------------------------
subroutine U_phi(i,ovlp)
use param
use mc_loop_param
use phi_param
use lattice_param
use model_param
use phiT_param
use method_param
use project_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(INOUT)::ovlp(Ntot,Ntot,Dtot)
integer::j,k
real(kind=8)::avgnu,avgnd
!test
!inver=zero
!end test

if(weight(i).le.0.d0) return

!Get the inverse of <phiT|phi>
!if(crn.LT.0.d0) call get_inverse(i,inver(1,1,1))

if(crn.GT.0.d0) then
  !set the background
  if(bgset.eq.0) then !use mean field
  else if(bgset.eq.1) then !use dynamic field from walker
    call modGS_i(i)
    do j=1,Nsite,1
       if(dtype.EQ.'c') then
         avgnu=zero
         do k=1,Ntot,1
            avgnu=avgnu+(abs(phi(j,k,i)))**2
         end do
         avgnd=zero
         do k=1,Ntot,1
            avgnd=avgnd+(abs(phi(j+Nsite,k,i)))**2
         end do
       else if(dtype.EQ.'d') then
         avgnu=zero
         do k=1,Nspin(1),1
            avgnu=avgnu+(abs(phi(j,k,i)))**2
         end do
         avgnd=zero
         do k=Nspin(1)+1,Ntot,1
            avgnd=avgnd+(abs(phi(j+Nsite,k,i)))**2
         end do
       end if
       if(kcrn.eq.3.or.kcrn.eq.1) then
          ng(j)=dble(avgnu-avgnd)
       else if(kcrn.eq.4.or.kcrn.eq.2) then
          ng(j)=dble(avgnu+avgnd)
       end if
    end do
  else
    write(*,*) "something is wrong with bgset:",bgset
    call mystop
  end if
end if


if(Nbands.eq.1)then
  do j=1,Nsite,1
     call single_U_phi(i,j,ovlp(1,1,1))
     if(weight(i).le.0.d0) exit
  end do
elseif(Nbands.eq.3)then
  do j=1,Nsite,1
    call single_U_phi_3b(i,j,ovlp(1,1,1))
    if(weight(i).le.0.d0) exit
  enddo
!  stop 'DEBUG'
else
  write(*,*) "something is wrong with Nbands:",Nbands
  call mystop
endif
end subroutine U_phi



!-------------------------------
!Get the inverse of <phiT|phi_i>
!-------------------------------
subroutine get_inverse(i,inver)
use phi_param
use lattice_param
use model_param
use phiT_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(OUT)::inver(Ntot,Ntot,Dtot)
integer::k
do k=1,Dtot,1
   !call deter_overlap(2*Nsite,Ntot,phiT(1,1,k),phi(1,1,i),inver(1,1,k))
   !call inverse(inver(1:Ntot,1:Ntot,k),Ntot)
   call over_lap_inv_dc(phiT(1,1,k),phi(1,1,i),inver(1,1,k)) 
end do
end subroutine get_inverse


!------------------------------------
!update exp_Uj to the wave function i
!------------------------------------
subroutine single_U_phi_3b(i,j,inver)
use param
use rand_num
use phi_param
use model_param
use lattice_param
use phiT_param
use project_param
use method_param
use mpi_serial_param
implicit none
integer,intent(IN)::i,j
complex(kind=8),intent(INOUT)::inver(Ntot,Ntot,Dtot)
complex(kind=8)::explr_up,explr_dn,lr_up,lr_de
complex(kind=8)::tot_weig
complex(kind=8)::nj(2,2,Dtot)
complex(kind=8)::rat_l,rat_h(2),ratio(2)!Used to update the weight
real(kind=8)::fR(2),prefac
!complex(kind=8)::inv(2,2)
complex(kind=8)::um(Ntot,2,Dtot),vm(2,Ntot,Dtot)
!complex(kind=8)::tmp(2,Ntot),tmp_p(Ntot,Ntot),tmp_pp(Ntot,Ntot)
real(kind=8)::x
real(kind=8)::p(2),Normp
real(kind=8)::U_loc
complex(kind=8)::e_up_loc(2),e_dn_loc(2)
complex(kind=8)::gama_loc
complex(kind=8)::r(2)
integer::aux
integer::k,m,n


   if(j.le.Nbravais)then
     U_loc=onsitUd
     gama_loc=gamad
     e_up_loc(:)=explnd_up(:)
     e_dn_loc(:)=explnd_dn(:)
   elseif(j.le.2*Nbravais)then
     U_loc=onsitUp
     gama_loc=gamax
     e_up_loc(:)=explnx_up(:)
     e_dn_loc(:)=explnx_dn(:)
   else
     U_loc=onsitUp
     gama_loc=gamay
     e_up_loc(:)=explny_up(:)
     e_dn_loc(:)=explny_dn(:)
   endif
  
   !Get the um and vm, i is a walker, j is a basis element
   call get_um_vm(i,j,um,vm,inver)

   !get the green's function matrix nj(2,2,Dtot)
   call get_green_matrix(um,vm,nj)

   !Get two possible new tot_imp
   rat_l=tot_imp(i)
   do m=1,2,1
      rat_h(m)=zero
      do k=1,Dtot,1
         rat_h(m)=rat_h(m)+conjg(coe_multi(k))*impfunc(k,i)  &
                 & *(e_up_loc(m)*e_dn_loc(m)*(nj(1,1,k)*nj(2,2,k)-nj(1,2,k)*nj(2,1,k))  &
                 &  +e_up_loc(m)*nj(1,1,k)+e_dn_loc(m)*nj(2,2,k)+one)
      end do
      ratio(m)=rat_h(m)/rat_l
   end do

   !Update the total normalization factor
   if(dcp.eq.'S') then
     do k=1,2,1
        fR(k)=dble(ratio(k))
        if(fR(k).le.0.d0) fR(k)=0.d0  !We can add mirror correction here. 
     end do
   else if(dcp.eq.'C') then
     !We need another weight change for the density HS transformation
     do k=1,2,1
        fR(k)=dble(ratio(k)*exp(dt*U_loc/2.d0-dble(gama_loc)*dble((-1)**k)))
        if(fR(k).le.0.d0) fR(k)=0.d0
     end do
   end if

   prefac=fR(1)+fR(2) !normalization factor,No the true normalization.Need to
                      !multipuly 0.5

   if(prefac.le.0.d0) then
      weight(i)=0.d0
      dlogw(i)=-1d100
      return
   end if

   weight(i)=weight(i)*prefac*0.5d0   !We must mutipuly 0.5d0 here so that 
                                      !the weight will not grow up too
                                      !large.
   dlogw(i)=dlogw(i)+dlog(prefac*0.5d0)


   !Choose the auxiliary flied.
    if(fR(1)/prefac.ge.rndm()) then
      aux=1
      x=1.d0
    else
      aux=2
      x=2.d0
    end if

    explr_up=e_up_loc(aux)
    explr_dn=e_dn_loc(aux)


   !Change the impfunction.
   do k=1,Dtot,1
      impfunc(k,i)=impfunc(k,i)*(explr_up*explr_dn*(nj(1,1,k)*nj(2,2,k)-nj(1,2,k)*nj(2,1,k)) &
                 & +explr_up*nj(1,1,k)+explr_dn*nj(2,2,k)+one)
   end do
   tot_imp(i)=rat_h(aux)

   if(j.LT.Nsite) then
     !Get inverse of the overlap
     do k=1,Dtot,1
       call update_inverse_3b(i,j,k,nj(1,1,k),aux,um(1,1,k),vm(1,1,k),inver(1,1,k))

     end do
   else if(j.eq.Nsite) then
     !We do not need to update the inverse
   else
     write(*,*) "Something is wrong with the onsit j input:",j
     call mystop
   end if


!Store the propogation Ising spin for back propogation.
 if(back_pro) then
    back_store(j,i_back,i)=x
 end if

!DEB
!   write(*,*)
!   write(*,*)'e_up_loc ....'
!   write(*,*)'basis element ',j
!   write(*,*)'U_loc= ',U_loc
!   write(*,*)'gama_loc= ',gama_loc
!   write(*,*)'e_up_loc=',e_up_loc
!   write(*,*)'e_dn_loc=',e_dn_loc
!   write(*,*)'explr_up+one=',explr_up+one
!   write(*,*)'explr_dn+one=',explr_dn+one
!   write(*,*)
!   write(*,*)



!get the new determinate
 if(dtype.EQ.'c') then
   do k=1,Ntot,1
      phi(j,k,i)=phi(j,k,i)*(explr_up+one)
      phi(j+Nsite,k,i)=phi(j+Nsite,k,i)*(explr_dn+one)
   end do
 else if(dtype.EQ.'d') then
   do k=1,Nspin(1),1
      phi(j,k,i)=phi(j,k,i)*(explr_up+one)
   end do
   do k=Nspin(1)+1,Ntot,1
      phi(j+Nsite,k,i)=phi(j+Nsite,k,i)*(explr_dn+one)
   end do
 end if


end subroutine single_U_phi_3b



!------------------------------------
!update exp_Uj to the wave function i
!------------------------------------
subroutine single_U_phi(i,j,inver)
use param
use rand_num
use phi_param
use model_param
use lattice_param
use phiT_param
use project_param
use method_param
use mpi_serial_param
implicit none
integer,intent(IN)::i,j
complex(kind=8),intent(INOUT)::inver(Ntot,Ntot,Dtot)
complex(kind=8)::explr_up,explr_dn,lr_up,lr_dn
complex(kind=8)::tot_weig
complex(kind=8)::nj(2,2,Dtot)
complex(kind=8)::rat_l,rat_h(2),ratio(2)!Used to update the weight
real(kind=8)::fR(2),prefac
!complex(kind=8)::inv(2,2)
complex(kind=8)::um(Ntot,2,Dtot),vm(2,Ntot,Dtot)
!complex(kind=8)::tmp(2,Ntot),tmp_p(Ntot,Ntot),tmp_pp(Ntot,Ntot)
real(kind=8)::x
real(kind=8)::p(2),Normp
complex(kind=8)::r(2)
integer::aux
integer::k,m,n


!Get the update explr_up and explr_dn with update the weight.
 if(crn.GT.0.d0) then

   if(kcrn.eq.1) then !Use discrete spin decouple
     x=rndm()

     p(1)=abs(exp(-gamaf*ng(j)))
     p(2)=abs(exp(gamaf*ng(j)))
     r(1)=exp(-gamaf*ng(j))/p(1)
     r(2)=exp(gamaf*ng(j))/p(2)
     Normp=p(1)+p(2)

     if(x.LT.p(1)/Normp) then
       lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
       lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       rx(i)=rx(i)*r(1)
       tot_weig=exp(gamaf*ng(j))
     else
       lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
       lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       rx(i)=rx(i)*r(2)
       tot_weig=exp(-gamaf*ng(j))
     end if

     weight(i)=weight(i)*Normp*0.5*abs(tot_weig)
     dlogw(i)=dlogw(i)+dlog(Normp*0.5)+dlog(abs(tot_weig))
     rx(i)=rx(i)*tot_weig/abs(tot_weig)

     explr_up=exp(lr_up)-one
     explr_dn=exp(lr_dn)-one

   else if(kcrn.eq.2) then !Use discrete charge decouple
     x=rndm()

     p(1)=abs(exp((1-ng(j))*gamaf))
     p(2)=abs(exp((ng(j)-1)*gamaf))
     r(1)=exp((1-ng(j))*gamaf)/p(1)
     r(2)=exp((ng(j)-1)*gamaf)/p(2)
     Normp=p(1)+p(2)

     if(x.LT.p(1)/Normp) then
       lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
       lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       rx(i)=rx(i)*r(1)
       tot_weig=exp(gamaf*ng(j))
     else
       lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
       lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       rx(i)=rx(i)*r(2)
       tot_weig=exp(-gamaf*ng(j))
     end if

     weight(i)=weight(i)*abs(exp(dt*OnsitU/2.d0))*Normp*0.5*abs(tot_weig)
     dlogw(i)=dlogw(i)+dt*OnsitU/2.d0+dlog(Normp*0.5)+dlog(abs(tot_weig))
     rx(i)=rx(i)*tot_weig/abs(tot_weig)

     explr_up=exp(lr_up)-one
     explr_dn=exp(lr_dn)-one

   else if(kcrn.eq.3) then !Use continous spin decouple
     x=gauss_rndm()
     !lr_up=x*sqrt(dcmplx(dt*onsitU))
     !lr_up=lr_up-(dt*onsitU/2.d0)

     !lr_dn=-x*sqrt(dcmplx(dt*onsitU))
     !lr_dn=lr_dn-(dt*onsitU/2.d0)

     lr_up=x*sqrt(dcmplx(dt*onsitU))
     lr_up=lr_up-(dt*onsitU/2.d0)+dt*onsitU*ng(j)

     lr_dn=-x*sqrt(dcmplx(dt*onsitU))
     lr_dn=lr_dn-(dt*onsitU/2.d0)-dt*onsitU*ng(j)

     tot_weig=-1.d0*x*sqrt(dcmplx(dt*onsitU))*ng(j)
     tot_weig=tot_weig-0.5d0*dt*onsitU*(ng(j)**2)

     weight(i)=weight(i)*abs(exp(tot_weig))
     dlogw(i)=dlogw(i)+dlog(abs(exp(tot_weig)))
     rx(i)=rx(i)*exp(tot_weig)/abs(exp(tot_weig))


     explr_up=exp(lr_up)-one
     explr_dn=exp(lr_dn)-one

   else if(kcrn.eq.4) then !Use continous charge decouple
     x=gauss_rndm()
     lr_up=x*sqrt(dcmplx(-1.d0*dt*onsitU))
     lr_up=lr_up+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

     lr_dn=x*sqrt(dcmplx(-1.d0*dt*onsitU))
     lr_dn=lr_dn+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

     tot_weig=-1.d0*x*sqrt(dcmplx(-1.d0*dt*onsitU))*ng(j)
     tot_weig=tot_weig+0.5d0*dt*onsitU*(ng(j)**2)

     weight(i)=weight(i)*abs(exp(tot_weig))
     dlogw(i)=dlogw(i)+dlog(abs(exp(tot_weig))) 
     rx(i)=rx(i)*exp(tot_weig)/abs(exp(tot_weig))  

     explr_up=exp(lr_up)-one
     explr_dn=exp(lr_dn)-one
   else
     write(*,*) "Something is wrong with kcrn input:",kcrn
     call mystop
   end if

 else

   !Get the um and vm, i is a walker, j is a basis element
   call get_um_vm(i,j,um,vm,inver)

   !get the green's function matrix nj(2,2,Dtot)
   call get_green_matrix(um,vm,nj)


   !Get <phiT|nj(sigma,sigma')|phi>
   !do k=1,Dtot,1
   !   !call cal_cidcj_withovlpinv(2*Nsite,Ntot,phiT(1,1,k),phi(1,1,i), &
   !   !       & inver(1,1,k),nj(1,1,k),j,j)
   !   !call cal_cidcj_withovlpinv(2*Nsite,Ntot,phiT(1,1,k),phi(1,1,i), &
   !   !       & inver(1,1,k),nj(1,2,k),j,j+Nsite)
   !   !call cal_cidcj_withovlpinv(2*Nsite,Ntot,phiT(1,1,k),phi(1,1,i), &
   !   !       & inver(1,1,k),nj(2,1,k),j+Nsite,j)
   !   !call cal_cidcj_withovlpinv(2*Nsite,Ntot,phiT(1,1,k),phi(1,1,i), &
   !   !       & inver(1,1,k),nj(2,2,k),j+Nsite,j+Nsite)
   !   !if(rank.eq.1) write(*,*) nj(1:2,1:2,k)
   !   call cal_cidcj_withovlpinv_dc(phiT(1,1,k),phi(1,1,i),inver(1,1,k),nj(1,1,k),j,j)
   !   call cal_cidcj_withovlpinv_dc(phiT(1,1,k),phi(1,1,i),inver(1,1,k),nj(1,2,k),j,j+Nsite)
   !   call cal_cidcj_withovlpinv_dc(phiT(1,1,k),phi(1,1,i),inver(1,1,k),nj(2,1,k),j+Nsite,j)
   !   call cal_cidcj_withovlpinv_dc(phiT(1,1,k),phi(1,1,i),inver(1,1,k),nj(2,2,k),j+Nsite,j+Nsite)
   !   !if(rank.eq.1) write(*,*) nj(1:2,1:2,k)
   !   !call mystop
   !end do
   !if(rank.eq.1) write(*,*) phiT(1:2*Nsite,1:Ntot,1)
   !call mystop
   

   !Get two possible new tot_imp
   rat_l=tot_imp(i)
   do m=1,2,1
      rat_h(m)=zero
      do k=1,Dtot,1
         rat_h(m)=rat_h(m)+conjg(coe_multi(k))*impfunc(k,i)  &
                 & *(expln_up(m)*expln_dn(m)*(nj(1,1,k)*nj(2,2,k)-nj(1,2,k)*nj(2,1,k))  &
                 & +expln_up(m)*nj(1,1,k)+expln_dn(m)*nj(2,2,k)+one)
      end do
      ratio(m)=rat_h(m)/rat_l
   end do
   
   !Update the total normalization factor
   if(dcp.eq.'S') then
     do k=1,2,1
        fR(k)=dble(ratio(k)) 
        if(fR(k).le.0.d0) fR(k)=0.d0  !We can add mirror correction here. 
     end do
   else if(dcp.eq.'C') then
     !We need another weight change for the density HS transformation
     do k=1,2,1
       !fR(k)=dble(ratio(k))
        fR(k)=dble(ratio(k)*exp(dt*OnsitU/2.d0-dble(gama)*dble((-1)**k)))
        if(fR(k).le.0.d0) fR(k)=0.d0
     end do
   end if

   prefac=fR(1)+fR(2) !normalization factor,No the true normalization.Need to
                      !multipuly 0.5

   if(prefac.le.0.d0) then
      weight(i)=0.d0
      dlogw(i)=-1d100
      return
   end if

   weight(i)=weight(i)*prefac*0.5d0   !We must mutipuly 0.5d0 here so that 
                                      !the weight will not grow up too
                                      !large ... no it is p(x) from HS
   dlogw(i)=dlogw(i)+dlog(prefac*0.5d0)


   !Choose the auxiliary flied.
    if(fR(1)/prefac.ge.rndm()) then
      aux=1
      x=1.d0
    else
      aux=2
      x=2.d0
    end if


    explr_up=expln_up(aux)
    explr_dn=expln_dn(aux)


   !Change the impfunction.
   do k=1,Dtot,1
      impfunc(k,i)=impfunc(k,i)*(expln_up(aux)*expln_dn(aux)*(nj(1,1,k)*nj(2,2,k)-nj(1,2,k)*nj(2,1,k))  &
                 & +expln_up(aux)*nj(1,1,k)+expln_dn(aux)*nj(2,2,k)+one)
   end do
   tot_imp(i)=rat_h(aux)

   if(j.LT.Nsite) then
     !Get inverse of the overlap
     do k=1,Dtot,1

        !old code
        !!Get the inv(2,2)
        !inv(1,1)=one+expln_up(aux)*nj(1,1,k)
        !inv(1,2)=expln_up(aux)*nj(2,1,k)
        !inv(2,1)=expln_dn(aux)*nj(1,2,k)
        !inv(2,2)=one+expln_dn(aux)*nj(2,2,k)
        !call inv2(inv(1,1))
        ! 
        !!Get vm(2,Ntot)
        !do m=1,Ntot,1
        !   vm(1,m)=phi(j,m,i)*expln_up(aux)
        !   vm(2,m)=phi(j+Nsite,m,i)*expln_dn(aux)
        !end do 
        ! 
        !!Get um(Ntot,2)
        !um(1:Ntot,1)=Conjg(phiT(j,1:Ntot,k))
        !um(1:Ntot,2)=Conjg(phiT(j+Nsite,1:Ntot,k))
        !
        !!tmp=inv.vm
        !call zgemm('N','N',2,Ntot,2,one,inv(1,1),2,vm(1,1),2,zero,tmp(1,1),2)
        !!tmp_p=um.tmp
        !call zgemm('N','N',Ntot,Ntot,2,one,um(1,1),Ntot,tmp(1,1),2,zero,tmp_p(1,1),Ntot)
        !!tmp_pp=tmp_p.inver
        !call zgemm('N','N',Ntot,Ntot,Ntot,one,tmp_p(1,1),Ntot,inver(1,1,k),Ntot,zero,tmp_pp(1,1),Ntot)
        !!inver-inver.tmp_pp
        !call zcopy(Ntot*Ntot,inver(1,1,k),1,tmp_p(1,1),1)
        !call zgemm('N','N',Ntot,Ntot,Ntot,(-1.d0,0.d0),tmp_p(1,1),Ntot,tmp_pp(1,1),Ntot,one,inver(1,1,k),Ntot)

        call update_inverse(i,j,k,nj(1,1,k),aux,um(1,1,k),vm(1,1,k),inver(1,1,k))
        
     end do
   else if(j.eq.Nsite) then
     !We do not need to update the inverse
   else
     write(*,*) "Something is wrong with the onsit j input:",j
     call mystop
   end if

 end if


!Store the propogation Ising spin for back propogation.
 if(back_pro) then
    back_store(j,i_back,i)=x
 end if


!get the new determinate
 if(dtype.EQ.'c') then
   do k=1,Ntot,1
      phi(j,k,i)=phi(j,k,i)*(explr_up+one)
      phi(j+Nsite,k,i)=phi(j+Nsite,k,i)*(explr_dn+one)
   end do
 else if(dtype.EQ.'d') then
   do k=1,Nspin(1),1
      phi(j,k,i)=phi(j,k,i)*(explr_up+one)
   end do
   do k=Nspin(1)+1,Ntot,1
      phi(j+Nsite,k,i)=phi(j+Nsite,k,i)*(explr_dn+one)
   end do
 end if


!test
! if(i.eq.1050) then
!    write(*,*) tot_imp(i),i,j
!    !write(*,*) inver(2,2,1) 
!    call get_test(i)
!    pause 
! end if
!end test

end subroutine single_U_phi


!The vm and um arrays are useful in calculating the green's function and updating the inverse.
!Here we set vm=(phi)_{j,1:Ntot}.Inverse, and the um=conjg(phiT(j,1:Ntot))
subroutine get_um_vm(i,j,um,vm,inver)
use param
use phiT_param
use model_param
use lattice_param
use project_param
use phi_param
implicit none
integer,intent(IN)::i,j
complex(kind=8),intent(IN)::inver(Ntot,Ntot,Dtot)
complex(kind=8),intent(OUT)::vm(2,Ntot,Dtot) 
complex(kind=8),intent(OUT)::um(Ntot,2,Dtot)
complex(kind=8)::vm_t(2,Ntot)
integer::k
if(dtype.EQ.'c') then
  !Get vm(2,Ntot)
   vm_t(1,1:Ntot)=phi(j,1:Ntot,i)
   vm_t(2,1:Ntot)=phi(j+Nsite,1:Ntot,i)
   do k=1,Dtot,1
      call zgemm('N','N',2,Ntot,Ntot,one,vm_t(1,1),2,inver(1,1,k),Ntot,zero,vm(1,1,k),2)
      um(1:Ntot,1,k)=Conjg(phiT(j,1:Ntot,k))
      um(1:Ntot,2,k)=Conjg(phiT(j+Nsite,1:Ntot,k))
   end do
else if(dtype.EQ.'d') then
  !Get vm(2,Ntot)
   vm_t(1,1:Nspin(1))=phi(j,1:Nspin(1),i)
   vm_t(2,(Nspin(1)+1):Ntot)=phi(j+Nsite,(Nspin(1)+1):Ntot,i)
   do k=1,Dtot,1
      call zgemv('T',Nspin(1),Nspin(1),one,inver(1,1,k),Ntot,vm_t(1,1),2,zero,vm(1,1,k),2)
      call zgemv('T',Nspin(2),Nspin(2),one,inver((Nspin(1)+1),(Nspin(1)+1),k),Ntot,vm_t(2,(Nspin(1)+1)),2, &
           zero,vm(2,(Nspin(1)+1),k),2)
      um(1:Nspin(1),1,k)=Conjg(phiT(j,1:Nspin(1),k))
      um((Nspin(1)+1):Ntot,2,k)=Conjg(phiT(j+Nsite,(Nspin(1)+1):Ntot,k))
   end do
end if
end subroutine get_um_vm


!This subroutine get the green's function matrix from nj=Transpose(vm.um)
subroutine get_green_matrix(um,vm,nj)
use param
use model_param
use phiT_param
implicit none
complex(kind=8),intent(IN)::vm(2,Ntot,Dtot)
complex(kind=8),intent(IN)::um(Ntot,2,Dtot)
complex(kind=8),intent(OUT)::nj(2,2,Dtot)
complex(kind=8),external::zdotu
integer::k
if(dtype.EQ.'c') then
  do k=1,Dtot,1
     !nj=(vm.um)^T=(um^T).(vm^T)
     call zgemm('T','T',2,2,Ntot,one,um(1,1,k),Ntot,vm(1,1,k),2,zero,nj(1,1,k),2)
  end do
else if(dtype.EQ.'d') then
  do k=1,Dtot,1
     nj(1,1,k)=zdotu(Nspin(1),vm(1,1,k),2,um(1,1,k),1)
     nj(2,2,k)=zdotu(Nspin(2),vm(2,(Nspin(1)+1),k),2,um((Nspin(1)+1),2,k),1)
     nj(1,2,k)=zero
     nj(2,1,k)=zero 
  end do
end if
end subroutine get_green_matrix


!Update the inverse matrix in single U to phi
subroutine update_inverse_3b(i,j,k,nj,aux,um,vm,inver)
use param
use rand_num
use phi_param
use model_param
use lattice_param
use phiT_param
use project_param
use method_param
implicit none
integer,intent(IN)::i,j,k,aux
complex(kind=8),intent(IN)::nj(2,2)
complex(kind=8),intent(IN)::um(Ntot,2),vm(2,Ntot)
complex(kind=8),intent(INOUT)::inver(Ntot,Ntot)
complex(kind=8)::inv(2,2)
complex(kind=8)::tmp_vm(2,Ntot),tmp_um(Ntot,2)
integer::m,n
complex(kind=8)::e_up_loc(2),e_dn_loc(2)

if(j.le.Nbravais)then
     e_up_loc(:)=explnd_up(:)
     e_dn_loc(:)=explnd_dn(:)
elseif(j.le.2*Nbravais)then
     e_up_loc(:)=explnx_up(:)
     e_dn_loc(:)=explnx_dn(:)
else
     e_up_loc(:)=explny_up(:)
     e_dn_loc(:)=explny_dn(:)
endif

if(dtype.EQ.'c') then

  !Get the inv(2,2)
  inv(1,1)=one+e_up_loc(aux)*nj(1,1)
  inv(1,2)=e_up_loc(aux)*nj(2,1)
  inv(2,1)=e_dn_loc(aux)*nj(1,2)
  inv(2,2)=one+e_dn_loc(aux)*nj(2,2)
  call inv2(inv(1,1))

  !inv=inv.Diag(expln_up(aux),expln_dn(aux))
  inv(1,1)=inv(1,1)*e_up_loc(aux)
  inv(2,1)=inv(2,1)*e_up_loc(aux)
  inv(1,2)=inv(1,2)*e_dn_loc(aux)
  inv(2,2)=inv(2,2)*e_dn_loc(aux)

  !tmp_um=inver.um
  call zgemm('N','N',Ntot,2,Ntot,one,inver(1,1),Ntot,um(1,1),Ntot,zero,tmp_um(1,1),Ntot)
  !tmp_vm=inv.vm (vm is defined before)
  call zgemm('N','N',2,Ntot,2,one,inv(1,1),2,vm(1,1),2,zero,tmp_vm(1,1),2)
  !inver-tmp_um.tmp_vm
  call zgemm('N','N',Ntot,Ntot,2,(-1.d0,0.d0),tmp_um(1,1),Ntot,tmp_vm(1,1),2,one,inver(1,1),Ntot)

else if(dtype.EQ.'d') then

  !Get the inv(2,2) with Diag(-expln_up(aux),-expln_dn(aux))
  inv(1,1)=(-1.d0,0.d0)*e_up_loc(aux)/(one+e_up_loc(aux)*nj(1,1))
  inv(2,2)=(-1.d0,0.d0)*e_dn_loc(aux)/(one+e_dn_loc(aux)*nj(2,2))


  !tmp_um=inver.um
  call zgemv('N',Nspin(1),Nspin(1),one,inver(1,1),Ntot,um(1,1),1,zero,tmp_um(1,1),1)
  call zgemv('N',Nspin(2),Nspin(2),one,inver((Nspin(1)+1),(Nspin(1)+1)),Ntot,um((Nspin(1)+1),2),1, &
           & zero,tmp_um((Nspin(1)+1),2),1)

  !inver-inv(i,i).tmp_um.vm
  !call
  !zgemm('N','N',Nspin(1),Nspin(1),1,inv(1,1),tmp_um(1,1),Ntot,vm(1,1),2,one,inver(1,1),Ntot)
  !call zgemm('N','N',Nspin(2),Nspin(2),1,inv(2,2),tmp_um(Nspin(1)+1,2),Ntot,  &
  !        & vm(2,Nspin(1)+1),2,one,inver(Nspin(1)+1,Nspin(1)+1),Ntot)
  call zgeru(Nspin(1),Nspin(1),inv(1,1),tmp_um(1,1),1,vm(1,1),2,inver(1,1),Ntot)
  call zgeru(Nspin(2),Nspin(2),inv(2,2),tmp_um(Nspin(1)+1,2),1,vm(2,Nspin(1)+1),2,inver(Nspin(1)+1,Nspin(1)+1),Ntot)




  !write(*,*) "here";call mystop
end if


end subroutine update_inverse_3b


!Update the inverse matrix in single U to phi
subroutine update_inverse(i,j,k,nj,aux,um,vm,inver)
use param
use rand_num
use phi_param
use model_param
use lattice_param
use phiT_param
use project_param
use method_param
implicit none
integer,intent(IN)::i,j,k,aux
complex(kind=8),intent(IN)::nj(2,2)
complex(kind=8),intent(IN)::um(Ntot,2),vm(2,Ntot)
complex(kind=8),intent(INOUT)::inver(Ntot,Ntot)
complex(kind=8)::inv(2,2)
complex(kind=8)::tmp_vm(2,Ntot),tmp_um(Ntot,2)
integer::m,n

if(dtype.EQ.'c') then

  !Get the inv(2,2)
  inv(1,1)=one+expln_up(aux)*nj(1,1)
  inv(1,2)=expln_up(aux)*nj(2,1)
  inv(2,1)=expln_dn(aux)*nj(1,2)
  inv(2,2)=one+expln_dn(aux)*nj(2,2)
  call inv2(inv(1,1))

  !inv=inv.Diag(expln_up(aux),expln_dn(aux))
  inv(1,1)=inv(1,1)*expln_up(aux)
  inv(2,1)=inv(2,1)*expln_up(aux)
  inv(1,2)=inv(1,2)*expln_dn(aux)
  inv(2,2)=inv(2,2)*expln_dn(aux)
  
  !tmp_um=inver.um
  call zgemm('N','N',Ntot,2,Ntot,one,inver(1,1),Ntot,um(1,1),Ntot,zero,tmp_um(1,1),Ntot)
  !tmp_vm=inv.vm (vm is defined before)
  call zgemm('N','N',2,Ntot,2,one,inv(1,1),2,vm(1,1),2,zero,tmp_vm(1,1),2)
  !inver-tmp_um.tmp_vm
  call zgemm('N','N',Ntot,Ntot,2,(-1.d0,0.d0),tmp_um(1,1),Ntot,tmp_vm(1,1),2,one,inver(1,1),Ntot)

else if(dtype.EQ.'d') then

  !Get the inv(2,2) with Diag(-expln_up(aux),-expln_dn(aux))
  inv(1,1)=(-1.d0,0.d0)*expln_up(aux)/(one+expln_up(aux)*nj(1,1))
  inv(2,2)=(-1.d0,0.d0)*expln_dn(aux)/(one+expln_dn(aux)*nj(2,2))


  !tmp_um=inver.um
  call zgemv('N',Nspin(1),Nspin(1),one,inver(1,1),Ntot,um(1,1),1,zero,tmp_um(1,1),1)
  call zgemv('N',Nspin(2),Nspin(2),one,inver((Nspin(1)+1),(Nspin(1)+1)),Ntot,um((Nspin(1)+1),2),1, &
           & zero,tmp_um((Nspin(1)+1),2),1)

  !inver-inv(i,i).tmp_um.vm
  !call zgemm('N','N',Nspin(1),Nspin(1),1,inv(1,1),tmp_um(1,1),Ntot,vm(1,1),2,one,inver(1,1),Ntot)
  !call zgemm('N','N',Nspin(2),Nspin(2),1,inv(2,2),tmp_um(Nspin(1)+1,2),Ntot,  &
  !        & vm(2,Nspin(1)+1),2,one,inver(Nspin(1)+1,Nspin(1)+1),Ntot)
  call zgeru(Nspin(1),Nspin(1),inv(1,1),tmp_um(1,1),1,vm(1,1),2,inver(1,1),Ntot)
  call zgeru(Nspin(2),Nspin(2),inv(2,2),tmp_um(Nspin(1)+1,2),1,vm(2,Nspin(1)+1),2,inver(Nspin(1)+1,Nspin(1)+1),Ntot)


  

  !write(*,*) "here";call mystop
end if

end subroutine update_inverse



!test
!subroutine get_test(k)
!use param
!use phiT_param
!use phi_param
!use model_param
!use mpi_serial_param
!use lattice_param
!implicit none
!integer,intent(IN)::k
!complex(kind=8)::imp,tot
!complex(kind=8)::inv(Ntot,Ntot)
!integer::m
!tot=zero
!do m=1,Dtot,1
!   !call deter_overlap(2*Nsite,Ntot,phiT(1,1,m),phi(1,1,k),inv(1,1))
!   !call inverse(inv(1,1),Ntot)
!   !if(m.eq.1) write(*,*) inv(2,2)
!   call deter_overlap_imp(2*Nsite,Ntot,phiT(1,1,m),phi(1,1,k),imp)
!   if(abs(imp).LT.1d-15) then
!      write(*,*) "Small overlap phiT with phi",k,m,rank
!   !  call mystop
!   end if
!   !impfunc(m,k)=imp
!   tot=tot+conjg(coe_multi(m))*imp
!end do
!   write(*,*) tot
!end subroutine get_test
!end test
