!This subroutine measure the quanity. It is combined with FW and BK measure.

!---------------------------
!The free projection measure
!---------------------------
subroutine freep_measure(i,j)
use param
use lattice_param
use meas_param
use adET_param
use mpi_serial_param
use mc_loop_param
use one_meas_param
implicit none
integer,intent(IN)::i,j  !i is the sample,j is the i_local
integer::sitei,sitej

 call measure()
 kin_l(i,j)=dble(K_one)
 v_l(i,j)=dble(V_one)
 e_l(i,j)=dble(E_one)
 var_l(i,j)=var_one
 nu_l(i,j)=dble(nu_one)
 nd_l(i,j)=dble(nd_one)
 do sitei=1,Nsite,1
    sisj_l(i,sitei,j)=dble(S_one(sitei))
 end do
 do sitei=1,2*Nsite,1
    cicj_l(i,sitei,j)=c_one(sitei)
 end do
 do sitei=1,2*Nsite,1
    do sitej=1,2*Nsite,1
       cicj_l_global(i,sitei,sitej,j)=c_one_global(sitei,sitej)
    enddo
 end do
 sig(j)=sig(j)+Sig_one
 absig(j)=absig(j)+AbSig_one

 if(j.lt.(max_ad*meastep)) then 
   ET=dble(E_one) !adjust the ET
   if(rank.eq.0) then
      write(*,*) "ADJUST ET FPMC:",ET
      write(*,*) ""
      write(*,*) ""
      write(*,*) ""
   end if
 end if
 !if(rank.eq.0) write(*,*) E_one
 !write(*,*) E_one;pause
end subroutine freep_measure


!--------------------------
!The cpmc and rcpmc measure
!--------------------------
subroutine add_measure(sample,i_release)
use param
use lattice_param
use meas_param
use one_meas_param
implicit none
integer,intent(IN)::sample,i_release
integer::sitei,sitej

  call measure()
  kin_l(sample,i_release)=kin_l(sample,i_release)+dble(K_one)
  v_l(sample,i_release)=v_l(sample,i_release)+dble(V_one)
  e_l(sample,i_release)=e_l(sample,i_release)+dble(E_one)
  var_l(sample,i_release)=var_l(sample,i_release)+var_one 
  nu_l(sample,i_release)=nu_l(sample,i_release)+dble(nu_one)
  nd_l(sample,i_release)=nd_l(sample,i_release)+dble(nd_one)
  do sitei=1,Nsite,1
     sisj_l(sample,sitei,i_release)=sisj_l(sample,sitei,i_release)+dble(S_one(sitei))
  end do
  do sitei=1,2*Nsite,1
     cicj_l(sample,sitei,i_release)=cicj_l(sample,sitei,i_release)+c_one(sitei)
  end do
  do sitei=1,2*Nsite,1
     do sitej=1,2*Nsite,1
        cicj_l_global(sample,sitei,sitej,i_release)=cicj_l_global(sample,sitei,sitej,i_release)+c_one_global(sitei,sitej)
     enddo
  end do
  sig(i_release)=sig(i_release)+Sig_one
  absig(i_release)=absig(i_release)+AbSig_one
end subroutine add_measure




!-----------------------------------------------------------------
!The measurement give out the K V E S C sig the result is only one
!measurement, it should work together with other add subroutines
!-----------------------------------------------------------------
subroutine measure()
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param

use io_module
implicit none
#ifdef MPI
include "mpif.h"
#endif

!denominator and the numerator of the measurement
complex(kind=8)::Kinm,Vinm,sisjm(Nsite),cicjm(2*Nsite),cicjm_global(2*Nsite,2*Nsite)
complex(kind=8)::htwom
complex(kind=8)::nu_m,nd_m
complex(kind=8)::denominator  !Determinate the sum of weight

!The local files: <phL|O_local|phr>
complex(kind=8)::Kinm_local(Dtot),Vinm_local(Dtot)
complex(kind=8)::htwom_local(Dtot)
complex(kind=8)::nu_local(Dtot),nd_local(Dtot)
complex(kind=8),allocatable::sisj_local(:,:)
complex(kind=8),allocatable::cicj_local(:,:)
complex(kind=8),allocatable::cicj_global(:,:,:)

!For the tmp wavefunction file
complex(kind=8),allocatable::phL(:,:,:)
complex(kind=8)::coe(Dtot)
complex(kind=8)::phR(2*Nsite,Ntot)
complex(kind=8),allocatable::ovp_local(:,:,:)
complex(kind=8)::imp_local(Dtot),tot_local
complex(kind=8),allocatable::Amat_local(:,:,:)

!For the weight meas
complex(kind=8)::w_meas,tot_meas

!The cycle parameter or tmp parameter
integer::i,j,k,sitei,sitej
complex(kind=8)::tmp


allocate(sisj_local(Nsite,Dtot))
allocate(cicj_local(2*Nsite,Dtot))
allocate(cicj_global(2*Nsite,2*Nsite,Dtot))
allocate(phL(2*Nsite,Ntot,Dtot))
allocate(ovp_local(Ntot,Ntot,Dtot))
allocate(Amat_local(2*Nsite,2*Nsite,Dtot))


Kinm=zero
Vinm=zero
Htwom=zero
nu_m=zero
nd_m=zero
sisjm=zero
cicjm=zero
cicjm_global=zero
denominator=zero
do i=1,Nwalkers,1
   !Get the phL,coe,phR and the ovp_local,imp_local,tot_local
   !also get the measure weight w_meas and tot_meas
   call get_meas_array(i,phL,coe,phR,ovp_local,imp_local,tot_local,w_meas,tot_meas)

   
   !To get the Amat(2*Nsite,2*Nsite,1:Dtot)=<phR(1:Dtot)|(ci^+)(cj)|phL>
   do k=1,Dtot,1
      !call cal_Amat_withovlpinv(2*Nsite,Ntot,phL(1,1,k),phR(1,1),ovp_local(1,1,k),Amat_local(1,1,k))
      call cal_Amat_withovlpinv_dc(phL(1,1,k),phR(1,1),ovp_local(1,1,k),Amat_local(1,1,k))
   end do



   !Kinectic energy

   !Get the Kinm_local(k)=<phL(k)|K|phR>/<phL(k)|phR>
   call get_klocal(Kinm_local,Amat_local)
   !Get the numerator Kinm=w_meas*rx_i*<phL|K|phR>/d_tmpi.
   !cpmc: d_tmpi=<phL|phR>
   !fpmc and rcpmc: d_tmpi=one
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,Kinm_local,Kinm)



   !Potential energy

   !Get the Vinm_local(k)=<phL(k)|V|phR>/<phL(k)|phR>
   call get_vlocal(Vinm_local,Amat_local)
   !Get the numerator Vinm=w_meas*rx_i*<phL|V|phR>/d_tmpi
   !cpmc: d_tmpi=<phL|phR>
   !fpmc and rcpmc: d_tmpi=one
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,Vinm_local,Vinm)



   !htwo energy

   !Get the htwom_local(k)=<phL(k)|H^2|phR>/<phL(k)|phR>
   !call get_htwolocal(htwom_local,Amat_local)
   !Get the numerator htwom=w_meas*rx_i*<phL|H^2|phR>/d_tmpi.
   !cpmc: d_tmpi=<phL|phR>
   !fpmc and rcpmc: d_tmpi=one
   !call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,htwom_local,htwom)




   !number of system
   !Get the nu_local(k)=<phL(k)|Nup|phR>/<phL(k)|phR>
   !        dn_local(k)=<phL(k)|Ndn|phR>/<phL(k)|phR>
   call get_nlocal(nu_local,nd_local,Amat_local)
   !Get the numerator nud_m=w_meas*rx_i*<phL|Nup,Ndn|phR>/d_tmpi
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,nu_local,nu_m)
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,nd_local,nd_m)


   !The spin correlation
   !Get the sizsjz_local(1:Nsite,k)=<phL(k)|S1Sj|phR>/<phL(k)|phR>
   call measure_sisj(Amat_local(1,1,1),sisj_local(1,1))
   do sitei=1,Nsite,1
      call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,sisj_local(sitei,1:Dtot),sisjm(sitei))
   enddo


   !The cicj correlation
   !Get the cicj_local(1:Nsite,k)=<phL(k)|c1cj|phR>/<phL(k)|phR>
   call measure_cicj(Amat_local(1,1,1),cicj_local(1,1))
   do sitei=1,2*Nsite,1
      call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,cicj_local(sitei,1:Dtot),cicjm(sitei))
   enddo
   
   !-------------------------------------
   !C^\dagger iCj global (cicj_local is induced by Transitional Symmetry)
   !-------------------------------------
   cicj_global=Amat_local;
   do sitei=1,2*Nsite,1
      do sitej=1,2*Nsite,1
         call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,cicj_global(sitei,sitej,1:Dtot),cicjm_global(sitei,sitej))
      enddo
   enddo


   call add_denominator(i,tot_meas,w_meas,denominator)  
end do

!Get the measurement value to physics_one=numerator/denominator
#ifdef MPI
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(denominator,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
denominator=tmp
call MPI_ALLREDUCE(Kinm,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
Kinm=tmp
call MPI_ALLREDUCE(Vinm,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
Vinm=tmp
call MPI_ALLREDUCE(htwom,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
htwom=tmp
call MPI_ALLREDUCE(nu_m,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
nu_m=tmp
call MPI_ALLREDUCE(nd_m,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
nd_m=tmp
do sitei=1,Nsite,1
   call MPI_ALLREDUCE(sisjm(sitei),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
   sisjm(sitei)=tmp
end do
do sitei=1,2*Nsite,1
   call MPI_ALLREDUCE(cicjm(sitei),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
   cicjm(sitei)=tmp
end do
do sitei=1,2*Nsite,1
   do sitej=1,2*Nsite,1
      call MPI_ALLREDUCE(cicjm_global(sitei,sitej),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
      cicjm_global(sitei,sitej)=tmp
   enddo
enddo
#endif


!if(abs(denominator).LT.1d-8) then
!  write(*,*) "denominator zero",denominator
!  call mystop
!end if


K_one=Kinm/denominator
V_one=Vinm/denominator
E_one=(Kinm+Vinm)/denominator
var_one=(htwom/denominator-E_one**2) !/(E_one**2)
nu_one=nu_m/denominator
nd_one=nd_m/denominator
do sitei=1,Nsite,1
   S_one(sitei)=sisjm(sitei)/denominator
end do
do sitei=1,2*Nsite,1
   c_one(sitei)=cicjm(sitei)/denominator
end do
do sitei=1,2*Nsite,1
   do sitej=1,2*Nsite,1
      c_one_global(sitei,sitej)=cicjm_global(sitei,sitej)/denominator
   enddo
end do
Sig_one=denominator
AbSig_one=abs(denominator)

deallocate(sisj_local,cicj_local,cicj_global,phL,ovp_local,Amat_local)

end subroutine measure








!------------------------------------------------------------
!We write the wave function after exp_mhalf_K, get the w_meas
!tot_meas, we also get the phL and coe, phR, the overlap 
!information between left and right wave functions.
!------------------------------------------------------------
subroutine get_meas_array(i,phL,coe,phR,ovp_local,imp_local,tot_local,w_meas,tot_meas)
use param
use lattice_param
use model_param
use project_param
use phiT_param
use phi_param
use method_param
use caldet_module
use mpi_serial_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(OUT)::phL(2*Nsite,Ntot,Dtot),coe(Dtot),phR(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::ovp_local(Ntot,Ntot,Dtot)
complex(kind=8),intent(OUT)::imp_local(Dtot),tot_local,w_meas,tot_meas
complex(kind=8)::ovp_tmp(Ntot,Ntot),tmp,phtmp(2*Nsite,Ntot)

integer::j,k

!test
! phR=zero
! phL=zero
!end test


!Get the phR
if(.not.back_pro) then
  !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi(1,1,i),2*Nsite,zero,phR,2*Nsite)
  call k_to_ph_dc(exp_halfK,phi(1,1,i),phR)
else if(back_pro) then
  !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi0(1,1,i),2*Nsite,zero,phR,2*Nsite)
  call k_to_ph_dc(exp_halfK,phi0(1,1,i),phR)
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if


!Get the phL, coe
if(.not.back_pro) then
  !call zcopy(2*Nsite*Ntot*Dtot,phiT(1,1,1),1,phL(1,1,1),1)
  call copy_wf_T_dc(phiT(1,1,1),phL(1,1,1))
  call zcopy(Dtot,coe_multi(1),1,coe(1),1)
else if(back_pro) then
  call back_prog_phL(i,phL,coe)
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if




!Get the ovp_local,imp_local,tot_local
tot_local=zero
do k=1,Dtot,1
   !call deter_overlap(2*Nsite,Ntot,phL(1,1,k),phR(1,1),ovp_local(1,1,k))
   call over_lap_dc(phL(1,1,k),phR(1,1),ovp_local(1,1,k))
   !call caldet_dc(ovp_local(1:Ntot,1:Ntot,k),imp_local(k))
   !call inverse_dc(ovp_local(1:Ntot,1:Ntot,k))
   call inverse_d_dc(ovp_local(1:Ntot,1:Ntot,k),imp_local(k))
   tot_local=tot_local+conjg(coe(k))*imp_local(k)
end do



!Get the tot_meas
if(.not.back_pro) then
  tot_meas=tot_local
else if(back_pro) then
  !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi(1,1,i),2*Nsite,zero,phtmp,2*Nsite)
  call k_to_ph_dc(exp_halfK,phi(1,1,i),phtmp)
  tot_meas=zero
  do k=1,Dtot,1
     !call deter_overlap(2*Nsite,Ntot,phiT(1,1,k),phtmp(1,1),ovp_tmp(1,1))
     call over_lap_dc(phiT(1,1,k),phtmp(1,1),ovp_tmp(1,1))
     !call caldet(Ntot,ovp_tmp(1:Ntot,1:Ntot),tmp)
     call caldet_dc(ovp_tmp(1:Ntot,1:Ntot),tmp)
     tot_meas=tot_meas+conjg(coe_multi(k))*tmp
  end do
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if


!Get the w_meas
if(crn.LT.0.d0) then
  w_meas=weight(i)*tot_meas/tot_imp(i)
else 
  w_meas=weight(i)*one
end if

end subroutine get_meas_array




!--------------------------------------------
!Get the kinm_local(k)=<phR|K|phL>/ <phR|phL>
!the kinm_local is should from 1~Dtot
!--------------------------------------------
subroutine get_klocal(Kinm_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(OUT)::Kinm_local(Dtot)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
integer::k,sitei,sitej
Kinm_local=zero
if(dtype.EQ.'c') then
   do k=1,Dtot,1
     do sitei=1,2*Nsite,1
        do sitej=1,2*Nsite,1
           Kinm_local(k)=Kinm_local(k)+Hzero(sitei,sitej)*Amat(sitei,sitej,k)
        enddo
     enddo
   end do
else if(dtype.EQ.'d') then
   do k=1,Dtot,1
     do sitei=1,Nsite,1
        do sitej=1,Nsite,1
           Kinm_local(k)=Kinm_local(k)+Hzero(sitei,sitej)*Amat(sitei,sitej,k)
        enddo
     enddo
     do sitei=Nsite+1,2*Nsite,1
        do sitej=Nsite+1,2*Nsite,1
           Kinm_local(k)=Kinm_local(k)+Hzero(sitei,sitej)*Amat(sitei,sitej,k)
        enddo
     enddo
   end do
end if
end subroutine get_klocal



!--------------------------------------------
!Get the Vinm_local(k)=<phR|V|phL>/ <phR|phL>
!the Vinm_local is should from 1~Dtot
!--------------------------------------------
subroutine get_vlocal(Vinm_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(OUT)::Vinm_local(Dtot)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
integer::k,sitei
Vinm_local=zero
if(dtype.EQ.'c') then
  do k=1,Dtot,1
     do sitei=1,Nsite,1
        Vinm_local(k)=Vinm_local(k)+Amat(sitei,sitei,k)*Amat(sitei+Nsite,sitei+Nsite,k) &
                   & -Amat(sitei,sitei+Nsite,k)*Amat(sitei+Nsite,sitei,k)
     end do
     Vinm_local(k)=dcmplx(onsitU)*Vinm_local(k)
  end do
else if(dtype.EQ.'d') then
  do k=1,Dtot,1
     do sitei=1,Nsite,1
        Vinm_local(k)=Vinm_local(k)+Amat(sitei,sitei,k)*Amat(sitei+Nsite,sitei+Nsite,k)
     end do
     Vinm_local(k)=dcmplx(onsitU)*Vinm_local(k)
  end do
end if
end subroutine get_vlocal


!-------------------------------------------
!measure the nu_local and nd_local from Amat
!-------------------------------------------
subroutine get_nlocal(nu_local,nd_local,Amat)
use param
use lattice_param
use phiT_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::nu_local(Dtot),nd_local(Dtot)
integer::i,j,k,m,n

nu_local=zero;nd_local=zero
do k=1,Dtot,1
   do i=1,Nsite,1
      nu_local(k)=nu_local(k)+Amat(i,i,k)
      nd_local(k)=nd_local(k)+Amat(i+Nsite,i+Nsite,k)
   end do
end do
end subroutine get_nlocal


!---------------------------------------------------------------
!We add the numberator to addm when measuring different quantity
!---------------------------------------------------------------
subroutine add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,add_local,addm)
use param
use phiT_param
use phi_param
use method_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(IN)::imp_local(Dtot),tot_local,coe(Dtot)
complex(kind=8),intent(IN)::w_meas,tot_meas
complex(kind=8),intent(IN)::add_local(Dtot)
complex(kind=8),intent(INOUT)::addm
complex(kind=8)::rat_h
integer::k


!rat_h=<phL|hat{O}|phR>
rat_h=zero
do k=1,Dtot,1
   rat_h=rat_h+conjg(coe(k))*imp_local(k)*add_local(k)
end do

!addm:
!rcpmc,fpmc: w_meas*rx_i*tot_meas*<phL|hat{O}|phR>/<phL|phR>
!cpmc: w_meas*rx_i*<phL|hat{O}|phR>/<phL|phR>
if(crn.GT.0.d0) then
  addm=addm+w_meas*rx(i)*tot_meas*rat_h/tot_local
else
  addm=addm+w_meas*rx(i)*rat_h/tot_local
end if
end subroutine add_numerator



!-------------------------------------------
!We add the denominator to addm in measuring
!-------------------------------------------
subroutine add_denominator(i,tot_meas,w_meas,addm)
use phi_param
use method_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(IN)::tot_meas,w_meas
complex(kind=8),intent(INOUT)::addm
if(crn.GT.0.d0) then
   addm=addm+w_meas*rx(i)*tot_meas
else
   addm=addm+w_meas*rx(i)
end if
end subroutine add_denominator
