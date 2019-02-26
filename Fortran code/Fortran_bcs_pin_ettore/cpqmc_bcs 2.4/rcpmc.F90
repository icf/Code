!program of rcpmc code
subroutine rcpmc()
use param
use mc_loop_param
use mpi_serial_param
use adET_param
use meas_param
use mpi_serial_param
use method_param
use meas_param
use model_param
implicit none
integer::sample !define the number of samples<Nsample.
integer::i_pop ! Use to define the step to do population control.
integer::i_GS  ! Use to define the step to do Modified GS.
integer::i_meas! Use to define the step to meas.
integer::i_local  !Use to define step to write each local physical quantity.

!allocate the arrays
 call allocate_meas()
 call allocate_one_meas()
 kin_l=0.d0;v_l=0.d0;e_l=0.d0;var_l=zero;nu_l=0.d0;nd_l=0.d0
 if(I_obdm.eq.1)obdm_l=zero
 if(ipinn.eq.0)then
  sisj_l=0.d0;ninj_l=0.d0;szsz_l=0.d0;sk_l=0.d0
!cacb_l=zero
 endif
 if(ipinn.eq.1)then
  nofr_l=zero
 endif
!cicj_l=zero;ck_l=0.d0;
 sig=zero;absig=zero
 if(Npair.gt.0)didj_l=zero
 if(Nbeta.gt.0)then
   if(I_onebf.eq.1)then
     cicj_t_l=zero
     cicjh_t_l=zero
   endif
!   rho_t_l=zero
   if(I_twob.eq.1)then
     nupnup_t_l=zero
     ndnnup_t_l=zero
   endif
   GreenP_t_l=zero
   GreenH_t_l=zero
 endif 

!for back propogation
 back_pro=.false.

!Inital the phi wave function, with weight and rx
 call get_phi()

!Thermalize the MC 
i_pop=0;i_GS=0
call therm(i_pop,i_GS)


do sample=1,Nsamples,1

   do i_local=1,max_local,1

      do i_meas=1,meastep,1

         if(rank.eq.0) then
           write(*,*) "sample:",sample
           write(*,*) "i_local:",i_local
           write(*,*) "i_meas:",i_meas
           write(*,*) ""
           write(*,*) ""
           write(*,*) ""
         end if

         call one_propagation(i_pop,i_GS)

      end do

      call rcpmc_meas(sample,i_pop,i_GS)

   end do
end do


!Get the mean value of different i_local
call mean_meas()

!data manipulate and write into the file
if(rank.eq.0) call data_cpmc_rcpmc()
end subroutine rcpmc


!-------------------------------------
!The subroutine for release constraint
!-------------------------------------
subroutine rcpmc_meas(sample,i_pop,i_GS)
use phi_param
use phiT_param
use lattice_param
use mc_loop_param
use model_param
use method_param
use mpi_serial_param
use backup_param
implicit none
integer,intent(IN)::sample
integer,intent(INOUT)::i_pop,i_GS
integer::i_release
integer::i,j,k


!backup:
call back_phi(i_pop,i_GS)


!prepare for the release
!if(rank.eq.0) write(*,*) "prepare for the release"
call Modified_GS();i_GS=0
crn=1.d0
do i=1,Nwalkers,1
   rx(i)=rx(i)/tot_imp(i)
end do



!release
!if(rank.eq.0) write(*,*) "start release"
if(fw_bk.EQ.'FW') call add_measure(sample,0)
do i_release=1,max_crn,1

   !if(rank.eq.0) write(*,*) i_release;pause
   
   call one_propagation(i_pop,i_GS)

   if(fw_bk.EQ.'FW') call add_measure(sample,i_release)

end do
if(fw_bk.EQ.'BK') then
  back_pro=.true.
  call back_propag(sample,i_pop,i_GS,0)
  back_pro=.false.
end if


!copy back
call recov_phi(i_pop,i_GS)
end subroutine rcpmc_meas
