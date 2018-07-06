!program of cpmc code
subroutine cpmc()
use param
use mc_loop_param
use mpi_serial_param
use adET_param
use meas_param
use mpi_serial_param
use method_param
use meas_param
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
 sisj_l=0.d0;sk_l=0.d0;cicj_l=zero;ck_l=0.d0
 sig=zero;absig=zero

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


      if(fw_bk.eq.'FW') then
        call add_measure(sample,0)
      else 
        back_pro=.true.
        call back_phi(i_pop,i_GS)
        call back_propag(sample,i_pop,i_GS,0)
        call recov_phi(i_pop,i_GS)
        back_pro=.false.
      end if

   end do
end do

!Get the mean value of different i_local
call mean_meas()

if(rank.eq.0) write(*,*) e_l

!data manipulate and write into the file
if(rank.eq.0) call data_cpmc_rcpmc()
end subroutine cpmc
