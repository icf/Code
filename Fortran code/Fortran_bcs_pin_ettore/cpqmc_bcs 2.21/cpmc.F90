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
use model_param

use sc_loop_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer::sample !define the number of samples<Nsample.
integer::i_pop ! Use to define the step to do population control.
integer::i_GS  ! Use to define the step to do Modified GS.
integer::i_meas! Use to define the step to meas.
integer::i_local  !Use to define step to write each local physical quantity.
!real(kind=8)::t0,t1,t2,t3
 
!allocate the arrays
if(sc_loop_flag .EQ. 0)then
   call allocate_meas()
   call allocate_one_meas()  
endif 

!initialization
call set_beginning_zero()

!for back propogation
back_pro=.false.

!Inital the phi wave function, with weight and rx
call get_phi()

!Thermalize the MC 
i_pop=0;i_GS=0

call therm(i_pop,i_GS)


do sample=1,Nsamples,1

   do i_local=1,max_local,1  !max_local= Neqblock*blockstep

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

        if(Nbk(1).gt.0)call back_phi(i_pop,i_GS)  !if the constraint is released 
        call back_propag(sample,i_pop,i_GS,0)
        if(Nbk(1).gt.0)call recov_phi(i_pop,i_GS)  !if the constraint is released

        back_pro=.false.
      end if

   end do
   call print_observables(sample)
end do

!Get the mean value of different i_local
call mean_meas()

if(rank.eq.0) write(*,*) e_l

!data manipulate and write into the file
if(rank.eq.0) call data_cpmc_rcpmc()

end subroutine cpmc
