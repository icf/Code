!--------------------
!program of fpmc code
!--------------------
subroutine step_qmc()
use param
use mc_loop_param
use mpi_serial_param
use adET_param
use meas_param
use mpi_serial_param
use method_param
use phi_param
use lattice_param
implicit none
integer::sample !define the number of samples<Nsample.
integer::i_pop ! Use to define the step to do population control.
integer::i_GS  ! Use to define the step to do Modified GS.
integer::i_meas! Use to define the step to meas.
integer::i_local  !Use to define step to write each local physical quantity.


!Measure the QMC code
call allocate_meas()
!call allocate_meas_step()
call allocate_one_meas()
sig=zero
absig=zero
back_pro=.false. !step the back to be false at the beginning


do sample=1,Nsamples,1
   !Inital the phi wave function, with weight
   call deallocate_phi()
   call get_phi()


   !Therm the QMC code
   i_pop=0;i_GS=0
   call therm(i_pop,i_GS)

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
        call freep_measure(sample,i_local)
      else

        !Set the back walk step
        Nstps_fwd=Thermblock*blockstep+i_local*meastep
        if(allocated(back_store)) deallocate(back_store)
        allocate(back_store(Nsite,Nstps_fwd,Nwalkers))
        if(crn.GT.0.d0) then
          Nbk(1)=Nstps_fwd
          Nbk(2)=0
        else if(crn.LT.0.d0) then
          Nbk(1)=0
          Nbk(2)=Nstps_fwd
        end if


        back_pro=.true.
        call back_phi(i_pop,i_GS)
        call back_propag(sample,i_pop,i_GS,i_local)
        call recov_phi(i_pop,i_GS)
        back_pro=.false.
      end if

   end do
end do

!manipulate the data and write into file
if(rank.eq.0) call data_step_mc()
end subroutine step_qmc
