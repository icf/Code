!main program of code
program main
use mpi_serial_param
use timing_module
use rand_num
use phiT_param
use lattice_param
use HF_param
use sc_loop_param
use model_param

use method_param
use mc_loop_param
use project_param
use phi_param
use adET_param

use io_module
implicit none
integer::i
integer::Neqblock_temp

!--------------------------------------------------
!Initial rank and Nsize, if mpi then inital the MPI
!--------------------------------------------------
 call start_code()

!-----------------------------------------------
!BeginTiming and EndTiming get the running time.
!-----------------------------------------------
 if(rank.eq.0) write(*,*) "Start to run the routine."
 call BeginTiming()
 call PrintTiming()


!----------------------------------------------
!initial parameter:random number and read param
!----------------------------------------------
 call init_genrand()
 call readparam()


!-----------------------------------------------------------
!Get the lattice information,Get Hzero matrix of the lattice
!-----------------------------------------------------------
 call set_lattice()

!--------------------------------------------------
!Set the mpi htype and phtype after set the lattice
!--------------------------------------------------
 call init_mpi_type()


!-------------------------------------
!Inital the projection part of K and V
!-------------------------------------
 call inital_k_v()


!-----------------------------------
!Get the try wave function, MD or SD
!-----------------------------------
 !Normal CPMC at first
 I_wavefun=1
 
 Neqblock_temp=Neqblock
 if(sc_initial_choose .EQ. 1) then
    Neqblock=1
    max_local=Neqblock*blockstep
 endif
 !----------------------------------
 !s.c. setting
 !----------------------------------
 sc_step_counter=0
 sc_loop_flag=0
 sc_ite_flag=0
 sc_loop_Na=Ntot
 !----------------------------------
 !HF setting
 !----------------------------------
 HF_num=1000
 HF_a=0.75
 allocate(cicj_sc_global(2*Nsite,2*Nsite))
 allocate(cicj_sc_global1(2*Nsite,2*Nsite))
 call get_phiT()

!-----------------------------------------------
!Print the method we want to use:CPMC,RCPMC,FPMC
!-----------------------------------------------

 call chose_method()
 
!step-qmc is choosed, icf, 1/2/2018
!SC_loop, icf, 1/7/2018
!-----------------------------------------------
 !BCS CPMC start at second step
 I_wavefun=2

 Neqblock=Neqblock_temp+1
 max_local=Neqblock*blockstep
 Dtot=1      !due to get_phiT coe_muti not Bcast
 sc_loop_flag=1;
 sc_loop_Na=Nspin(1)
 allocate(phi_sc(Nsite,Ntot,Dtot))
 allocate(BCS_sc(Nsite,Nsite))
 do i=1,sc_loop_num,1
     sc_step_counter=i
     call get_phiT()
     call chose_method()
     sc_ite_flag=0;
     call openUnit('step_counter',100,'R')
     write(100,*)'step_counter: ',i
     close(100)
 enddo

 call deallocate_sc_loop_param
!--------------------
!Prepare for the stop
!--------------------
 call clean
 
end program main




!-----------------------------------------------
!Print the method we want to use:CPMC,RCPMC,FPMC
!-----------------------------------------------
subroutine chose_method()
use mpi_serial_param
use method_param
implicit none
!old code
!if(crn.LT.0.d0) then
!  if(max_crn.NE.0) then
!    write(*,*) "Max_crn must be zero when we do not release CP."
!    call mystop
!  else
!    if(rank.eq.0) write(*,*) "Running CPQMC for the result."
!    call cpmc()
!  end if
!
!else if(crn.GT.0.d0) then
!  if(max_crn.LT.0) then
!    if(rank.eq.0) write(*,*) "Running Free Projectiong for the result."
!    call fpmc()
!  else if(max_crn.GT.0) then
!    if(rank.eq.0) write(*,*) "Running Rlease Constraint for the result."
!    if(rank.eq.0) write(*,*) "Release step is:",max_crn
!    crn=-1.d0 !Run CPQMC first
!    call rcpmc()
!  else
!    write(*,*) "max_crn can not be zero when crn GT 0.d0",max_crn
!    call mystop
!  end if
!end if


if(max_crn.LT.0) then !Step by step run
  if(rank.eq.0) write(*,*) "Running step by step QMC measurement."
  if(crn.GT.0.d0) then
    if(rank.eq.0) write(*,*) "This is a free projection run."
  else if(crn.LT.0.d0) then
    if(rank.eq.0) write(*,*) "This is a constraint run."
  else
    if(rank.eq.0) write(*,*) "Something is wrong with the crn input:",crn
    call mystop
  end if
  call step_qmc()
else if(max_crn.EQ.0) then !CPMC run
  if(crn.LT.0.d0) then
    if(rank.eq.0) write(*,*) "Running CPQMC for the result."
!    stop 'DEBUG'
    call cpmc()
  else
    if(rank.eq.0) write(*,*) "When max_crn EQ 0, we are using cpmc, crn must be LT 0:",crn
    call mystop
  end if  
else if(max_crn.GT.0) then !RCPMC run
  if(crn.GT.0.d0) then
    if(rank.eq.0) write(*,*) "Running Released Constraint for the result."
    if(rank.eq.0) write(*,*) "Release step is:",max_crn
    crn=-1.d0 !Run CPQMC first
    call rcpmc()  
  else
    if(rank.eq.0) write(*,*) "When max_crn GT 0, we are use rcpmc, crn must be GT 0",crn
    call mystop
  end if
end if


if(rank.eq.0) write(*,*) ""
if(rank.eq.0) write(*,*) ""
if(rank.eq.0) write(*,*) ""
end subroutine chose_method




!-------------------------------------------------------
!The therm subroutine in the QMC for both FPMC CPMC RPMC
!-------------------------------------------------------
subroutine therm(i_pop,i_GS)
use param
use mc_loop_param
use mpi_serial_param
use adET_param
use meas_param
use method_param
!DEBUG
use phi_param
implicit none
integer,intent(INOUT)::i_pop,i_GS
logical::ad_ET
integer::i_therm,max_therm

max_therm=Thermblock*blockstep
ad_ET=.false.;E_sum=0.d0;E_sum_old=0.d0;i_ad=0

!to avoid problems if PopContrlstep=1
if(PopContrlstep.eq.1)then 
  m_w=1.d0;m_w_old=m_w
endif

!call Modified_GS();call mystop

call K_phi_W(1)

!*******************
!We can also read ET
!******************* 
call init_ET()
!call mystop

if(rank.eq.0) write(*,*) "Start therm:"
do i_therm=1,max_therm,1
   if(rank.eq.0) write(*,*) i_therm

   i_pop=i_pop+1;i_GS=i_GS+1

   !propagate one step in MC
   call K_V_update()


   !Do periodically Modified GS
   !When free-projection it is done
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
     ad_ET=.true.
   end if
   call adjustET(ad_ET)


   !write phi to file
   !if(i_therm.eq.1) then
   !  call init_ET()
   !  call Modified_GS()
   !  call PopControl()
   !  call write_phi_rank()
   !  call mystop
   !end if
   
end do

!write phi_rank
!if(rank.eq.0) then
!  write(*,*) "write the phi of different rank into file"
!  write(*,*) "We'd better write the phi just after population control"
!  write(*,*) "The next run with the phi might not do pop contrl at the beginning."
!  write(*,*) ""
!  write(*,*) ""
!  write(*,*) ""
!end if
!call write_phi_rank()
!call init_ET()
!call mystop

end subroutine therm
