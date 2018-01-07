module twist
implicit none
integer::Ntt=20
real(kind=8),allocatable::eavg(:)
real(kind=8),allocatable::eavgcp(:)
end module twist

!main program of code
program main
use mpi_serial_param
use timing_module
use rand_num
use phiT_param
use lattice_param
use twist
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer::i

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
 call readparam()



allocate(eavg(Ntt),eavgcp(Ntt))
open(unit=55,file='twist.dat',status='old')
do i=1,Ntt,1

call init_genrand()
read(55,'(2f15.8)') kbound(1),kbound(2)
if(rank.eq.0) write(*,*) "------------------------------"
if(rank.eq.0) write(*,*) "------------------------------"
if(rank.eq.0) write(*,*) "------------------------------"
if(rank.eq.0) write(*,*) i
if(rank.eq.0) write(*,*) kbound(1),kbound(2)
if(rank.eq.0) write(*,*) "------------------------------"
if(rank.eq.0) write(*,*) "------------------------------"
if(rank.eq.0) write(*,*) "------------------------------"

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
 call get_phiT()



! call test()



!-----------------------------------------------
!Print the method we want to use:CPMC,RCPMC,FPMC
!-----------------------------------------------
 call chose_method()

! call init_ET()
! call write_phi_rank()


call end_genrand()
call deallocatearray()
#ifdef MPI
call MPI_TYPE_FREE(htype,ierr)
call MPI_TYPE_FREE(phtype,ierr)
#endif

if(rank.eq.0)  call get_energy(i)
end do

close(55)


if(rank.eq.0) call twist_data()
!--------------------
!Prepare for the stop
!--------------------
! call clean
deallocate(eavg,eavgcp)
call EndTiming()
#ifdef MPI
call MPI_Finalize(ierr)
#endif
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
    call cpmc()
  else
    if(rank.eq.0) write(*,*) "When max_crn EQ 0, we are using cpmc, crn must be LT 0:",crn
    call mystop
  end if  
else if(max_crn.GT.0) then !RCPMC run
  if(crn.GT.0.d0) then
    if(rank.eq.0) write(*,*) "Running Rlease Constraint for the result."
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
implicit none
integer,intent(INOUT)::i_pop,i_GS
logical::ad_ET
integer::i_therm,max_therm

max_therm=Thermblock*blockstep
ad_ET=.false.;E_sum=0.d0;E_sum_old=0.d0;i_ad=0

!call K_phi(1)
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


subroutine get_energy(i)
use twist
use mc_loop_param
use io_module
use meas_param
implicit none
integer,intent(IN)::i
integer::j
real(kind=8)::pt,e_t,e_e,v_t,v_e,k_t,k_e,var_t,var_e,sgn
real(kind=8)::var_tmp
call get_filename()
call openUnit(EnergyName,10,'O')
var_tmp=1d100
do j=1,max_local,1
   read(10,'(10f15.8)') pt,k_t,k_e,v_t,v_e,e_t,e_e,var_t,var_e,sgn
   if(abs(var_t).LT.abs(var_tmp))  then
     eavg(i)=e_t
     var_tmp=var_t 
   end if
end do
eavgcp(i)=e_t
close(10)
end subroutine get_energy



subroutine twist_data
use twist
use mpi_serial_param
implicit none
real(kind=8)::e_t,e_e
real(kind=8)::e_tcp,e_ecp
call err_anal(eavg(1:Ntt),Ntt,e_t,e_e)
call err_anal(eavgcp(1:Ntt),Ntt,e_tcp,e_ecp)
  open(unit=10,file='twist-result.dat',status='replace')
      write(10,*) e_t,e_e,e_tcp,e_ecp
  close(10)
end subroutine twist_data
