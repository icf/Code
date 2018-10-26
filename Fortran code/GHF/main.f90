program main
use mpi_serial_param
use timing_module
use rand_num
use lattice_param
implicit none
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


!Get the xhf calculation
 call xhf()


!--------------------
!Prepare for the stop
!--------------------
 call clean

end program main
