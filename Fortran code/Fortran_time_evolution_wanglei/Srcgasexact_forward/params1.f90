!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: prec
! NOTICE :
! TYPE   : module
! PURPOSE: For defining the real kind value Rkind. 
! I/O    :
! COMMENT: This module is for RS and Param. Here since main program did not
!          use Rkind, so Rkind=8 should always be used, can not use 16.
!          otherwise there will be problem in transfering data to RS().
!========+=========+=========+=========+=========+=========+=========+=$
   module prec
   implicit none
   integer,parameter::Rkind=8
   end module prec
!========================================================================


!========+=========+=========+=========+=========+=========+=========+=$
 module hubbard_param
   use prec
   implicit none

  
   !real(kind=Rkind), parameter:: Pi=3.1415926535_Rkind     !Pi
   !complex(kind=Rkind), parameter:: Xi=(0.0_Rkind, 1.0_Rkind)    !imaginary unit!    
      
!-------------lattice size----------------------------------------------
   integer:: Nsite=6        !total site number =Nrank *Ncol.
   integer(kind=4):: Nsmax=400 !853776 !2944656 !11778624 !63504 ! !4900 ! !  !2944656  !  ! !400   ! !213444    !  !63504 ! the largest space 
                                           !largest subspace (2,2).
!--- model parameters -------------
   real(kind=Rkind):: Thop=1.0_RKind !hopping amplitude
   real(kind=Rkind):: U=1.0_RKind    !on-site U 
   real(kind=Rkind):: Mu=0.0_RKind     !chemical potential 
   real(kind=Rkind):: Vext=0.0_RKind     !external potential
   real(kind=Rkind):: stagger_h=0.0_RKind     !staggered magnetic field
   real(kind=Rkind),allocatable:: Tmn(:,:)              !hopping matrix.
!--- variables ---------------------
   integer(kind=4),allocatable::Table(:) !configuration(numbering)
   integer(kind=4), allocatable::invTable(:) !numbering(configuration)

   integer(kind=4),allocatable:: Nstate(:,:)        !containing the dimension D(nup, ndo)                                  
  end module hubbard_param                                            
!==============================================================



!========================================================================================
module lanc_param                          !parameters 
use prec
!----------parameter for the modified Lanczos method-----------------------
real (kind=Rkind),parameter:: Pre=1d-6  ! control of the convergence
integer,parameter :: Lanc_size=100   ! size of the Lanczos Matrix
integer, parameter:: Noflow=1
integer:: stepforcycle=1
end module lanc_param
!===================================================================================================

!==============================================================================
module flags_param                          !flags

!character(len=132):: Dir
!real(kind=8):: dt=1d-1
!integer:: TrotterOrder=2
!real(kind=8):: Maxtime=20d0
!integer:: Maxsteps=100
integer, parameter:: Max_threads=8
integer, parameter:: Iterprint=0              ! control print in lanczos subroutine 
integer(kind=4), parameter:: Condition_omp=5000
integer, parameter:: Schmit_flag=1
integer,parameter:: invTable_flag=1  ! 1 produce invTable, or look up in th e Table
end module flags_param
!=============================================================================

