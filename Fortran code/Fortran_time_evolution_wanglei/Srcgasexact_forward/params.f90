module system_parameters
implicit none

complex(kind=8),parameter:: Xi=cmplx(0d0,1d0)

integer:: Nwalkers=1

!integer:: Nsite=2
integer:: nup=1
integer:: ndo=1
integer:: numupdo(2)

real(kind=8):: dt=1d-2
real(kind=8):: Maxshift=1d0
real(kind=8):: Upperbound=2d0
real(kind=8):: Lowerbound=0.5d0
integer::      Maxsteps=100
integer::      Totalsteps=100
integer::      PopContrlsteps=10
character(len=132):: Dir='../Dat/'
integer::      Nsamples=10
!logical::      ImpSamp=.true.
integer ::      Shift=0
logical ::      Phaseless=.true.
logical ::      Local_energy=.false.

!integer::      StepforMeasure=10


!real(kind=8):: Thop=1d0
!real(kind=8):: U=0d0

complex(kind=8):: gamma=0d0
integer:: TrotterOrder=2
integer:: HS = 1
integer:: IniState=1


!complex(kind=8),allocatable:: phiT(:,:,:)  ! Nsite,num,spin
complex(kind=8),allocatable:: phi_time(:)  ! exact wf
!complex(kind=8),allocatable:: phiexact(:)  ! exact wf
complex(kind=8),allocatable:: phi_ini(:,:,:)  ! Nsite,num,spin ! initial state
complex(kind=8),allocatable:: phi(:,:,:,:) ! Nsite,num,spin,walker

!complex(kind=8),allocatable:: phiHF(:,:,:,:)  ! Nsite,num,spin
!complex(kind=8),allocatable:: phiHF_persite(:,:,:,:,:)  ! Nsite,num,spin,step

!real(kind=8), allocatable:: occT(:,:) ! HF occupation
!real(kind=8), allocatable:: occHF(:,:,:) ! HF occupation

real(kind=8),allocatable:: weight(:) !
!real(kind=8),allocatable:: wscaling(:) !
real(kind=8),allocatable:: wphase(:) !

complex(kind=8),allocatable:: impfunc0(:),impfunc1(:)!

!complex(kind=8),allocatable:: Bmats(:,:,:,:,:) !Nsite,Nsite,2,timestep,walker

real(kind=8), allocatable:: Kmat(:,:)    ! Nsite,Nsite
complex(kind=8),allocatable:: expK(:,:)  ! Nsite,Nsite

real(kind=8),allocatable:: Ham(:,:)
complex(kind=8), allocatable:: ExpHam(:,:)
real(kind=8), allocatable:: EigenHam(:)

end module
