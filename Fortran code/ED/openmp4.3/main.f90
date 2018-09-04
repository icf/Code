!------------------------------------------------------------------------------
!Program: Lanczos Hubbard model
!AUTHOR:  Hao Shi
!VERSION: 6-July-2010
!NOTICE:
!TYPE:    Serial code (Can be paralleled by MKL in one node)
!PURPOSE: Get the ground and some exact of the hubbard model as below by mixed 
!         Lanczos Method.
!         H=kin*sum((ci_sigma^+) *(cj_sigma),i,j,sigma)+sum(onsitU*(ni_up)*(ni_down),i)
!         (1).lattice:
!           The code support any kind of lattice model and any kind of boundary
!         (2).Method:
!           The code use Mixed lanczos method, the user should Change the random
!           number seed to make sure the result is not changed (make sure the input wave
!           function is not orthogonality with the gound state.
!
!           Enforce orthogonolization is used in each lanczos step to make sure the result
!           is exact enough.The user can choose to Note it when the gound state is only needed.
!           The enforce language is in the lanczos.f90 subroutine lanmatrix
!           (call schmidt(k),surround by:make schmidt orthogonalize ). 
!           Comment it will make the routine run fast, but maybe not exact enough.
!           Only ground state: Comment all of it.(Lit can be larger,depend on user's need)
!           Some excited state: choose call schmidt2()
!           Want more precision: choose call schmidt()==>When NlanM is larger
!           than 100 schmidt is needed.
!
!           The lanczos matrix is set by hand in param--LanM. The user can change LanM to find
!           The best LanM to get the result fast. LanM=2, it will become the modified lanczos.
!
!         (3).The lattice number is 32 and the up down particle is 16 is the Limit condition of this
!             code. If larger, the hilbert space will be larger than 2**63-1(Double integer limit)
!
!         (4).Input
!             set: 1,2,3 (integer)
!                  set=1:user should write hop document by hand. In document, the first and the
!                  the second two column is the hoping term (i,j) in hubbard model. The third 
!                  column is the kinetic energy between i,j.So different lattice model have different
!                  hoping document. User can write a code to get there hoping document
!                  set=2:The period boundary condition is used. 
!                  set=3:The twist boundary condition is used.
!                  when set=1, there is no need to give Ndimension,Nl(1),Nl(2),Nl(3),kin
!                  when set=2,3 there is no need to give Nhop,Nlattice.The hoping is write by Ndimension
!                               Nl(1:3) and kin    
!             Nlattice: GE 1 (integer)
!                  Only when set=1, Nlattice is used. It denote the lattice number in the model.
!             Nhop: (integer)
!                  Only when set=1, Nhop is used. It denote the number of hoping terms in the model
!             Ndimension: GE 1 (integer)
!                  Only when set=2,3 Ndimension is used. It denote the dimension of the model
!             Nl(1:Ndimension):(integer)
!                  Only when set=2,3 Nl is used. Of course Nl(1:dimension) is useful. It denote different number
!                  in square lattice Nx,Ny,Nz. 
!             kbound(1:Ndimension) (integer)
!                  Only when set=3, kbound is used. It denote the twist number of the TBC. Kinetic*exp(i*kbound(l)*Pi)            
!             Nup,Ndown:0~Nlattice (integer)
!                 The number of spin up particles and spin down particles
!             Kin: (Double Complex)
!                 Kinetic enegry of the Hamiltionian. Notice: in the code, the Hamiltonian do not has a minus simbol.
!                 usually Kin=(-1.d0,0.d0)
!             OnsitU:(Double Real) 
!                 Onsit energy of the Hamiltonian.
!             LanM:(integer)
!                 It denote the size of the lanczos matrix.When LanM=2 the code is modified lanczos.If the Enforce 
!                 orthogonolization is used, the LanM large will lead more cputime. LanM larger, the converge will
!                 be faster, and the memory cost larger. User can chose different LanM by experience.
!             LanMnin,LanMmax,Lanadd,Nlanmatrix_set:(integer)
!                 Sometime we find that the eigenvalue can not be converge to
!                 real state,so add the four parameters.If Nlanmatrix is larger
!                 than Nlanmatrix_set, Nlanmatrix is set back to 1.And LanM is
!                 added by Lanadd. IF LanM is larger than LanMmax, LanM is set
!                 form [LanMmin,LanMmax-1].Usually, LanMmax should not be larger
!                 than 100, or the whole orthogonalize is needed.Initial LanM
!                 should be from [LanMmin,LanMmax-1].
!             Lit:(Double Real)
!                 The little number decide to exit the cycle. When Norm(H/phi>-<phi/H/phi>)<lit, we think find a eigenstate.
!                 Lit should not be large, or it will not good for excited state. Lit should not little than 1.d-14,
!                 Where numerical error is large.If only the gound state is needed, lit can be larger---depend on the precision
!                 the user needed.
!             seed:(integer)
!                 The seed for random number generator. Seed usually is odd number.(RNG is 19937)
!             Nexcit:(integer)
!                 The number of excited state to be calculated.
!             wstate:1,0(integer)
!                 wstate=1,write the eigenstate into the disk
!                 wstate=0,do not write.
!Email:   Boruoshihao@gmail.com
!COMMENT:

!-------------------------------------------------------------------------------------------
!The main program of the routine************************************************************
!-------------------------------------------------------------------------------------------
program main
use param
use timing_module
use rand_num
implicit none
integer::mi,mj,mk

!-----------------------------------------------
!BeginTiming and EndTiming get the running time.
!-----------------------------------------------
 write(*,*) "Start to run the routine."
 call BeginTiming()
 call PrintTiming()



!---------------------------------------
!set the param used in the hubbard model
!---------------------------------------
 Pi=2.d0*ASIN(1.d0)
 call readparam()
 !Set this first in the subroutine
 if(set.NE.1) then
   Nlattice=1
   do mi=1,Dimen,1
    Nlattice=Nlattice*Nl(mi)
   end do
   Nhop=Nlattice*2*Dimen
 end if

 !check the Nup and Ndown with Nlattice
 if(Nup.GT.Nlattice.OR.Ndn.GT.Nlattice) then
   write(*,*) "Something is wrong with the Nup and Ndown-->one of them larger than Nlattice"
   call mystop
 end if


!-------------------------------
!give the rndm seed of the sprng
!-------------------------------
 call init_genrand()


!add SS
 two_SS=Nup+Ndn
 !two_rSS is read from initial_end
 !two_rSS=ABS(Nup-Ndn)

! write(*,*) set,Nlattice,Nhop
! write(*,*) Dimen,Nl(1),Nl(2),Nl(3)
! write(*,*) kbound(1),kbound(2),kbound(3)
! write(*,*) Nup,Ndown,kin,onsitU;pause

! test
!  do
!  call mintomax()
!  write(*,*) LanM;pause
!  end do
! end test



!-----------------------------------------------
!set the binomail number table and hilbert space
!-----------------------------------------------
 !We cannot allocate arrays before get Nhilbert
 allocate(Binomial(0:Nlattice,0:Nlattice))
 call set_binomial() 

 Nhup=Binomial(Nup,Nlattice);Nhdn=Binomial(Ndn,Nlattice)
 Nhilbert=Binomial(Nup,Nlattice)*Binomial(Ndn,Nlattice)
 if(Nexcit.GT.Nhilbert) then
   write(*,*) "Nexcit should no larger than Nhilbert."
   call mystop
 end if

 write(*,*) "Nhilbert=",Nhilbert
 write(*,*) "set the Nhilbert and the Basis of the hilbert space------done"
 call PrintTiming()

!---------------------------------------------
!Allocate all the arrays after we get Nhilbert
!---------------------------------------------
 call allocatearray()


!-----------------------------------------------------------------
!set Nhop and hopt(Nhop) sit(Nhop,2) in three different conditions
!-----------------------------------------------------------------
 call set_lattice_hop()
 write(*,*) "get the hop information-------done"
 call PrintTiming()

 call set_Tmatrix()


!---------------------------
!Set the matrix of up and dn
!---------------------------
 call get_Hm()
 write(*,*) "get the H and T information-------done"
 call PrintTiming()


!-----------------------------------------
!use lanczos method to get the gound state
!-----------------------------------------
 call lanczos()

 call mystop
end program main
