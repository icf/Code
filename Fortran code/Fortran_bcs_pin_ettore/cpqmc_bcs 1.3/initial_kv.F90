subroutine inital_k_v()
use lattice_param
use mpi_serial_param
use project_param
use model_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

!allocated exp_K and exp_halfK and ng
call allocate_project_array()

!since different thread might get different result due to the degeneracy in
!Hzero, we only do rank.eq.0 and Bcast the result.
if(dtype.EQ.'c') then
  if(rank.eq.0) call initial_K(2*Nsite,Hzero,dt,exp_K,exp_halfK,exp_mhalfK)
else if(dtype.EQ.'d'.or.dtype.eq.'m') then
  if(rank.eq.0) call initial_K_d()
else
  write(*,*) "Something is wrong with dtype input:",dtype
  call mystop
end if


#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_BCAST(exp_K(1,1),1,htype,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(exp_halfK(1,1),1,htype,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(exp_mhalfK(1,1),1,htype,0,MPI_COMM_WORLD,IERR)
#endif


!get the gama,gamaf and ng(Nsite)
call initial_V() 

end subroutine inital_k_v



!-----------------------------------------------------------
!Input the h0,dt
!Output the exp(-dt*h0/2.d0),exp(-dt*h0) and exp(dt*h0/2.d0)
!-----------------------------------------------------------
subroutine initial_K(nl,h0,dt,exp_h0,exp_half_h0,exp_mhalf_h0)
use mpi_serial_param
implicit none
integer,intent(IN)::nl
complex(kind=8),intent(IN)::h0(nl,nl)
real(kind=8),intent(IN)::dt
complex(kind=8),intent(OUT)::exp_h0(nl,nl)
complex(kind=8),intent(OUT)::exp_half_h0(nl,nl)
complex(kind=8),intent(OUT)::exp_mhalf_h0(nl,nl)

complex(kind=8)::hu(nl,nl)
real(kind=8)::ev(nl)
integer::i,j,k


call check_Hermite_c(h0,nl)
call zcopy(nl*nl,h0,1,hu,1)
call eigen(Hu,nl,ev)

!if(rank.eq.0)then
! open(2,file='bands.info',status='unknown')
! do i=1,nl,1
!   write(2,*)i,ev(i)
! enddo
! close(2)
!endif

!-------------------------------------------------------
!Try to find the blacs code or use openmp to parallelize
!-------------------------------------------------------
exp_half_h0=dcmplx(0.d0,0.d0)
exp_mhalf_h0=dcmplx(0.d0,0.d0)
exp_h0=dcmplx(0.d0,0.d0)

do i=1,nl,1
   do j=1,nl,1
      do k=1,nl,1
         exp_half_h0(i,j)=exp_half_h0(i,j)+hu(i,k)*dcmplx(exp(-dt/2.d0*ev(k)))*conjg(hu(j,k))
         exp_mhalf_h0(i,j)=exp_mhalf_h0(i,j)+hu(i,k)*dcmplx(exp(dt/2.d0*ev(k)))*conjg(hu(j,k))
         exp_h0(i,j)=exp_h0(i,j)+hu(i,k)*dcmplx(exp(-dt*ev(k)))*conjg(hu(j,k))
      end do
   end do
end do


!DEBUG
!write(*,*)
!write(*,*)'exp_half_h0 '
!do j=1,nl,1
!  write(*,*)(exp_half_h0(j,i),i=1,nl)
!enddo
!write(*,*)
!write(*,*)

end subroutine initial_K


!--------------------------------------------------
!used to call initial_K for the decoupled condition
!--------------------------------------------------
subroutine initial_K_d
use lattice_param
use mpi_serial_param
use project_param
use model_param
implicit none
complex(kind=8)::h0(Nsite,Nsite)
complex(kind=8)::exp_h0(Nsite,Nsite)
complex(kind=8)::exp_half_h0(Nsite,Nsite)
complex(kind=8)::exp_mhalf_h0(Nsite,Nsite)
!complex(kind=8),allocatable::h0(:,:)
!complex(kind=8),allocatable::exp_h0(:,:)
!complex(kind=8),allocatable::exp_half_h0(:,:)
!complex(kind=8),allocatable::exp_mhalf_h0(:,:)
integer::i,j,k
integer        :: i_deb,j_deb  !DEBUG


!allocate(h0(Nsite,Nsite),exp_h0(Nsite,Nsite),exp_half_h0(Nsite,Nsite),exp_mhalf_h0(Nsite,Nsite))


h0(1:Nsite,1:Nsite)=Hzero(1:Nsite,1:Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),Hzero(1,1),2*Nsite,h0(1,1),Nsite)
call initial_K(Nsite,h0,dt,exp_h0,exp_half_h0,exp_mhalf_h0)
exp_K(1:Nsite,1:Nsite)=exp_h0(1:Nsite,1:Nsite)
exp_halfK(1:Nsite,1:Nsite)=exp_half_h0(1:Nsite,1:Nsite)
exp_mhalfK(1:Nsite,1:Nsite)=exp_mhalf_h0(1:Nsite,1:Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_h0(1,1),Nsite,exp_K(1,1),2*Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_half_h0(1,1),Nsite,exp_halfK(1,1),2*Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_mhalf_h0(1,1),Nsite,exp_mhalfK(1,1),2*Nsite)


h0(1:Nsite,1:Nsite)=Hzero((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),Hzero(Nsite+1,Nsite+1),2*Nsite,h0(1,1),Nsite)
call initial_K(Nsite,h0,dt,exp_h0,exp_half_h0,exp_mhalf_h0)
exp_K((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=exp_h0(1:Nsite,1:Nsite)
exp_halfK((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=exp_half_h0(1:Nsite,1:Nsite)
exp_mhalfK((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=exp_mhalf_h0(1:Nsite,1:Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_h0(1,1),Nsite,exp_K(Nsite+1,Nsite+1),2*Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_half_h0(1,1),Nsite,exp_halfK(Nsite+1,Nsite+1),2*Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_mhalf_h0(1,1),Nsite,exp_mhalfK(Nsite+1,Nsite+1),2*Nsite)

!DEBUG
!write(*,*)'Let us begin ...'
!write(*,*)
!write(*,*)'exp_halfK'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,exp_halfK(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,exp_halfK(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)

!stop 'DEBUG'



!deallocate(h0,exp_h0,exp_half_h0,exp_mhalf_h0)
end subroutine initial_K_d

!----------------------
!test with mathematica
!subroutine test()
!implicit none
!integer,parameter::nl=2
!complex(kind=8)::h0(2*nl,2*nl)
!real(kind=8)::dt=0.01d0
!complex(kind=8)::exp_h0(2*nl,2*nl),exp_half_h0(2*nl,2*nl),exp_mhalf_h0(2*nl,2*nl)
!h0(1,1)=dcmplx(1.d0,0.d0);h0(1,2)=dcmplx(0.d0,0.d0);h0(1,3)=dcmplx(1.d0,-1.d0);h0(1,4)=dcmplx(2.d0,-1.d0)
!h0(2,1)=dcmplx(0.d0);h0(2,2)=dcmplx(2.d0);h0(2,3)=dcmplx(0);h0(2,4)=dcmplx(3.d0,1.d0)
!h0(3,1)=dcmplx(1.d0,1.d0);h0(3,2)=dcmplx(0);h0(3,3)=dcmplx(3.d0);h0(3,4)=dcmplx(2.d0)
!h0(4,1)=dcmplx(2.d0,1.d0);h0(4,2)=dcmplx(3.d0,-1.d0);h0(4,3)=dcmplx(2.d0);h0(4,4)=dcmplx(4.d0)
!call initial_K(2*nl,h0,dt,exp_h0,exp_half_h0,exp_mhalf_h0)
!write(*,*) exp_mhalf_h0;call mystop
!end subroutine test
!end test
!-----------------------


!-----------------------------------------------------------------
!We initial the gama and gamaf which is used in H-S transformation
!Also initial the ng(2*Nsite) in the free projection
!-----------------------------------------------------------------
subroutine initial_V()
use param
use lattice_param
use model_param
use method_param
use project_param
implicit none
real(kind=8)::y,cosha
integer::i,j


if(Nbands.eq.1)then
!For the gama in CPMC
  if(dcp.eq.'S') then
    y=dt*onsitU/2.d0
    if(onsitU.ge.0.d0)then
      call get_gama(gama,y)
    else
      cosha = exp(0.5d0*dabs(onsitU)*dt)      
      gama=dcmplx(0.d0,dacos(1.d0/cosha))
    endif
  else if(dcp.eq.'C') then
    y=-1.d0*dt*onsitU/2.d0
    call get_gama(gama,y)
  else
    write(*,*) "Do not know what kind of decouple method."
    call mystop
  end if
!  call get_gama(gama,y)


!For expln in cpmc
  if(dcp.eq.'S') then
    !1 x=-1; 2 x=1
    expln_up(1)=exp(-1.d0*(dt*onsitU*0.5d0+gama))-one
    expln_up(2)=exp(-1.d0*(dt*onsitU*0.5d0-gama))-one

    expln_dn(1)=exp(-1.d0*(dt*onsitU*0.5d0-gama))-one
    expln_dn(2)=exp(-1.d0*(dt*onsitU*0.5d0+gama))-one
  else if(dcp.eq.'C') then
  !1 x=-1; 2 x=1
    expln_up(1)=exp(-1.d0*(dt*onsitU*0.5d0+gama))-one
    expln_up(2)=exp(-1.d0*(dt*onsitU*0.5d0-gama))-one

    expln_dn(1)=exp(-1.d0*(dt*onsitU*0.5d0+gama))-one
    expln_dn(2)=exp(-1.d0*(dt*onsitU*0.5d0-gama))-one
  else
    write(*,*) "Do not know what kind of decouple method."
    call mystop
  end if

elseif(Nbands.eq.3)then
  
!d orbitals
  if(dcp.eq.'S') then
    y=dt*onsitUd/2.d0
  else if(dcp.eq.'C') then
    y=-1.d0*dt*onsitUd/2.d0
  else
    write(*,*) "Do not know what kind of decouple method."
    call mystop
  endif
  call get_gama(gamad,y)

  if(dcp.eq.'S') then
    !1 x=-1; 2 x=1
    explnd_up(1)=exp(-1.d0*(dt*onsitUd*0.5d0+gamad))-one
    explnd_up(2)=exp(-1.d0*(dt*onsitUd*0.5d0-gamad))-one

    explnd_dn(1)=exp(-1.d0*(dt*onsitUd*0.5d0-gamad))-one
    explnd_dn(2)=exp(-1.d0*(dt*onsitUd*0.5d0+gamad))-one
  else if(dcp.eq.'C') then
  !1 x=-1; 2 x=1
    explnd_up(1)=exp(-1.d0*(dt*onsitUd*0.5d0+gamad))-one
    explnd_up(2)=exp(-1.d0*(dt*onsitUd*0.5d0-gamad))-one

    explnd_dn(1)=exp(-1.d0*(dt*onsitUd*0.5d0+gamad))-one
    explnd_dn(2)=exp(-1.d0*(dt*onsitUd*0.5d0-gamad))-one
  else
    write(*,*) "Do not know what kind of decouple method."
    call mystop
  end if

!px and py orbitals
  if(dcp.eq.'S') then
    y=dt*onsitUp/2.d0
  else if(dcp.eq.'C') then
    y=-1.d0*dt*onsitUp/2.d0
  else
    write(*,*) "Do not know what kind of decouple method."
    call mystop
  endif
  call get_gama(gamax,y)
  gamay=gamax

  if(dcp.eq.'S') then
    !1 x=-1; 2 x=1
    explnx_up(1)=exp(-1.d0*(dt*onsitUp*0.5d0+gamax))-one
    explnx_up(2)=exp(-1.d0*(dt*onsitUp*0.5d0-gamax))-one
    explny_up(1)=exp(-1.d0*(dt*onsitUp*0.5d0+gamay))-one
    explny_up(2)=exp(-1.d0*(dt*onsitUp*0.5d0-gamay))-one

    explnx_dn(1)=exp(-1.d0*(dt*onsitUp*0.5d0-gamax))-one
    explnx_dn(2)=exp(-1.d0*(dt*onsitUp*0.5d0+gamax))-one
    explny_dn(1)=exp(-1.d0*(dt*onsitUp*0.5d0-gamay))-one
    explny_dn(2)=exp(-1.d0*(dt*onsitUp*0.5d0+gamay))-one
  else if(dcp.eq.'C') then
  !1 x=-1; 2 x=1
    explnx_up(1)=exp(-1.d0*(dt*onsitUp*0.5d0+gamax))-one
    explnx_up(2)=exp(-1.d0*(dt*onsitUp*0.5d0-gamax))-one
    explny_up(1)=exp(-1.d0*(dt*onsitUp*0.5d0+gamay))-one
    explny_up(2)=exp(-1.d0*(dt*onsitUp*0.5d0-gamay))-one

    explnx_dn(1)=exp(-1.d0*(dt*onsitUp*0.5d0+gamax))-one
    explnx_dn(2)=exp(-1.d0*(dt*onsitUp*0.5d0-gamax))-one
    explny_up(1)=exp(-1.d0*(dt*onsitUp*0.5d0+gamay))-one
    explny_up(2)=exp(-1.d0*(dt*onsitUp*0.5d0-gamay))-one
  else
    write(*,*) "Do not know what kind of decouple method."
    call mystop
  end if


endif

!For the gamaf in Release or Free projection
if(kcrn.eq.1) then
  y=dt*onsitU/2.d0
else if(kcrn.eq.2) then
  y=-1.d0*dt*onsitU/2.d0
else if(kcrn.eq.3.or.kcrn.eq.4) then
  y=0.d0
else
  write(*,*) "Do not know what kind of decouple method in gama free projection."
  call mystop
end if
call get_gama(gamaf,y)


!For the ng background in Release or Free projection
if(kcrn.eq.3.or.kcrn.eq.1) then
  do i=1,Nsite,1
     !ng(i)=0.d0
     ng(i)=dble(Nspin(1)-Nspin(2))/dble(Nsite)
  end do
else if(kcrn.eq.4.or.kcrn.eq.2) then
  do i=1,Nsite,1
     ng(i)=dble(Ntot)/dble(Nsite)
  end do
end if

end subroutine initial_V



!---------------------------------
!solve the equation cosh(x)=exp(y)
!---------------------------------
subroutine get_gama(x,y)
implicit none
real(kind=8),intent(IN)::y
complex(kind=8),intent(OUT)::x
complex(kind=8)::ey,rdummy
if(abs(y)<1d-10) then
   x=dcmplx(0.d0,0.d0)
else
   ey=dcmplx(exp(y))
   x=log(ey-sqrt(ey**2-1.d0))
endif

!Check x
rdummy=log((exp(x)+exp(-x))/2d0)
if(abs(rdummy-y)<1d-10) then
  !print *, 'gamma checked!'
else
  print *, 'wrong x!'
  call mystop
endif
end subroutine get_gama
