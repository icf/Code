subroutine green_particles()
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param
use project_param
use method_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

integer::Qx,Qy

!denominator and the numerator of the measurement
complex(kind=8)::GreenPm_t(0:Nbeta)
complex(kind=8)::denominator  !Determinate the sum of weight

!The local files: <phL|O_local|phr>
complex(kind=8),allocatable::GreenP_t_local(:,:)

!For the tmp wavefunction file
complex(kind=8),allocatable::phL(:,:,:)
complex(kind=8)::coe(Dtot)
complex(kind=8)::phR(2*Nsite,Ntot),phRstart(2*Nsite,Ntot)
complex(kind=8)::phRbeta(2*Nsite,Ntot,0:Nbeta),MuBeta(2*Nsite,0:Nbeta)
complex(kind=8)::ovp(Ntot,Ntot)


complex(kind=8)::imp_local(Dtot,0:Nbeta),tot_local(0:Nbeta)

!For the weight meas
complex(kind=8)::w_meas,tot_meas

!The cycle parameter or tmp parameter
integer::icall,flag_zero
integer::i,j,k,sitei,sitej,aux,ib,jb,kb,i_beta,alpha,site,ip,ichoice
integer:: i_deb,j_deb,k_deb !DEBUG
real(kind=8)::x,temp
real(kind=8)::detD,normr,small
complex(kind=8)::tmp,norm,ph
complex(kind=8)::imp
complex(kind=8)::c00,proj0
complex(kind=8)::mu_temp(2*Nsite),mu_tilde(2*Nsite),dnpu(0:Nbeta)
complex(kind=8)::ui(2*Nsite),proj(Ntot),projr(Ntot),projl(Ntot)
complex(kind=8)::Rmat(Ntot,Ntot)
complex(kind=8),external::zdotc
data icall/0/
save icall

small=1.d-8

!allocate arrays for local estimator
allocate(GreenP_t_local(0:Nbeta,Dtot))
allocate(phL(2*Nsite,Ntot,Dtot))

!wave vector and orbital for the dynamical correlation function
if(icall.eq.0)then
  if(rank.eq.0)then
    open(10,file='orbital',status='old')
      read(10,*)ichoice
      read(10,*)Qx,Qy  
    close(10)
  endif
  call MPI_BCAST(ichoice,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Qx,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Qy,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  if(ichoice.eq.1)then !k space
    do i=1,Nbravais
      ph=2.d0*Pi*(dble(Qx*(coor(i,1)-1))/dble(Nl(1))+dble(Qy*(coor(i,2)-1))/dble(Nl(2)))
      norm=dcmplx(1.0/dsqrt(dble(Nbravais)),0)
      mu(i)=norm*exp(Xi*ph)
    enddo
    do i=Nbravais+1,2*Nsite
      mu(i)=zero
    enddo
    mu0(:)=mu(:)

  elseif(ichoice.eq.0)then !r space

    mu(:)=zero
    mu(Qy)=one
    mu0(:)=zero
    mu0(Qx)=one

  elseif(ichoice.eq.2)then  !read orbital from file
    if(rank.eq.0)then
      open(10,file='orbital_components',status='old')
!       write(115,*)'Open file with orbital '
!       flush(115)
       do i=1,2*Nsite
         read(10,*)mu(i)
!         write(115,*)mu(i)
       enddo
      close(10)      
    endif
!    flush(115)
!    write(115,*)'MPI_BCAST ',2*Nsite
!    flush(115)
    call MPI_BCAST(mu,2*Nsite,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
    mu0(:)=mu(:)
  endif
!  write(115,*)'I am ready '
!  flush(115)

  icall=1

endif

flag_zero=0
GreenPm_t=zero
denominator=zero

do i=1,Nwalkers,1   !loop over walkers

!Collect the information about phL and phR, to build the matrix element
  call get_phiR_phiL(i,phL,coe,phR)
  call modGS(phL(1,1,1),2*Nsite,Ntot,detD,Rmat) !Alert a GS: Phi = U V D, detD = 1/det(D)
  call modGS(phR,2*Nsite,Ntot,detD,Rmat)


!Backup for phR
  phRstart(:,:)=phR(:,:)

!zero imaginary time
  call copy_wf_dc(phR(1,1),phRbeta(1,1,0))

 !make it orthogonal to other orbitals
  do ip=1,Ntot
    ui(:)=phR(:,ip)
    proj(ip)=zdotc(2*Nsite,ui,1,mu,1)
  enddo
  do j=1,2*Nsite
    mu_tilde(j)=mu(j)
    do ip=1,Ntot
      mu_tilde(j)=mu_tilde(j)-proj(ip)*phR(j,ip)
    enddo
  enddo
  normr=dsqrt(dble(zdotc(2*Nsite,mu_tilde,1,mu_tilde,1)))
  norm=dcmplx(normr,0.d0)  
 
  if(normr.lt.small)then
    MuBeta(:,0)=zero
    dnpu(0)=zero
    flag_zero=1
  else 
    mu_tilde(:)=mu_tilde(:)/norm
    MuBeta(:,0)=mu_tilde(:)
    dnpu(0)=zdotc(2*Nsite,mu_tilde,1,mu,1)
  endif

  if(flag_zero.eq.1)then
    GreenP_t_local(:,:)=zero
    go to 1
  endif
 

!finite imaginary time
  do i_beta=1,Nbeta    !loop over imaginary time instants

!propagate phR
    call copy_wf_dc(phRbeta(1,1,i_beta-1),phR(1,1))
    call one_step_left_bare(phR,i_beta,i)
    call modGS(phR,2*Nsite,Ntot,detD,Rmat)
    call copy_wf_dc(phR(1,1),phRbeta(1,1,i_beta))

!propagate the orbital
    mu_temp(:)=MuBeta(:,i_beta-1)
    call propagate_orbital(mu_temp,i_beta,i)

!make it orthogonal to other orbitals
    do ip=1,Ntot
      ui(:)=phRbeta(:,ip,i_beta)
      proj(ip)=zdotc(2*Nsite,ui,1,mu_temp,1)
    enddo
    do j=1,2*Nsite
      mu_tilde(j)=mu_temp(j)
      do ip=1,Ntot
        mu_tilde(j)=mu_tilde(j)-proj(ip)*phRbeta(j,ip,i_beta)
      enddo
    enddo
    normr=dsqrt(dble(zdotc(2*Nsite,mu_tilde,1,mu_tilde,1)))
    norm=dcmplx(normr,0.d0)

    if(normr.lt.small)then
      MuBeta(:,i_beta)=zero
      dnpu(i_beta)=zero
      flag_zero=1
    else
      mu_tilde(:)=mu_tilde(:)/norm
      MuBeta(:,i_beta)=mu_tilde(:)
      dnpu(i_beta)=zdotc(2*Nsite,mu_tilde,1,mu_temp,1)*dnpu(i_beta-1)
    endif


  enddo !end of loop over imaginary time instants


  do i_beta=Nbeta,0,-1    !loop over imaginary time instants


!information for the weighted average
! imp_local(k,..) = <phL_k | phR(i_beta)>
! tot_local(..) = sum_k imp_local(k,..)
! tot_meas = <phiT | phi(at the end of the BP path) >
! w_meas = weight(walker) * tot_meas/tot_imp(walker) 
    call get_useful_params(i,phL,coe,phRstart,                      &
                           imp_local(1,i_beta),tot_local(i_beta),w_meas,tot_meas)


    do ip=1,Ntot
      ui(:)=phRbeta(:,ip,i_beta)
      projr(ip)=zdotc(2*Nsite,mu0,1,ui,1)
    enddo
    mu_temp(:)=MuBeta(:,i_beta)  
    proj0=zdotc(2*Nsite,mu0,1,mu_temp,1)

    do k=1,Dtot,1

      do ip=1,Ntot
        ui(:)=phL(:,ip,k)
        projl(ip)=zdotc(2*Nsite,ui,1,mu_temp,1)
      enddo

      call over_lap_dc(phL(1,1,k),phRbeta(1,1,i_beta),ovp)
      call inverse_d_dc(ovp,imp)

      c00=zero
      do ip=1,Ntot
        do j=1,Ntot
          c00=c00+projr(ip)*ovp(ip,j)*projl(j)
        enddo
      enddo 
       
      if(abs(dnpu(i_beta)).lt.small)then
        GreenP_t_local(i_beta,k)=zero
      else
        GreenP_t_local(i_beta,k)=dnpu(i_beta)*(proj0-c00)
      endif
      
    enddo

    if(i_beta.gt.0)call one_step_right(phL,i_beta,i,coe)

  enddo  !end of loop over imaginary time instant

1 continue

  do i_beta=0,Nbeta,1
      call add_numerator(i,imp_local(1,i_beta),tot_local(i_beta),                    &
             coe,w_meas,tot_meas,GreenP_t_local(i_beta,1:Dtot),GreenPm_t(i_beta))
  enddo

  call add_denominator(i,tot_meas,w_meas,denominator)

enddo !end of loop over walkers



!t0=MPI_Wtime()

#ifdef MPI
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(denominator,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
denominator=tmp
!do i_beta=0,Nbeta,1
!  do sitei=1,Nbravais,1
!    call MPI_ALLREDUCE(rhom_t(sitei,i_beta),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
!    rhom_t(sitei,i_beta)=tmp
!  enddo
!enddo
do i_beta=0,Nbeta,1
    call MPI_ALLREDUCE(GreenPm_t(i_beta),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
    GreenPm_t(i_beta)=tmp
enddo

#endif

!do i_beta=0,Nbeta,1
!  do sitei=1,Nbravais,1
!    RHO_one(sitei,i_beta)=rhom_t(sitei,i_beta)/denominator
!  enddo
!enddo

do i_beta=0,Nbeta,1
    GreenP_one(i_beta)=GreenPm_t(i_beta)/denominator
!    write(98,*)i_beta*dt,GreenP_one(i_beta)
enddo

deallocate(GreenP_t_local)


end subroutine green_particles




!---------------------------------------------------------------------
subroutine propagate_orbital(orb_mu,i_beta,i_walker)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
use mc_loop_param
implicit none
integer,intent(IN)::i_beta,i_walker
complex(kind=8),intent(INOUT)::orb_mu(2*Nsite)

integer::i_b
real(kind=8)::x
complex(kind=8)::mu_temp(2*Nsite)
complex(kind=8)::explr_up,explr_dn,lr_up,lr_dn
integer::j,aux

!write(*,*)2*Nsite
mu_temp(:)=orb_mu(:)
!write(*,*)mu_temp
call zgemm('N','N',2*Nsite,1,2*Nsite,one,exp_halfK,2*Nsite,mu_temp,2*Nsite,zero,orb_mu,2*Nsite)
!write(*,*)
!write(*,*)
!write(*,*)mu_temp
!orb_mu(:)=mu_temp(:)

i_b=i_beta
do j=1,Nsite,1

  x=back_store(j,i_beta,i_walker)

  if(i_b.LE.Nbk(1).AND.i_b.GE.1) then
     if(kcrn.eq.1) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.2) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.3) then
       lr_up=x*sqrt(dcmplx(dt*onsitU))
       lr_up=lr_up-(dt*onsitU/2.d0)

       lr_dn=-x*sqrt(dcmplx(dt*onsitU))
       lr_dn=lr_dn-(dt*onsitU/2.d0)

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.4) then
       lr_up=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_up=lr_up+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       lr_dn=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_dn=lr_dn+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else
       write(*,*) "Something is wrong with kcrn input:",kcrn
       call mystop
     end if
  else if(i_b.GT.Nbk(1).AND.i_b.LE.Nstps_fwd) then
     aux=NINT(x)

     if(Nbands.eq.1)then

       explr_up=expln_up(aux)
       explr_dn=expln_dn(aux)

     elseif(Nbands.eq.3)then

       if(j.le.Nbravais)then
         explr_up=explnd_up(aux)
         explr_dn=explnd_dn(aux)
       elseif(j.le.2*Nbravais)then
         explr_up=explnx_up(aux)
         explr_dn=explnx_dn(aux)
       else
         explr_up=explny_up(aux)
         explr_dn=explny_dn(aux)
       endif

     endif

  else
     write(*,*) "Something is wrong of i_b propagate_orbital ",i_b
     call mystop
  end if

  orb_mu(j)=orb_mu(j)            *(explr_up+one)
  orb_mu(j+Nsite)=orb_mu(j+Nsite)*(explr_dn+one)


!  write(*,*)'i_b = ',i_b
!  write(*,*)'Auxiliary field ',x,aux
!  write(*,*)j,orb_mu(j),explr_up+one,explr_dn+one
!  write(*,*)
!  write(*,*)
enddo
!stop'deb'

mu_temp(:)=orb_mu(:)
call zgemm('N','N',2*Nsite,1,2*Nsite,one,exp_halfK,2*Nsite,mu_temp,2*Nsite,zero,orb_mu,2*Nsite)
!orb_mu(:)=mu_temp(:)

end subroutine propagate_orbital
