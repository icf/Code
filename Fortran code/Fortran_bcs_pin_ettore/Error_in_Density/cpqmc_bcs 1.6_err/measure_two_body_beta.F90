subroutine measure_density_beta()
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

!sites for the two-body correlation function
integer::Ri,Rj  

!denominator and the numerator of the measurement
complex(kind=8)::rhom_t(2*Nsite,0:Nbeta)
complex(kind=8)::denominator  !Determinate the sum of weight

!The local files: <phL|O_local|phr>
complex(kind=8),allocatable::rho_t_local(:,:,:)

!For the tmp wavefunction file
complex(kind=8),allocatable::phL(:,:,:)
complex(kind=8)::coe(Dtot)
complex(kind=8)::phR(2*Nsite,Ntot),phR_prime(2*Nsite,Ntot),phR_sec(2*Nsite,Ntot),phRstart(2*Nsite,Ntot)
complex(kind=8)::phRbeta(2*Nsite,Ntot,0:Nbeta),phRbeta_prime(2*Nsite,Ntot,0:Nbeta),phRbeta_sec(2*Nsite,Ntot,0:Nbeta)
complex(kind=8)::ovp(Ntot,Ntot),ovp_prime(Ntot,Ntot),ovp_sec(Ntot,Ntot)
complex(kind=8)::wfa_n(Ntot,2*Nsite),wfb_n(Ntot,2*Nsite),wfc_n(Ntot,2*Nsite)


complex(kind=8)::imp_local(Dtot,0:Nbeta),tot_local(0:Nbeta)

!For the weight meas
complex(kind=8)::w_meas,tot_meas

!The cycle parameter or tmp parameter
integer::i,j,k,sitei,sitej,aux,ib,jb,kb,i_beta,ip,alpha,site
integer:: i_deb,j_deb,k_deb !DEBUG
real(kind=8)::x,temp
real(kind=8)::detD,detD_prime,detD_sec
complex(kind=8)::explr_up,explr_dn,tmp
complex(kind=8)::factor,zeta_prime,zeta_sec,imp,imp_prime,imp_sec
complex(kind=8)::part_a_up,part_a_dn,part_b_up,part_b_dn,part_c_dn,part_c_up
complex(kind=8)::ratio_prime(0:Nbeta),ratio_sec(0:Nbeta)
complex(kind=8)::Rmat(Ntot,Ntot)



allocate(rho_t_local(2*Nsite,0:Nbeta,Dtot))
allocate(phL(2*Nsite,Ntot,Dtot))

!reference point for the dynamical correlation function
Ri=1


rhom_t=zero
denominator=zero

do i=1,Nwalkers,1   !loop over walkers

!Collect the information about phL and phR, to build the matrix element
  call get_phiR_phiL(i,phL,coe,phR)
  call modGS(phL(1,1,1),2*Nsite,Ntot,detD,Rmat)
!Alert a GS: Phi = U V D, detD = 1/det(D)  
  call modGS(phR,2*Nsite,Ntot,detD,Rmat)

!Backup for phR
  phRstart(:,:)=phR(:,:)

!exp(n_{Ri,up}) applied to phR
  phR_prime(:,:)=phR(:,:)
  do ip=1,Ntot
    phR_prime(Ri,ip)=phR_prime(Ri,ip)*Nepero
  enddo

!exp(n_{Ri,dn}) applied to phR
  phR_sec(:,:)=phR(:,:)
  do ip=1,Ntot
    phR_sec(Ri+Nsite,ip)=phR_sec(Ri+Nsite,ip)*Nepero
  enddo

!DEBUG--------------------------------------------------
  write(*,*)'Nsite, Nbravais',Nsite, Nbravais
  write(*,*)
  write(*,*)'PhiR     '
  write(*,*)
  do i_deb=1,2*Nsite
    write(*,*)(phR(i_deb,j_deb),j_deb=1,Ntot)
  enddo
  write(*,*)
  write(*,*)'PhiR_prime     '
  write(*,*)
  do i_deb=1,2*Nsite
    write(*,*)(phR_prime(i_deb,j_deb),j_deb=1,Ntot)
  enddo
  write(*,*)
  write(*,*)'PhiR_sec     '
  write(*,*)
  do i_deb=1,2*Nsite
    write(*,*)(phR_sec(i_deb,j_deb),j_deb=1,Ntot)
  enddo
!DEBUG--------------------------------------------------

!zero imaginary time
  ratio_prime(0)=one !dcmplx(detD_prime,0.d0)
  ratio_sec(0)=one! dcmplx(detD_sec,0.d0)
  call copy_wf_dc(phR(1,1),phRbeta(1,1,0))
  call copy_wf_dc(phR_prime(1,1),phRbeta_prime(1,1,0))
  call copy_wf_dc(phR_sec(1,1),phRbeta_sec(1,1,0))

!finite imaginary time
  do i_beta=1,Nbeta    !loop over imaginary time instants

!propagate phR
    call copy_wf_dc(phRbeta(1,1,i_beta-1),phR(1,1))
    call one_step_left_new(phR,i_beta,detD,i)
    call copy_wf_dc(phR(1,1),phRbeta(1,1,i_beta))

!propagate phR_prime
    call copy_wf_dc(phRbeta_prime(1,1,i_beta-1),phR_prime(1,1))
    call one_step_left_new(phR_prime,i_beta,detD_prime,i)
    call copy_wf_dc(phR_prime(1,1),phRbeta_prime(1,1,i_beta))
    factor=dcmplx(detD/detD_prime,0.d0)
    ratio_prime(i_beta)=ratio_prime(i_beta-1)*factor

!propagate phR_sec
    call copy_wf_dc(phRbeta_sec(1,1,i_beta-1),phR_sec(1,1))
    call one_step_left_new(phR_sec,i_beta,detD_sec,i)
    call copy_wf_dc(phR_sec(1,1),phRbeta_sec(1,1,i_beta))
    factor=dcmplx(detD/detD_sec,0.d0)
    ratio_sec(i_beta)=ratio_sec(i_beta-1)*factor


!DEBUG--------------------------------------------------
    write(*,*)'********** Time instant  ',i_beta
    write(*,*)
    write(*,*)'PhiR propagated     '
    write(*,*)
    do i_deb=1,2*Nsite
      write(*,*)(phR(i_deb,j_deb),j_deb=1,Ntot)
    enddo
    write(*,*)
    write(*,*)'PhiR_prime propagated     '
    write(*,*)
    do i_deb=1,2*Nsite
      write(*,*)(phR_prime(i_deb,j_deb),j_deb=1,Ntot)
    enddo
    write(*,*)
    write(*,*)'PhiR_sec propagated     '
    write(*,*)
    do i_deb=1,2*Nsite
      write(*,*)(phR_sec(i_deb,j_deb),j_deb=1,Ntot)
    enddo
!DEBUG--------------------------------------------------
   

  enddo !end of loop over imaginary time instants



  do i_beta=Nbeta,0,-1    !loop over imaginary time instants


!information for the weighted average
! imp_local(k,..) = <phL_k | phR>
! tot_local(..) = sum_k imp_local(k,..)
! tot_meas = <phiT | phi(at the end of the BP path) >
! w_meas = weight(walker) * tot_meas/tot_imp(walker) 
    call get_useful_params(i,phL,coe,phRstart,                      &
                         imp_local(1,i_beta),tot_local(i_beta),w_meas,tot_meas)


!DEBUG--------------------------------------------------
    write(*,*)
    write(*,*)'************************************'
    write(*,*)'Imaginary time instant ',i_beta
    write(*,*)
    write(*,*)
    write(*,*)'PhL    i_beta ',i_beta
    write(*,*)
    do i_deb=1,2*Nsite
      write(*,*)(phL(i_deb,j_deb,1),j_deb=1,Ntot)
    enddo
    write(*,*)'PhR    i_beta ',i_beta
    write(*,*)
    do i_deb=1,2*Nsite
      write(*,*)(phRbeta(i_deb,j_deb,i_beta),j_deb=1,Ntot)
    enddo
    write(*,*)
    write(*,*)'PhR_prime    i_beta ',i_beta
    write(*,*)
    do i_deb=1,2*Nsite
      write(*,*)(phRbeta_prime(i_deb,j_deb,i_beta),j_deb=1,Ntot)
    enddo
    write(*,*)
    write(*,*)'PhR_sec    i_beta ',i_beta
    write(*,*)
    do i_deb=1,2*Nsite
      write(*,*)(phRbeta_sec(i_deb,j_deb,i_beta),j_deb=1,Ntot)
    enddo
    write(*,*)
    write(*,*)'********** Time instant  ',i_beta
    write(*,*)'ratio_prime(i_beta) ',ratio_prime(i_beta)
    write(*,*)'ratio_sec(i_beta) ',ratio_sec(i_beta)
    write(*,*)'imp_local(1,i_beta) = ',imp_local(1,i_beta)
    write(*,*)'tot_local(i_beta) = ',tot_local(i_beta)
    write(*,*)'w_meas = ',w_meas
    write(*,*)'tot_meas = ',tot_meas
    write(*,*)'rx(i) = ',rx(i)
    write(*,*)'tot_imp(i) = ',tot_imp(i),i
    write(*,*)'weight(i) = ',weight(i)
    write(*,*)'weight(i)*tot_meas/tot_imp(i) = ',weight(i)*tot_meas/tot_imp(i)
    write(*,*)
!DEBUG--------------------------------------------------



    do k=1,Dtot,1

!part B
      call over_lap_dc(phL(1,1,k),phRbeta(1,1,i_beta),ovp)
      call inverse_d_dc(ovp,imp)
      call ZGEMM('N','c',Ntot,2*Nsite,Ntot,one,ovp,Ntot,phL(1,1,k),2*Nsite,zero,wfb_n,Ntot)

!DEBUG--------------------------------------------------
      write(*,*)'<phL | phR> ',imp
!DEBUG--------------------------------------------------

!part A
      call over_lap_dc(phL(1,1,k),phRbeta_prime(1,1,i_beta),ovp_prime)
      call inverse_d_dc(ovp_prime,imp_prime)
      zeta_prime=imp_prime/imp
      call ZGEMM('N','c',Ntot,2*Nsite,Ntot,one,ovp_prime,Ntot,phL(1,1,k),2*Nsite,zero,wfa_n,Ntot)

!DEBUG--------------------------------------------------
      write(*,*)'<phL | phR_prime> ',imp_prime
       write(*,*)'<phL | phR_prime> /<phL | phR>',zeta_prime
!DEBUG--------------------------------------------------


!part C
      call over_lap_dc(phL(1,1,k),phRbeta_sec(1,1,i_beta),ovp_sec)
      call inverse_d_dc(ovp_sec,imp_sec)
      zeta_sec=imp_sec/imp
      call ZGEMM('N','c',Ntot,2*Nsite,Ntot,one,ovp_sec,Ntot,phL(1,1,k),2*Nsite,zero,wfc_n,Ntot)

!DEBUG--------------------------------------------------
      write(*,*)'<phL | phR_prime> ',imp_sec
      write(*,*)'<phL | phR_prime>/<phL | phR> ',zeta_sec
!DEBUG--------------------------------------------------


      do site=1,Nbravais


!DEBUG--------------------------------------------------
        write(*,*)'&&&&&&&&&&&&&&&&&&&&&&&&&'
        write(*,*)'site ',site
        write(*,*)
!DEBUG--------------------------------------------------
     


        part_b_up=zero
        part_b_dn=zero
        do alpha=1,Ntot
          part_b_up=part_b_up+phRbeta(site,alpha,i_beta)*wfb_n(alpha,site)
          part_b_dn=part_b_dn+phRbeta(site+Nsite,alpha,i_beta)*wfb_n(alpha,site+Nsite)
        enddo
!DEBUG--------------------------------------------------
        write(*,*)'<phL | n(site,up) phR>/ <phL | phR> ',part_b_up
        write(*,*)'<phL | n(site,dn) phR>/ <phL | phR> ',part_b_dn
!DEBUG--------------------------------------------------

        part_a_up=zero
        part_a_dn=zero
        do alpha=1,Ntot
          part_a_up=part_a_up+phRbeta_prime(site,alpha,i_beta)*wfa_n(alpha,site)
          part_a_dn=part_a_dn+phRbeta_prime(site+Nsite,alpha,i_beta)*wfa_n(alpha,site+Nsite)
        enddo
!DEBUG--------------------------------------------------
        write(*,*)'<phL | n(site,up) phR_prime>/ <phL | phR_prime> ',part_a_up
        write(*,*)'<phL | n(site,dn) phR_prime>/ <phL | phR_prime> ',part_a_dn
!DEBUG--------------------------------------------------
        part_a_up=part_a_up*ratio_prime(i_beta)*zeta_prime
        part_a_dn=part_a_dn*ratio_prime(i_beta)*zeta_prime


        part_c_up=zero
        part_c_dn=zero
        do alpha=1,Ntot
          part_c_up=part_c_up+phRbeta_sec(site,alpha,i_beta)*wfc_n(alpha,site)
          part_c_dn=part_c_dn+phRbeta_sec(site+Nsite,alpha,i_beta)*wfc_n(alpha,site+Nsite)
        enddo
!DEBUG--------------------------------------------------
        write(*,*)'<phL | n(site,up) phR_sec>/ <phL | phR_sec> ',part_c_up
        write(*,*)'<phL | n(site,dn) phR_sec>/ <phL | phR_sec> ',part_c_dn
!DEBUG--------------------------------------------------
        part_c_up=part_c_up*ratio_sec(i_beta)*zeta_sec
        part_c_dn=part_c_dn*ratio_sec(i_beta)*zeta_sec


        rho_t_local(site,i_beta,k)=(part_a_up + part_c_up - 2.d0*part_b_up +            &
                                    part_a_dn + part_c_dn - 2.d0*part_b_dn)             &
                                   /(Nepero-one)

!DEBUG--------------------------------------------------
        write(*,*)
        write(*,*)'site, i_beta ',site,i_beta
        write(*,*)'rho(site,i_beta) ',rho_t_local(site,i_beta,k)
        write(*,*)
!DEBUG--------------------------------------------------
      
      enddo

    enddo

    call one_step_right(phL,i_beta,i,coe)

  enddo  !end of loop over imaginary time instant


  do i_beta=0,Nbeta,1
    do sitei=1,2*Nsite,1
      call add_numerator(i,imp_local(1,i_beta),tot_local(i_beta),                  &
             coe,w_meas,tot_meas,rho_t_local(sitei,i_beta,1:Dtot),rhom_t(sitei,i_beta))
    enddo
  enddo

  call add_denominator(i,tot_meas,w_meas,denominator)

enddo !end of loop over walkers



!t0=MPI_Wtime()

#ifdef MPI
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(denominator,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
denominator=tmp
do i_beta=0,Nbeta,1
  do sitei=1,Nbravais,1
    call MPI_ALLREDUCE(rhom_t(sitei,i_beta),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
    rhom_t(sitei,i_beta)=tmp
  enddo
enddo
#endif

do i_beta=0,Nbeta,1
  do sitei=1,Nbravais,1
    RHO_one(sitei,i_beta)=rhom_t(sitei,i_beta)/denominator
  enddo
enddo

!t2=MPI_Wtime()
!if(rank.eq.0) then
!      write(*,*)
!      write(*,*) 'GP_one time = ',t2-t0
!      write(*,*)
!  endif


!t0=MPI_Wtime()

deallocate(rho_t_local)
deallocate(phL)

!t2=MPI_Wtime()
!if(rank.eq.0) then
!      write(*,*)
!      write(*,*) 'deallocazione time = ',t2-t0
!      write(*,*)
!  endif


!tff=MPI_Wtime()
!if(rank.eq.0) then
!      write(*,*)
!      write(*,*)
!      write(*,*)
!      write(*,*) 'ALL THE SUBROUTINE time = ',tff-t00
!      write(*,*)
!      write(*,*)
!endif


!DEBUG------------------
!stop 'DEB'
!DEBUG------------------


end subroutine measure_density_beta




