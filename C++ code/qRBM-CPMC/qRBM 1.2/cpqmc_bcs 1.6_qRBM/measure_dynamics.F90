subroutine measure_one_body_dynamics()
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

!denominator and the numerator of the measurement
complex(kind=8)::cicjm_t(2*Nsite,Nbeta)
complex(kind=8)::cicjhm_t(2*Nsite,Nbeta)
complex(kind=8)::denominator  !Determinate the sum of weight

!The local files: <phL|O_local|phr>
complex(kind=8),allocatable::cicj_t_local(:,:,:)
complex(kind=8),allocatable::cicjh_t_local(:,:,:)

!For the tmp wavefunction file
complex(kind=8),allocatable::phL(:,:,:)
complex(kind=8)::coe(Dtot)
complex(kind=8)::phR(2*Nsite,Ntot),phRstart(2*Nsite,Ntot)
complex(kind=8)::phRbeta(2*Nsite,Ntot,0:Nbeta)
complex(kind=8)::B_matrix(2*Nsite,2*Nsite),Gp_matrix(2*Nsite,2*Nsite)
complex(kind=8)::Bh_matrix(2*Nsite,2*Nsite),Gh_matrix(2*Nsite,2*Nsite)
complex(kind=8)::Green_dynamic(2*Nsite,2*Nsite,Nbeta,Dtot)
complex(kind=8)::Green_dynamich(2*Nsite,2*Nsite,Nbeta,Dtot)


complex(kind=8)::imp_local(Dtot,Nbeta),tot_local(Nbeta)
complex(kind=8),allocatable::Amat_local(:,:,:)

!For the weight meas
complex(kind=8)::w_meas,tot_meas

!The cycle parameter or tmp parameter
integer::i,j,k,sitei,sitej,aux,ib,jb,kb,i_beta
integer:: i_deb,j_deb,k_deb !DEBUG
real(kind=8)::x
complex(kind=8)::explr_up,explr_dn,tmp
!real(kind=8)::t0,t2,t3,t4,t00,tff


allocate(cicj_t_local(2*Nsite,Nbeta,Dtot))
allocate(cicjh_t_local(2*Nsite,Nbeta,Dtot))

allocate(phL(2*Nsite,Ntot,Dtot))
allocate(Amat_local(2*Nsite,2*Nsite,Dtot))

cicjm_t=zero
cicjhm_t=zero
denominator=zero

do i=1,Nwalkers,1   !loop over walkers

  call get_phiR_phiL(i,phL,coe,phR) 
!Backup for phR
  phRstart(:,:)=phR(:,:)
!zero imaginary time
  call copy_wf_dc(phR(1,1),phRbeta(1,1,0)) 

  do i_beta=1,Nbeta    !loop over imaginary time instants

    call copy_wf_dc(phRbeta(1,1,i_beta-1),phR(1,1)) 
    call one_step_left(phR,i_beta,i)
    call copy_wf_dc(phR(1,1),phRbeta(1,1,i_beta))
     


!DEB
!    write(*,*)
!    write(*,*)'phRbeta  ',i_beta
!    write(*,*)
!    do i_deb=1,2*Nsite
!     write(*,*)(phRbeta(i_deb,j_deb,i_beta),j_deb=1,Ntot)
!    enddo



  enddo !end of loop over imaginary time instants

!  stop'DEB'

  do i_beta=Nbeta,1,-1    !loop over imaginary time instants

    call get_useful_params(i,phL,coe,phRstart,                      &
                         imp_local(1,i_beta),tot_local(i_beta),w_meas,tot_meas)

    call calculate_B(i_beta,i,B_matrix)
    call calculate_Bh(i_beta,i,Bh_matrix)



!DEB
!    if(dtype.EQ.'c') then
!    write(*,*)
!    write(*,*)'Bh_matrix  '
!    write(*,*)
!    do i_deb=1,2*Nsite
!      do j_deb=1,2*Nsite
!        write(*,*)i_deb,j_deb,Bh_matrix(i_deb,j_deb)
!      enddo
!    enddo
!    endif
!    stop'DEB'




    do k=1,Dtot,1


      call calculate_Gp(phL(1,1,k),phRbeta(1,1,i_beta),Gp_matrix(1,1))
      if(dtype.EQ.'c') then
        do ib=1,2*Nsite
          do jb=1,2*Nsite
            Gh_matrix(ib,jb)=-Gp_matrix(jb,ib)
          enddo
          Gh_matrix(ib,ib)=Gh_matrix(ib,ib)+one
        enddo
      else if(dtype.EQ.'d') then
        do ib=1,Nsite
          do jb=1,Nsite
            Gh_matrix(ib,jb)=-Gp_matrix(jb,ib)
            Gh_matrix(ib+Nsite,jb+Nsite)=-Gp_matrix(jb+Nsite,ib+Nsite)
          enddo
          Gh_matrix(ib+Nsite,ib+Nsite)=Gh_matrix(ib+Nsite,ib+Nsite)+one
        enddo
      endif

!DEB
!      write(*,*)
!      write(*,*)'Gh_matrix  '
!      write(*,*)
!      do i_deb=1,2*Nsite
!        do j_deb=1,2*Nsite
!          write(*,*)i_deb,j_deb,Gh_matrix(i_deb,j_deb)
!        enddo
!      enddo

 


      if(dtype.EQ.'c') then
        call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,Gp_matrix(1,1),2*Nsite,B_matrix(1,1)                 &
     &            ,2*Nsite,zero,Green_dynamic(1,1,i_beta,k),2*Nsite)

        call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,Bh_matrix(1,1),2*Nsite,Gh_matrix(1,1)                &
     &            ,2*Nsite,zero,Green_dynamich(1,1,i_beta,k),2*Nsite)
      else if(dtype.EQ.'d') then
        call zgemm('N','N',Nsite,Nsite,Nsite,one,Gp_matrix(1,1),2*Nsite,B_matrix(1,1)                       &
     &            ,2*Nsite,zero,Green_dynamic(1,1,i_beta,k),2*Nsite)
        call zgemm('N','N',Nsite,Nsite,Nsite,one,Gp_matrix(Nsite+1,Nsite+1),2*Nsite                         &
                  ,B_matrix(Nsite+1,Nsite+1),2*Nsite,zero,Green_dynamic(Nsite+1,Nsite+1,i_beta,k),2*Nsite)

        call zgemm('N','N',Nsite,Nsite,Nsite,one,Bh_matrix(1,1),2*Nsite,Gh_matrix(1,1)                       &
     &            ,2*Nsite,zero,Green_dynamich(1,1,i_beta,k),2*Nsite)
        call zgemm('N','N',Nsite,Nsite,Nsite,one,Bh_matrix(Nsite+1,Nsite+1),2*Nsite                         &
                  ,Gh_matrix(Nsite+1,Nsite+1),2*Nsite,zero,Green_dynamich(Nsite+1,Nsite+1,i_beta,k),2*Nsite)
        
      endif


    enddo


    call one_step_right(phL,i_beta,i,coe)


  enddo  !end of loop over imaginary time instant


  call measure_green_particles(Green_dynamic(1,1,1,1),cicj_t_local(1,1,1))
  call measure_green_holes(Green_dynamich(1,1,1,1),cicjh_t_local(1,1,1))

!DEB
!  write(*,*)'cicjh_t_local '
!  do i_deb=1,Nbeta
!    do j_deb=1,1
!      write(*,*)i_deb*dt,dble(cicjh_t_local(j_deb,i_deb,1)),j_deb
!    enddo
!  enddo

!  stop'DEB'





  do i_beta=1,Nbeta,1
    do sitei=1,2*Nsite,1
      call add_numerator(i,imp_local(1,i_beta),tot_local(i_beta),                  &
             coe,w_meas,tot_meas,cicj_t_local(sitei,i_beta,1:Dtot),cicjm_t(sitei,i_beta))

      call add_numerator(i,imp_local(1,i_beta),tot_local(i_beta),                  &
             coe,w_meas,tot_meas,cicjh_t_local(sitei,i_beta,1:Dtot),cicjhm_t(sitei,i_beta))
    enddo
  enddo

  call add_denominator(i,tot_meas,w_meas,denominator)

enddo !end of loop over walkers



!t0=MPI_Wtime()

#ifdef MPI
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(denominator,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
denominator=tmp

do i_beta=1,Nbeta,1
  do sitei=1,2*Nsite,1
    call MPI_ALLREDUCE(cicjm_t(sitei,i_beta),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
    cicjm_t(sitei,i_beta)=tmp
  enddo
enddo

do i_beta=1,Nbeta,1
  do sitei=1,2*Nsite,1
    call MPI_ALLREDUCE(cicjhm_t(sitei,i_beta),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
    cicjhm_t(sitei,i_beta)=tmp
  enddo
enddo


#endif

!t2=MPI_Wtime()
!  if(rank.eq.0) then
!      write(*,*)
!      write(*,*) 'MPI_ALLREDUCE time = ',t2-t0
!      write(*,*)
!  endif

!t0=MPI_Wtime()


do i_beta=1,Nbeta,1
  do sitei=1,2*Nsite,1
    GP_one(sitei,i_beta)=cicjm_t(sitei,i_beta)/denominator
  enddo
enddo

do i_beta=1,Nbeta,1
  do sitei=1,2*Nsite,1
    GH_one(sitei,i_beta)=cicjhm_t(sitei,i_beta)/denominator
  enddo
enddo

!t2=MPI_Wtime()
!if(rank.eq.0) then
!      write(*,*)
!      write(*,*) 'GP_one time = ',t2-t0
!      write(*,*)
!  endif


!t0=MPI_Wtime()

deallocate(cicj_t_local)
deallocate(cicjh_t_local)
deallocate(phL)
deallocate(Amat_local)

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



end subroutine measure_one_body_dynamics







!------------------------------------------------------------
!We write the wave function after exp_mhalf_K, get the w_meas
!tot_meas, we also get the phL and coe, phR, the overlap 
!information between left and right wave functions.
!------------------------------------------------------------
subroutine get_phiR_phiL(i,phL,coe,phR)
use param
use lattice_param
use model_param
use project_param
use phiT_param
use phi_param
use method_param
use caldet_module
use mpi_serial_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(OUT)::phL(2*Nsite,Ntot,Dtot),coe(Dtot),phR(2*Nsite,Ntot)

if(.not.back_pro) then
  write(*,*)"Please use back-propagation if you want dynamical properties"
  stop
endif


!Get the phR
call k_to_ph_dc(exp_halfK,phi0(1,1,i),phR)

!Get the phL, coe
call back_prog_phL2(i,phL,coe)

end subroutine get_phiR_phiL


!---------------------------------------------------------------------
subroutine get_other_params(i,phL,coe,phR,ovp_local,imp_local,tot_local,w_meas,tot_meas)
use param
use lattice_param
use model_param
use project_param
use phiT_param
use phi_param
use method_param
use caldet_module
use mpi_serial_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(IN)::phL(2*Nsite,Ntot,Dtot),coe(Dtot),phR(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::ovp_local(Ntot,Ntot,Dtot)
complex(kind=8),intent(OUT)::imp_local(Dtot),tot_local,w_meas,tot_meas
complex(kind=8)::ovp_tmp(Ntot,Ntot),tmp,phtmp(2*Nsite,Ntot)

integer::j,k

!Get the ovp_local,imp_local,tot_local
tot_local=zero
do k=1,Dtot,1
   !call deter_overlap(2*Nsite,Ntot,phL(1,1,k),phR(1,1),ovp_local(1,1,k))
   call over_lap_dc(phL(1,1,k),phR(1,1),ovp_local(1,1,k))
   !call caldet_dc(ovp_local(1:Ntot,1:Ntot,k),imp_local(k))
   !call inverse_dc(ovp_local(1:Ntot,1:Ntot,k))
   call inverse_d_dc(ovp_local(1:Ntot,1:Ntot,k),imp_local(k))
   tot_local=tot_local+conjg(coe(k))*imp_local(k)
end do



!Get the tot_meas
if(.not.back_pro) then
  tot_meas=tot_local
else if(back_pro) then
  !call
  !zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi(1,1,i),2*Nsite,zero,phtmp,2*Nsite)
  call k_to_ph_dc(exp_halfK,phi(1,1,i),phtmp)
  tot_meas=zero
  do k=1,Dtot,1
     !call deter_overlap(2*Nsite,Ntot,phiT(1,1,k),phtmp(1,1),ovp_tmp(1,1))
     call over_lap_dc(phiT(1,1,k),phtmp(1,1),ovp_tmp(1,1))
     !call caldet(Ntot,ovp_tmp(1:Ntot,1:Ntot),tmp)
     call caldet_dc(ovp_tmp(1:Ntot,1:Ntot),tmp)
     tot_meas=tot_meas+conjg(coe_multi(k))*tmp
  end do
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if


!Get the w_meas
if(crn.LT.0.d0) then
  w_meas=weight(i)*tot_meas/tot_imp(i)
!  write(*,*)'tot_meas/tot_imp(i) ',tot_meas/tot_imp(i),'AAAAAAAAAAAAAAAAAAAAAAA'
else
  w_meas=weight(i)*one
end if



end subroutine get_other_params



!----------------------------------------
!This subroutine propagate the phi to phL
!---------------------------------------- 
subroutine back_prog_phL2(i,phL,coe)
use param
use lattice_param
use model_param
use method_param
use phiT_param
use mc_loop_param
use phi_param
use project_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(OUT)::phL(2*Nsite,Ntot,Dtot)
complex(kind=8),intent(OUT)::coe(Dtot)
complex(kind=8)::phtmp(2*Nsite,Ntot)
integer::i_gs
integer::j,k,i_b

do k=1,Dtot,1
   !call
   !zgemm('C','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phiT(1,1,k),2*Nsite,zero,phL(1,1,k),2*Nsite)
   call dk_to_ph_dc(exp_halfK,phiT(1,1,k),phL(1,1,k))
end do
call zcopy(Dtot,coe_multi(1),1,coe(1),1)


i_gs=0
do i_b=Nstps_fwd,Nbeta+1,-1
   call back_one(i_b,i,phL)
   i_gs=i_gs+1
   if(i_gs.eq.StepforGram) then
     i_gs=0
     call gs_phL(phL,coe)
   end if
end do

do k=1,Dtot,1
   call copy_wf_dc(phL(1,1,k),phtmp(1,1))
   call dk_to_ph_dc(exp_mhalfK,phtmp(1,1),phL(1,1,k))
end do
end subroutine back_prog_phL2

