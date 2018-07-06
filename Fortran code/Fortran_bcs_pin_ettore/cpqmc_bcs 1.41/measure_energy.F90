subroutine measure_energy()
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param
use method_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

!denominator and the numerator of the measurement
complex(kind=8)::KEinm,VEinm
complex(kind=8)::denominator  !Determinate the sum of weight

!The local files: <phL|O_local|phr>
complex(kind=8)::KEinm_local(Dtot),VEinm_local(Dtot)

!For the tmp wavefunction file
complex(kind=8),allocatable::phL(:,:,:)
complex(kind=8)::coe(Dtot)
complex(kind=8)::phR(2*Nsite,Ntot)
complex(kind=8),allocatable::ovp_local(:,:,:)
complex(kind=8)::imp_local(Dtot),tot_local
complex(kind=8),allocatable::Amat_local(:,:,:)
complex(kind=8),allocatable::Amat_pure(:,:)

!For the weight meas
complex(kind=8)::w_meas,tot_meas

!The cycle parameter or tmp parameter
integer::i,j,k,sitei,sitej,sitek,ipair,ib,jb,alpha,beta,i_deb,j_deb
complex(kind=8)::tmp,dummy

if(I_wavefun.eq.1)then
  allocate(phL(2*Nsite,Ntot,Dtot))
  allocate(Amat_local(2*Nsite,2*Nsite,Dtot))
elseif(I_wavefun.eq.2)then
  allocate(Amat_pure(2*Nsite,2*Nsite))
endif
allocate(ovp_local(Ntot,Ntot,Dtot))

KEinm=zero
VEinm=zero
denominator=zero

do i=1,Nwalkers,1
   !Get the phL,coe,phR and the ovp_local,imp_local,tot_local
   !also get the measure weight w_meas and tot_meas
   if(I_wavefun.eq.1)then
     call get_meas_arrayE(i,phL,coe,phR,ovp_local,imp_local,tot_local,w_meas,tot_meas)
   elseif(I_wavefun.eq.2)then
     call bcs_get_meas_arrayE(i,phR,coe,w_meas,tot_meas)
   endif


   !To get the Amat(2*Nsite,2*Nsite,1:Dtot)=<phR(1:Dtot)|(ci^+)(cj)|phL>
   if(I_wavefun.eq.1)then
     do k=1,Dtot,1
      !call cal_Amat_withovlpinv(2*Nsite,Ntot,phL(1,1,k),phR(1,1),ovp_local(1,1,k),Amat_local(1,1,k))
        call cal_Amat_withovlpinv_dc(phL(1,1,k),phR(1,1),ovp_local(1,1,k),Amat_local(1,1,k))
     end do
   elseif(I_wavefun.eq.2)then
     call bcs_green_pure(i,phR,Amat_pure,.false.)
   endif


!   write(*,*)
!   write(*,*)
!   write(*,*)'Amat pure measure energy '
!   write(*,*)
!   do i_deb=1,2*Nsite
!     do j_deb=1,2*Nsite
!       write(*,*)j_deb,i_deb,Amat_pure(j_deb,i_deb)
!     enddo
!   enddo
   !stop 'deb'

   !Kinectic energy

   if(I_wavefun.eq.1)then
     call get_klocal(KEinm_local,Amat_local)
   elseif(I_wavefun.eq.2)then
     call bcs_get_klocal(KEinm_local,Amat_pure)
   endif
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,KEinm_local,KEinm)



   !Potential energy
   
   if(I_wavefun.eq.1)then
     call get_vlocal(VEinm_local,Amat_local)
   elseif(I_wavefun.eq.2)then
     call bcs_get_vlocal(VEinm_local,Amat_pure)
   endif
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,VEinm_local,VEinm)

   call add_denominator(i,tot_meas,w_meas,denominator)  
end do




!Get the measurement value to physics_one=numerator/denominator
#ifdef MPI
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(denominator,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
denominator=tmp
call MPI_ALLREDUCE(KEinm,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
KEinm=tmp
call MPI_ALLREDUCE(VEinm,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
VEinm=tmp
#endif

KE_one=KEinm/denominator
VE_one=VEinm/denominator
EE_one=(KEinm+VEinm)/denominator

deallocate(ovp_local)
if(allocated(Amat_local))deallocate(Amat_local)
if(allocated(Amat_pure)) deallocate(Amat_pure)
if(allocated(phL))       deallocate(phL)

end subroutine measure_energy



subroutine bcs_get_meas_arrayE(i,phR,coe,w_meas,tot_meas)
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
complex(kind=8),intent(OUT)::coe(Dtot),phR(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::w_meas,tot_meas
complex(kind=8)::ovp_tmp(Nspin(1),Nspin(1)),tmp,phtmp(2*Nsite,Ntot)
integer::j,k,sitei,sitej

!Get the phR
call k_to_ph_dc(exp_halfK,phi(1,1,i),phR)
coe(:)=one

tot_meas=zero
do k=1,Dtot,1
    if(Nzeta.eq.0)then
      call bcs_over_lap_dc(Fpairing,phR(1,1),ovp_tmp(1:Nspin(1),1:Nspin(1)))
    else
      call unp_bcs_over_lap_dc(Fpairing,Dunpaired,phR(1,1),ovp_tmp(1:Nspin(1),1:Nspin(1)))
    endif
    call inverse_d(ovp_tmp(1:Nspin(1),1:Nspin(1)),Nspin(1),tmp)
    tot_meas=tot_meas+conjg(coe_multi(k))*tmp
enddo

!Get the w_meas
if(crn.LT.0.d0) then
  w_meas=weight(i)*tot_meas/tot_imp(i)
else
  w_meas=weight(i)*one
end if


end subroutine


!------------------------------------------------------------
!We write the wave function after exp_mhalf_K, get the w_meas
!tot_meas, we also get the phL and coe, phR, the overlap 
!information between left and right wave functions.
!------------------------------------------------------------
subroutine get_meas_arrayE(i,phL,coe,phR,ovp_local,imp_local,tot_local,w_meas,tot_meas)
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
complex(kind=8),intent(OUT)::ovp_local(Ntot,Ntot,Dtot)
complex(kind=8),intent(OUT)::imp_local(Dtot),tot_local,w_meas,tot_meas
complex(kind=8)::ovp_tmp(Ntot,Ntot),tmp,phtmp(2*Nsite,Ntot)

integer::j,k,sitei,sitej

!Get the phR
call k_to_ph_dc(exp_halfK,phi(1,1,i),phR)


!Get the phL, coe
call copy_wf_T_dc(phiT(1,1,1),phL(1,1,1))
call zcopy(Dtot,coe_multi(1),1,coe(1),1)

!Get the ovp_local,imp_local,tot_local
tot_local=zero
do k=1,Dtot,1
   !call deter_overlap(2*Nsite,Ntot,phL(1,1,k),phR(1,1),ovp_local(1,1,k))
   call over_lap_dc(phL(1,1,k),phR(1,1),ovp_local(1,1,k))
   call inverse_d_dc(ovp_local(1:Ntot,1:Ntot,k),imp_local(k))
   tot_local=tot_local+conjg(coe(k))*imp_local(k)
end do



!Get the tot_meas
tot_meas=tot_local


!Get the w_meas
if(crn.LT.0.d0) then
  w_meas=weight(i)*tot_meas/tot_imp(i)
else 
  w_meas=weight(i)*one
end if

end subroutine get_meas_arrayE

