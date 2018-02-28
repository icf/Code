
!--------------------------------------------------
!Initial rank and Nsize, if mpi then inital the MPI
!--------------------------------------------------
subroutine start_code()
use mpi_serial_param
use io_module
implicit none
#ifdef MPI
include "mpif.h"
#endif

#ifdef MPI
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,Nsize,ierr)
#else
 rank=0
 Nsize=1
#endif
end subroutine start_code


!----------------------------------------------------
!This subroutine is used to read parameter from param
!at the beginning of running the code++++++++++++++++
!----------------------------------------------------
subroutine readparam
use lattice_param
use model_param
use method_param
use mc_loop_param
use project_param
use phiT_param
use phi_param
use adET_param
use mpi_serial_param
use sc_loop_param
implicit none
open(unit=10,file='param',status='old')

!read the lattice parameter
read(10,*) set
read(10,*) Nsite
read(10,*) Nhop
read(10,*) Dimen
read(10,*) Nl(1)
read(10,*) Nl(2)
read(10,*) Nl(3)
read(10,*) kbound(1)
read(10,*) kbound(2)
read(10,*) kbound(3)

!read the model parameter
read(10,*) t1
read(10,*) onsitU
read(10,*) Ntot
read(10,*) dtype
read(10,*) Nspin(1)
read(10,*) Nspin(2)

if(dtype.EQ.'d') then
  Ntot=Nspin(1)+Nspin(2)
end if



!The kind of methods
read(10,*) crn
read(10,*) ccoe
read(10,*) max_crn
read(10,*) dcp
read(10,*) kcrn
read(10,*) bgset
read(10,*) fw_bk
read(10,*) Nbk(1)
read(10,*) Nbk(2)
Nstps_fwd=Nbk(1)+Nbk(2)
!Number of loop steps
read(10,*) sc_loop_num

!MC loop parameter
read(10,*) Nsamples
read(10,*) Nwalkers
#ifdef MPI
 MNwalkers=Nwalkers*Nsize
#endif
read(10,*) blockstep
read(10,*) Thermblock
read(10,*) Neqblock
read(10,*) PopContrlstep
read(10,*) StepforGram
read(10,*) meastep
max_local=Neqblock*blockstep

!Project parameter
read(10,*) dt

!PhiT parameter
read(10,*) PT

!Phi parameter
read(10,*) PP

!adET parameter
read(10,*) max_ad

close(10)
end subroutine readparam


!------------------------------------------------
!This subroutine initial the mpi htype and phtype
!------------------------------------------------
subroutine init_mpi_type()
use model_param
use mpi_serial_param
use lattice_param
implicit none
integer::bl(2*Nsite)
integer::disp(2*Nsite)
integer::i
#ifdef MPI
include "mpif.h"
#endif


#ifdef MPI
 if(dtype.EQ.'c') then
   call MPI_TYPE_CONTIGUOUS(4*Nsite*Nsite,MPI_DOUBLE_COMPLEX,htype,ierr)
   call MPI_TYPE_CONTIGUOUS(2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,phtype,ierr)
 else if(dtype.EQ.'d') then

   bl=Nsite
   do i=1,Nsite,1
      disp(i)=(i-1)*2*Nsite
   end do
   do i=1,Nsite,1
      disp(i+Nsite)=(Nsite+i-1)*2*Nsite+Nsite
   end do
   call MPI_TYPE_INDEXED(2*Nsite,bl,disp,MPI_DOUBLE_COMPLEX,htype,ierr)


   bl=Nsite
   do i=1,Nspin(1),1
      disp(i)=(i-1)*2*Nsite
   end do
   do i=1,Nspin(2),1
      disp(i+Nspin(1))=(Nspin(1)+i-1)*2*Nsite+Nsite
   end do
   call MPI_TYPE_INDEXED(Ntot,bl,disp,MPI_DOUBLE_COMPLEX,phtype,ierr)

 end if
 call MPI_TYPE_COMMIT(htype,ierr)
 call MPI_TYPE_COMMIT(phtype,ierr)
#else
return
#endif

end subroutine init_mpi_type


!-----------------------------------------------------------------
!This subroutine allocate the arrays we need to use in set lattice
!-----------------------------------------------------------------
subroutine allocate_lattice_array()
use lattice_param
implicit none
allocate(hopt(Nhop))
allocate(sit(Nhop,2))
allocate(coor(Nsite,Dimen))
allocate(Hzero(2*Nsite,2*Nsite))
allocate(Tmatrix(Nsite,Dimen))
end subroutine allocate_lattice_array


!-------------------------------------------------------------------
!This subroutine deallocate the arrays we need to use in set lattice
!-------------------------------------------------------------------
subroutine deallocate_lattice_array()
use lattice_param
implicit none
if(allocated(hopt)) deallocate(hopt)
if(allocated(sit)) deallocate(sit)
if(allocated(coor)) deallocate(coor)
if(allocated(Hzero)) deallocate(Hzero)
if(allocated(Tmatrix)) deallocate(Tmatrix)
end subroutine deallocate_lattice_array



!-------------------------------------------------------------
!This subroutine allocate the arrays in the projection K and V
!------------------------------------------------------------- 
subroutine allocate_project_array()
use lattice_param
use project_param
implicit none
allocate(exp_halfK(2*Nsite,2*Nsite))
allocate(exp_mhalfK(2*Nsite,2*Nsite))
allocate(exp_K(2*Nsite,2*Nsite))
allocate(ng(Nsite))
end subroutine allocate_project_array



!---------------------------------------------------------------
!This subroutine deallocate the arrays in the projection K and V
!--------------------------------------------------------------- 
subroutine deallocate_project_array()
use project_param
implicit none
if(allocated(exp_halfK)) deallocate(exp_halfK)
if(allocated(exp_mhalfK)) deallocate(exp_mhalfK)
if(allocated(exp_K)) deallocate(exp_K)
if(allocated(ng)) deallocate(ng)
end subroutine deallocate_project_array


!---------------------------------------------------------
!This subroutine allocate the arrays when we need for phiT
!---------------------------------------------------------
subroutine allocate_phiT()
use lattice_param
use model_param
use phiT_param
implicit none
allocate(phiT(2*Nsite,Ntot,Dtot))
allocate(coe_multi(Dtot))
end subroutine allocate_phiT


!------------------------------------------------------------
!This subroutine deallocate the arrays where we need for phiT
!------------------------------------------------------------
subroutine deallocate_phiT()
use phiT_param
implicit none
if(allocated(phiT)) deallocate(phiT)
if(allocated(coe_multi))  deallocate(coe_multi)
end subroutine deallocate_phiT



!--------------------------------------
!This subroutine allocate arrays in phi
!--------------------------------------
subroutine allocate_phi()
use lattice_param
use model_param
use mc_loop_param
use method_param
use phiT_param
use phi_param
implicit none
allocate(phi(2*Nsite,Ntot,Nwalkers))
!allocate(invop(Ntot,Ntot,Nwalkers))
allocate(impfunc(Dtot,Nwalkers))
allocate(tot_imp(Nwalkers))
allocate(weight(Nwalkers))
allocate(dlogw(Nwalkers))
allocate(rx(Nwalkers))

if(fw_bk.NE.'FW') then
  if(max_crn.GE.0) allocate(back_store(Nsite,Nstps_fwd,Nwalkers))
  allocate(phi0(2*Nsite,Ntot,Nwalkers))
end if
end subroutine allocate_phi


!----------------------------------------
!This subroutine deallocate arrays in phi
!----------------------------------------
subroutine deallocate_phi()
use phi_param
implicit none
if(allocated(phi)) deallocate(phi)
!if(allocated(invop)) deallocate(invop)
if(allocated(impfunc)) deallocate(impfunc)
if(allocated(tot_imp)) deallocate(tot_imp)
if(allocated(weight)) deallocate(weight)
if(allocated(Mweight)) deallocate(Mweight)
if(allocated(dlogw)) deallocate(dlogw)
if(allocated(rx)) deallocate(rx)
if(allocated(back_store)) deallocate(back_store)
if(allocated(phi0)) deallocate(phi0)
end subroutine deallocate_phi


!---------------------------------------------
!This subroutine allocate the arrays in backup
!---------------------------------------------
subroutine allocate_backup()
use backup_param
use lattice_param
use model_param
use mc_loop_param
use phiT_param
implicit none
allocate(phi_tmp(2*Nsite,Ntot,Nwalkers))
allocate(weight_tmp(Nwalkers))
allocate(dlogw_tmp(Nwalkers))
allocate(rx_tmp(Nwalkers))
allocate(impfunc_tmp(Dtot,Nwalkers))
allocate(tot_imp_tmp(Nwalkers))
end subroutine allocate_backup


!-----------------------------------------------
!This subroutine deallocate the arrays in backup
!-----------------------------------------------
subroutine deallocate_backup()
use backup_param
implicit none
if(allocated(phi_tmp)) deallocate(phi_tmp)
if(allocated(weight_tmp)) deallocate(weight_tmp)
if(allocated(dlogw_tmp)) deallocate(dlogw_tmp)
if(allocated(rx_tmp)) deallocate(rx_tmp)
if(allocated(impfunc_tmp)) deallocate(impfunc_tmp)
if(allocated(tot_imp_tmp)) deallocate(tot_imp_tmp)
end subroutine deallocate_backup



!-----------------------------------------
!allocate the arrays in measure subroutine
!-----------------------------------------
subroutine allocate_meas()
use mc_loop_param
use method_param
use lattice_param
use meas_param
implicit none
integer::max_n
if(max_crn.GE.0) then
  max_n=max(max_crn,0)
  if(fw_bk.NE.'FW') max_n=0 !for release back propogation.
  allocate(kin_l(Nsamples,0:max_n))
  allocate(v_l(Nsamples,0:max_n))
  allocate(e_l(Nsamples,0:max_n))
  allocate(var_l(Nsamples,0:max_n))
  allocate(nu_l(Nsamples,0:max_n))
  allocate(nd_l(Nsamples,0:max_n))
  allocate(sisj_l(Nsamples,Nsite,0:max_n))
  allocate(sk_l(Nsamples,Nsite,0:max_n))
  allocate(cicj_l(Nsamples,2*Nsite,0:max_n))
  allocate(cicj_l_global(Nsamples,2*Nsite,2*Nsite,0:max_n))
  allocate(ck_l(Nsamples,2*Nsite,0:max_n))
  allocate(sig(0:max_n))
  allocate(absig(0:max_n))
else if(max_crn.LT.0) then
  !allocate(e_l(Nsamples,max_local))
  !allocate(var_l(Nsamples,max_local))
  !allocate(nu_l(Nsamples,max_local))
  !allocate(nd_l(Nsamples,max_local))
  !allocate(sig(max_local))
  !allocate(absig(max_local))
  call allocate_meas_step()
else
  write(*,*) "Something is wrong with max_crn input."
  call mystop
end if
end subroutine allocate_meas



!---------------------------------------------------------
!allocate the arrays in free projection measure subroutine
!---------------------------------------------------------
subroutine allocate_meas_step()
use mc_loop_param
use method_param
use lattice_param
use meas_param
implicit none
  allocate(kin_l(Nsamples,max_local))
  allocate(v_l(Nsamples,max_local))
  allocate(e_l(Nsamples,max_local))
  allocate(var_l(Nsamples,max_local))
  allocate(nu_l(Nsamples,max_local))
  allocate(nd_l(Nsamples,max_local))
  allocate(sisj_l(Nsamples,Nsite,max_local))
  allocate(sk_l(Nsamples,Nsite,max_local))
  allocate(cicj_l(Nsamples,2*Nsite,max_local))
  allocate(cicj_l_global(Nsamples,2*Nsite,2*Nsite,max_local))
  allocate(ck_l(Nsamples,2*Nsite,max_local))
  allocate(sig(max_local))
  allocate(absig(max_local))
end subroutine allocate_meas_step


!-------------------------------------------
!deallocate the arrays in measure subroutine
!-------------------------------------------
subroutine deallocate_meas()
use meas_param
implicit none
if(allocated(kin_l)) deallocate(kin_l)
if(allocated(v_l)) deallocate(v_l)
if(allocated(e_l)) deallocate(e_l)
if(allocated(var_l)) deallocate(var_l)
if(allocated(nu_l)) deallocate(nu_l)
if(allocated(nd_l)) deallocate(nd_l)
if(allocated(sisj_l)) deallocate(sisj_l)
if(allocated(sk_l)) deallocate(sk_l)
if(allocated(cicj_l)) deallocate(cicj_l)
if(allocated(cicj_l_global)) deallocate(cicj_l_global)
if(allocated(ck_l)) deallocate(ck_l)
if(allocated(sig)) deallocate(sig)
if(allocated(absig)) deallocate(absig)
end subroutine deallocate_meas


!--------------------------------------------------
!This subroutine allocate the arrays in one measure
!--------------------------------------------------
subroutine allocate_one_meas()
use one_meas_param
use lattice_param
implicit none
allocate(S_one(Nsite))
allocate(c_one(2*Nsite))
allocate(c_one_global(2*Nsite,2*Nsite))
end subroutine allocate_one_meas



!----------------------------------------------------
!This subroutine deallocate the arrays in one measure
!----------------------------------------------------
subroutine deallocate_one_meas()
use one_meas_param
implicit none
if(allocated(S_one)) deallocate(S_one)
if(allocated(c_one)) deallocate(c_one)
if(allocated(c_one_global)) deallocate(c_one_global)
end subroutine deallocate_one_meas


!--------------------------------------------------
!This subroutine allocate the S.C. 
!--------------------------------------------------
subroutine allocate_sc_loop_param()
use sc_loop_param
use lattice_param
use model_param
use phiT_param
implicit none
allocate(cicj_sc_global(2*Nsite,2*Nsite))
allocate(cicj_sc_global1(2*Nsite,2*Nsite))
allocate(BCS_sc(Nsite,Nsite))
allocate(phi_sc(Nsite,Ntot,Dtot))
end subroutine allocate_sc_loop_param

!----------------------------------------------------
!This subroutine deallocate the S.C.
!----------------------------------------------------
subroutine deallocate_sc_loop_param()
use sc_loop_param
use lattice_param
use model_param
use phiT_param
implicit none

if(allocated(cicj_sc_global)) then
   deallocate(cicj_sc_global)
endif

if(allocated(cicj_sc_global1)) then
   deallocate(cicj_sc_global1)
endif

if(allocated(BCS_sc)) then
   deallocate(BCS_sc)
endif

if(allocated(phi_sc)) then 
   deallocate(phi_sc)
endif

end subroutine deallocate_sc_loop_param


!-----------------------------------------------------------------
!This subroutine deallocate all the arrays at the end of the code.
!-----------------------------------------------------------------
subroutine deallocatearray()
implicit none
call deallocate_lattice_array()
call deallocate_project_array()
call deallocate_phiT()
call deallocate_phi()
call deallocate_backup()
call deallocate_meas()
call deallocate_one_meas()
end subroutine deallocatearray



!-------------------------------------------------
!This subroutine clean call the things before stop
!-------------------------------------------------
subroutine clean
!------------------------------------------------
! without sprng_f.h. skip rand_num, icf, 1/1/2018
!------------------------------------------------
use timing_module
!use rand_num
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
!call end_genrand()
call deallocatearray()
call EndTiming()
#ifdef MPI
call MPI_TYPE_FREE(htype,ierr)
call MPI_TYPE_FREE(phtype,ierr)
call MPI_Finalize(ierr)
#endif
end subroutine clean



!----------------------------------------------------------------
!This subroutine end the program before deallocate all the arrays
!----------------------------------------------------------------
subroutine mystop
implicit none
call clean()
stop
end subroutine mystop
