
!--------------------------------------------------
!Initial rank and Nsize, if mpi then inital the MPI
!--------------------------------------------------
subroutine start_code()
use mpi_serial_param
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
implicit none
integer::temp
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
read(10,*) boundary

!read the model parameter
read(10,*) I_ReadH0
read(10,*) Nbands
if(Nbands.eq.1)then
  read(10,*) t1
  read(10,*) Hpinn
  ipinn=1
  if(abs(Hpinn).lt.0.000001)then
    ipinn=0
    if(rank.eq.0)write(*,*)'No pinning field'
    Hpinn=0.d0
  else
    if(rank.eq.0)write(*,*)'A pinning field is applied'
  endif
  read(10,*) pinndir
  read(10,*) onsitU
elseif(Nbands.eq.3)then
  read(10,*) tpp
  read(10,*) tpd
  read(10,*) epsp
  read(10,*) epsd
  read(10,*) Hpinn 
  ipinn=1
  if(abs(Hpinn).lt.0.000001)then
    ipinn=0    
    if(rank.eq.0)write(*,*)'No pinning field'
    Hpinn=0.d0
  else
    if(rank.eq.0)write(*,*)'A pinning field is applied'
  endif 
  read(10,*) pinndir
  read(10,*) onsitUp
  read(10,*) onsitUd
  read(10,*) Vpd
else
  write(*,*)'Problems with Nbands'
  stop
endif
read(10,*) I_Vext
if(I_Vext.eq.1)ipinn=1
read(10,*) Ntot
read(10,*) dtype
read(10,*) Nspin(1)
read(10,*) Nspin(2)
if(Nspin(2).gt.Nspin(1))then
  temp=Nspin(1)
  Nspin(1)=Nspin(2)
  Nspin(2)=temp
endif
Nzeta=Nspin(1)-Nspin(2)

if(dtype.EQ.'d'.or.dtype.EQ.'m') then
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
read(10,*) StepforPure
read(10,*) meastep
max_local=Neqblock*blockstep

!Project parameter
read(10,*) dt

!Kind of trial wave function, HF or BCS
read(10,*) I_wavefun
!PhiT parameter
read(10,*) PT
!Phi parameter
read(10,*) PP

!adET parameter
read(10,*) max_ad

read(10,*) I_obdm

read(10,*) Npair

read(10,*) Nbeta

read(10,*) I_twob

read(10,*) I_onebf

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
   call MPI_TYPE_CONTIGUOUS(2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,phtype2,ierr)
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
   call MPI_TYPE_INDEXED(Ntot,bl,disp,MPI_DOUBLE_COMPLEX,phtype2,ierr)

 else if(dtype.EQ.'m') then

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

   call MPI_TYPE_CONTIGUOUS(2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,phtype2,ierr)
   
 end if
 call MPI_TYPE_COMMIT(htype,ierr)
 call MPI_TYPE_COMMIT(phtype,ierr)
 call MPI_TYPE_COMMIT(phtype2,ierr)
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
if(I_wavefun.eq.2)then
  allocate(FPairing(Nsite,Nsite))
  if(Nzeta.gt.0)then
    allocate(DUnpaired(Nsite,Nzeta))
  endif
endif
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
use model_param
use meas_param
implicit none
integer::max_n
if(max_crn.GE.0) then
  max_n=max(max_crn,0)
  if(fw_bk.NE.'FW') max_n=0 !for release back propogation.
  allocate(ZetaN_l(Nsamples,Dimen,0:max_n))
  allocate(kin_l(Nsamples,0:max_n))
  allocate(v_l(Nsamples,0:max_n))
  allocate(e_l(Nsamples,0:max_n))
  allocate(eE_l(Nsamples,0:max_n))
  allocate(var_l(Nsamples,0:max_n))
  allocate(nu_l(Nsamples,0:max_n))
  allocate(nd_l(Nsamples,0:max_n))
  if(I_obdm.eq.1)allocate(obdm_l(Nsamples,2*Nsite,2*Nsite,0:max_n))
!  if(ipinn.eq.0)allocate(cacb_l(Nsamples,Nbravais,Nbands,Nbands,0:max_n))
  if(ipinn.eq.1)allocate(nofr_l(Nsamples,Nsite,2,2,0:max_n))
  if(ipinn.eq.0)allocate(sisj_l(Nsamples,Nbravais,0:max_n))
  if(ipinn.eq.0)allocate(ninj_l(Nsamples,Nbravais,Nbands,Nbands,0:max_n))
  if(ipinn.eq.0)allocate(szsz_l(Nsamples,Nbravais,Nbands,Nbands,0:max_n))
  if(ipinn.eq.0)allocate(sk_l(Nsamples,Nsite,0:max_n))
!  allocate(cicj_l(Nsamples,2*Nsite,0:max_n))
!  allocate(ck_l(Nsamples,2*Nsite,0:max_n))
  if(Npair.gt.0)then
   allocate(didj_l(Nsamples,Nbravais,Npair,0:max_n))
  endif
  if(Nbeta.gt.0)then
   if(I_onebf.eq.1)then
     allocate(cicj_t_l(Nsamples,2*Nsite,Nbeta,0:max_n))
     allocate(cicjh_t_l(Nsamples,2*Nsite,Nbeta,0:max_n))
   endif
   allocate(GreenP_t_l(Nsamples,0:Nbeta,0:max_n))
   allocate(GreenH_t_l(Nsamples,0:Nbeta,0:max_n))
 !  allocate(rho_t_l(Nsamples,Nbravais,0:Nbeta,0:max_n))
   if(I_twob.eq.1)then
     allocate(nupnup_t_l(Nsamples,Nbravais,0:Nbeta,0:max_n))
     allocate(ndnnup_t_l(Nsamples,Nbravais,0:Nbeta,0:max_n))
   endif
  endif
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



subroutine set_beginning_zero()
use param
use mc_loop_param
use mpi_serial_param
use adET_param
use meas_param
use mpi_serial_param
use method_param
use meas_param
use model_param
implicit none

!properties to be always measured
 kin_l=0.d0;v_l=0.d0;e_l=0.d0;eE_l=zero;var_l=zero;nu_l=0.d0;nd_l=0.d0
 sig=zero;absig=zero

 if(I_obdm.eq.1)obdm_l=zero
 ZetaN_l=zero
!propertied for periodic systems
 if(ipinn.eq.0)then
  sisj_l=0.d0;ninj_l=0.d0;szsz_l=0.d0;sk_l=0.d0
!cacb_l=zero
 endif
!inhomogenerous systems, external fields
 if(ipinn.eq.1)then
  nofr_l=zero
 endif
 !cicj_l=zero;ck_l=0.d0
!pairing correlations
 if(Npair.gt.0)didj_l=zero
!dynamical properties
 if(Nbeta.gt.0)then
   if(I_onebf.eq.1)then
     cicj_t_l=zero
     cicjh_t_l=zero
   endif
  ! rho_t_l=zero
   if(I_twob.eq.1)then
     nupnup_t_l=zero
     ndnnup_t_l=zero
   endif
   GreenP_t_l=zero
   GreenH_t_l=zero
 endif

end subroutine set_beginning_zero




!---------------------------------------------------------
!allocate the arrays in free projection measure subroutine
!---------------------------------------------------------
subroutine allocate_meas_step()
use mc_loop_param
use method_param
use lattice_param
use meas_param
implicit none
  allocate(ZetaN_l(Nsamples,Dimen,max_local))
  allocate(kin_l(Nsamples,max_local))
  allocate(v_l(Nsamples,max_local))
  allocate(e_l(Nsamples,max_local))
  allocate(eE_l(Nsamples,max_local))
  allocate(var_l(Nsamples,max_local))
  allocate(nu_l(Nsamples,max_local))
  allocate(nd_l(Nsamples,max_local))
  if(I_obdm.eq.1)allocate(obdm_l(Nsamples,2*Nsite,2*Nsite,max_local))
!  allocate(cacb_l(Nsamples,Nbravais,Nbands,Nbands,max_local))
  allocate(nofr_l(Nsamples,Nsite,2,2,max_local))
  allocate(sisj_l(Nsamples,Nbravais,max_local))
  allocate(sk_l(Nsamples,Nsite,max_local))
!  allocate(cicj_l(Nsamples,2*Nsite,max_local))
!  allocate(ck_l(Nsamples,2*Nsite,max_local))
  if(Npair.gt.0)then
   allocate(didj_l(Nsamples,Nbravais,Npair,max_local))
  endif
  if(Nbeta.gt.0)then
   allocate(cicj_t_l(Nsamples,2*Nsite,Nbeta,max_local))
   allocate(cicjh_t_l(Nsamples,2*Nsite,Nbeta,max_local))
   allocate(GreenP_t_l(Nsamples,0:Nbeta,max_local))
!   allocate(rho_t_l(Nsamples,Nbravais,0:Nbeta,max_local))
   allocate(nupnup_t_l(Nsamples,Nbravais,0:Nbeta,max_local))
   allocate(ndnnup_t_l(Nsamples,Nbravais,0:Nbeta,max_local))
  endif
  allocate(sig(max_local))
  allocate(absig(max_local))
end subroutine allocate_meas_step


!-------------------------------------------
!deallocate the arrays in measure subroutine
!-------------------------------------------
subroutine deallocate_meas()
use meas_param
implicit none
if(allocated(ZetaN_l)) deallocate(ZetaN_l)
if(allocated(kin_l)) deallocate(kin_l)
if(allocated(v_l)) deallocate(v_l)
if(allocated(e_l)) deallocate(e_l)
if(allocated(eE_l)) deallocate(eE_l)
if(allocated(var_l)) deallocate(var_l)
if(allocated(nu_l)) deallocate(nu_l)
if(allocated(nd_l)) deallocate(nd_l)
if(allocated(obdm_l)) deallocate(obdm_l)
if(allocated(cacb_l)) deallocate(cacb_l)
if(allocated(nofr_l)) deallocate(nofr_l)
if(allocated(sisj_l)) deallocate(sisj_l)
if(allocated(ninj_l)) deallocate(ninj_l)
if(allocated(szsz_l)) deallocate(szsz_l)
if(allocated(sk_l)) deallocate(sk_l)
if(allocated(cicj_l)) deallocate(cicj_l)
if(allocated(ck_l)) deallocate(ck_l)
if(allocated(cicj_t_l)) deallocate(cicj_t_l)
if(allocated(cicjh_t_l)) deallocate(cicjh_t_l)
if(allocated(GreenP_t_l)) deallocate(GreenP_t_l)
!if(allocated(rho_t_l)) deallocate(rho_t_l)
if(allocated(nupnup_t_l)) deallocate(nupnup_t_l)
if(allocated(ndnnup_t_l)) deallocate(ndnnup_t_l)
if(allocated(sig)) deallocate(sig)
if(allocated(absig)) deallocate(absig)
end subroutine deallocate_meas


!--------------------------------------------------
!This subroutine allocate the arrays in one measure
!--------------------------------------------------
subroutine allocate_one_meas()
use one_meas_param
use lattice_param
use method_param
implicit none

allocate(ZetaN_one(Dimen))
!allocate(Cab_one(Nbravais,Nbands,Nbands))
if(I_obdm.eq.1)allocate(Obdm_one(2*Nsite,2*Nsite))
allocate(Nr_one(Nsite,2,2))
allocate(S_one(Nbravais))
allocate(N_one(Nbravais,Nbands,Nbands))
allocate(Sz_one(Nbravais,Nbands,Nbands))
allocate(c_one(2*Nsite))
if(Npair.gt.0)then
  allocate(D_one(Nbravais,Npair))
endif
if(Nbeta.gt.0)then
  if(I_onebf.eq.1)then
    allocate(GP_one(2*Nsite,Nbeta))
    allocate(GH_one(2*Nsite,Nbeta))
  endif
!  allocate(RHO_one(Nbravais,0:Nbeta))
  if(I_twob.eq.1)then
    allocate(nupnup_one(Nbravais,0:Nbeta))
    allocate(ndnnup_one(Nbravais,0:Nbeta))
  endif
  allocate(mu0(2*Nsite),mu(2*Nsite))
  allocate(GreenP_one(0:Nbeta))
  allocate(GreenH_one(0:Nbeta))
endif
end subroutine allocate_one_meas



!----------------------------------------------------
!This subroutine deallocate the arrays in one measure
!----------------------------------------------------
subroutine deallocate_one_meas()
use one_meas_param
use method_param
implicit none
if(allocated(ZetaN_one)) deallocate(ZetaN_one)
if(allocated(Cab_one)) deallocate(Cab_one)
if(allocated(Obdm_one)) deallocate(Obdm_one)
if(allocated(Nr_one)) deallocate(Nr_one)
if(allocated(S_one)) deallocate(S_one)
if(allocated(N_one)) deallocate(N_one)
if(allocated(Sz_one)) deallocate(Sz_one)
if(allocated(c_one)) deallocate(c_one)
if(allocated(D_one)) deallocate(D_one)
if(allocated(GP_one)) deallocate(GP_one)
if(allocated(GH_one)) deallocate(GH_one)
!if(allocated(RHO_one)) deallocate(RHO_one)
if(allocated(nupnup_one)) deallocate(nupnup_one)
if(allocated(ndnnup_one)) deallocate(ndnnup_one)
if(allocated(GreenP_one)) deallocate(GreenP_one)
if(allocated(GreenH_one)) deallocate(GreenH_one)
end subroutine deallocate_one_meas


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
use timing_module
use rand_num
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
call end_genrand()
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
