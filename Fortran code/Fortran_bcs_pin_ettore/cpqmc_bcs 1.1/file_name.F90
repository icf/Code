!-------------------------------------------------------
!This subroutine get the filename needs to be wrote done
!-------------------------------------------------------
subroutine get_filename()
use io_module
use lattice_param
use model_param
use method_param
use mc_loop_param
use project_param
use phiT_param
use phi_param
use adET_param
use mpi_serial_param
use meas_param
implicit none
call createFileName(basename,'hubb_gdmc_')
!call appendBaseName(basename,'se_',set)
!call appendBaseName(basename,'ns_',Nsite)
!call appendBaseName(basename,'nh_',Nhop)
!call appendBaseName(basename,'d_',Dimen)
!call appendBaseName(basename,'nl_',Nl(1))
!call appendBaseName(basename,'_',Nl(2))
!call appendBaseName(basename,'_',Nl(3))
!call appendBaseName(basename,'kb_',3,kbound(1))
!call appendBaseName(basename,'_',3,kbound(2))
!call appendBaseName(basename,'_',3,kbound(3))
!call appendBaseName(basename,'t_',3,dble(t1))
!call appendBaseName(basename,'U_',3,onsitU)
!call appendBaseName(basename,'N_',Ntot)
!call appendBaseName(basename,'_')
!call appendBaseName(basename,dtype)
!call appendBaseName(basename,'_')
!call appendBaseName(basename,'Ns1_',Nspin(1))
!call appendBaseName(basename,'Ns2_',Nspin(2))
!call appendBaseName(basename,'cr_',3,crn)
!call appendBaseName(basename,'ce_',3,ccoe)
!call appendBaseName(basename,'mn_',max_crn)
!call appendBaseName(basename,'_')
!call appendBaseName(basename,dcp)
!call appendBaseName(basename,'_')
!call appendBaseName(basename,'kc_',kcrn)
!call appendBaseName(basename,'_')
!call appendBaseName(basename,fw_bk)
!call appendBaseName(basename,'_')
!call appendBaseName(basename,'nk1_',Nbk(1))
!call appendBaseName(basename,'nk2_',Nbk(2))
!call appendBaseName(basename,'nsm_',Nsamples)
!call appendBaseName(basename,'nw_',Nwalkers)
!call appendBaseName(basename,'blp_',blockstep)
!call appendBaseName(basename,'thm_',Thermblock)
!call appendBaseName(basename,'eqb_',Neqblock)
!call appendBaseName(basename,'pop_',PopContrlstep)
!call appendBaseName(basename,'mgs_',StepforGram)
!call appendBaseName(basename,'msp_',meastep)
!call appendBaseName(basename,'dt_',3,dt)
!call appendBaseName(basename,'pt_',PT)
!call appendBaseName(basename,'pp_',PP)
!call appendBaseName(basename,'mad_',max_ad)
!call appendBaseName(basename,'size_',Nsize)

call copyName(BaseName,EnergyName)
call appendBaseName(EnergyName,'_energy.dat')

call copyName(BaseName,ZetaNName)
call appendBaseName(ZetaNName,'_zetaN.dat')

call copyName(BaseName,numName)
call appendBaseName(numName,'_num.dat')

call copyName(BaseName,cabName)
call appendBaseName(cabName,'_cacb.dat')

call copyName(BaseName,obdmName)
call appendBaseName(obdmName,'_obdm.dat')

call copyName(BaseName,nofrName)
call appendBaseName(nofrName,'_nofr.dat')

call copyName(BaseName,ScorrName)
call appendBaseName(ScorrName,'_sisj.dat')

call copyName(BaseName,NcorrName)
call appendBaseName(NcorrName,'_ninj.dat')

call copyName(BaseName,SzcorrName)
call appendBaseName(SzcorrName,'_szsz.dat')

call copyName(BaseName,SkName)
call appendBaseName(SkName,'_sk.dat')

call copyName(BaseName,ckName)
call appendBaseName(ckName,'_ck.dat')

call copyName(Basename,DcorrName)
call appendBaseName(DcorrName,'_didj.dat')

call copyName(BaseName,cijtName)
call appendBaseName(cijtName,'_cic+jt.dat')

call copyName(BaseName,cijhtName)
call appendBaseName(cijhtName,'_c+icjt.dat')

!call copyName(BaseName,rhotName)
!call appendBaseName(rhotName,'_rhot.dat')

call copyName(BaseName,nupnuptName)
call appendBaseName(nupnuptName,'_nupnupt.dat')

call copyName(BaseName,ndnnuptName)
call appendBaseName(ndnnuptName,'_ndnnupt.dat')

call copyName(BaseName,GreenPtName)
call appendBaseName(GreenPtName,'_cnuc+mut.dat')

call copyName(BaseName,GreenHtName)
call appendBaseName(GreenHtName,'_c+nucmut.dat')

end subroutine get_filename
