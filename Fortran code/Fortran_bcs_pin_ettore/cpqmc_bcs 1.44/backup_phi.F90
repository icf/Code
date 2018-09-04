!----------------------------------------------------
!This subroutine backup the wave function information 
!as well as the i_pop and i_GS
!----------------------------------------------------
subroutine back_phi(i_pop,i_GS)
use backup_param
use lattice_param
use model_param
use mc_loop_param
use phiT_param
use phi_param
use method_param
implicit none
integer,intent(IN)::i_pop,i_GS
call allocate_backup()

i_pop_tmp=i_pop;i_GS_tmp=i_GS

crn_tmp=crn

!call zcopy(2*Nsite*Ntot*Nwalkers,phi(1,1,1),1,phi_tmp(1,1,1),1)
call copy_wf_W_dc(phi(1,1,1),phi_tmp(1,1,1))
call dcopy(Nwalkers,weight(1),1,weight_tmp(1),1)
call dcopy(Nwalkers,dlogw(1),1,dlogw_tmp(1),1)
call zcopy(Nwalkers,rx(1),1,rx_tmp(1),1)
call zcopy(Dtot*Nwalkers,impfunc(1,1),1,impfunc_tmp(1,1),1)
call zcopy(Nwalkers,tot_imp(1),1,tot_imp_tmp(1),1)
end subroutine back_phi



!-----------------------------------------------------
!This subroutine recover the wave function information 
!as well as the i_pop and i_GS
!-----------------------------------------------------
subroutine recov_phi(i_pop,i_GS)
use backup_param
use lattice_param
use model_param
use mc_loop_param
use phiT_param
use phi_param
use method_param
implicit none
integer,intent(OUT)::i_pop,i_GS

i_pop=i_pop_tmp;i_GS=i_GS_tmp

crn=crn_tmp

call copy_wf_W_dc(phi_tmp(1,1,1),phi(1,1,1))
call dcopy(Nwalkers,weight_tmp(1),1,weight(1),1)
call dcopy(Nwalkers,dlogw_tmp(1),1,dlogw(1),1)
call zcopy(Nwalkers,rx_tmp(1),1,rx(1),1)
call zcopy(Dtot*Nwalkers,impfunc_tmp(1,1),1,impfunc(1,1),1)
call zcopy(Nwalkers,tot_imp_tmp(1),1,tot_imp(1),1)

call deallocate_backup()
end subroutine recov_phi
