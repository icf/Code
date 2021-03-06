!copy ph(:,:) to phn(:,:)
subroutine copy_wf_dc(ph,phn)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::phn(2*Nsite,Ntot)

if(dtype.EQ.'c') then
  call zcopy(2*Nsite*Ntot,ph(1,1),1,phn(1,1),1)
else if(dtype.EQ.'d') then
  phn(1:Nsite,1:Nspin(1))=ph(1:Nsite,1:Nspin(1))
  phn((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot)=ph((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot)
end if
end subroutine copy_wf_dc



!copy ph(:,:,:) to phn(:,:,:) of phiT
subroutine copy_wf_T_dc(ph,phn)
use param
use lattice_param
use model_param
use phiT_param
implicit none
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot,Dtot)
complex(kind=8),intent(OUT)::phn(2*Nsite,Ntot,Dtot)

if(dtype.EQ.'c') then
  call zcopy(2*Nsite*Ntot*Dtot,ph(1,1,1),1,phn(1,1,1),1)
else if(dtype.EQ.'d') then
  phn(1:Nsite,1:Nspin(1),1:Dtot)=ph(1:Nsite,1:Nspin(1),1:Dtot)
  phn((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,1:Dtot)=ph((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,1:Dtot)
end if
end subroutine copy_wf_T_dc


!copy ph(:,:,:) to phn(:,:,:) of Nwalkers phi
subroutine copy_wf_W_dc(ph,phn)
use param
use lattice_param
use model_param
use mc_loop_param
implicit none
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot,Nwalkers)
complex(kind=8),intent(OUT)::phn(2*Nsite,Ntot,Nwalkers)

if(dtype.EQ.'c') then
  call zcopy(2*Nsite*Ntot*Nwalkers,ph(1,1,1),1,phn(1,1,1),1)
else if(dtype.EQ.'d') then
  phn(1:Nsite,1:Nspin(1),1:Nwalkers)=ph(1:Nsite,1:Nspin(1),1:Nwalkers)
  phn((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,1:Nwalkers)=ph((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,1:Nwalkers)
end if
end subroutine copy_wf_W_dc


!get the <pht|ph>=ovp 
subroutine over_lap_dc(pht,ph,ovp)
use param
use lattice_param
use model_param

use sc_loop_param
implicit none
complex(kind=8),intent(IN)::pht(2*Nsite,Ntot)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::ovp(Ntot,Ntot)
if(dtype.EQ.'c') then
  call deter_overlap(2*Nsite,Ntot,pht(1,1),ph(1,1),ovp(1,1))
else if(dtype.EQ.'d') then   
  call zgemm('C','N',Nspin(1),Nspin(1),Nsite,one,pht(1,1),2*Nsite,ph(1,1),2*Nsite,zero,ovp(1,1),Ntot)
  call zgemm('C','N',Nspin(2),Nspin(2),Nsite,one,pht(Nsite+1,Nspin(1)+1),2*Nsite,ph(Nsite+1,Nspin(1)+1), &
           & 2*Nsite,zero,ovp((Nspin(1)+1),(Nspin(1)+1)),Ntot)
end if
end subroutine over_lap_dc


!get the <BCS|ph>=ovp 
subroutine over_lap_BCS_dc(ph,ovp)
use param
use lattice_param
use model_param

use sc_loop_param
use mpi_serial_param
implicit none
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::ovp(Nspin(1),Nspin(1))
complex(kind=8)::temp_BCS_sc(Nspin(1),Nsite)
if(dtype.EQ.'c') then
  write(*,*)'no c in BCS'
  call mystop
else if(dtype.EQ.'d') then   
  !ovp=matmul(matmul(transpose(ph(1:Nsite,1:Nspin(1))),BCS_sc),ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot))
temp_BCS_sc=0
ovp=0
call zgemm('T','N',Nspin(1),Nsite,Nsite,one,ph(1:Nsite,1:Nspin(1)),Nsite,BCS_sc,Nsite,zero,temp_BCS_sc,Nspin(1))
call zgemm('N','N',Nspin(1),Nspin(1),Nsite,one,temp_BCS_sc,Nspin(1),ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,ovp,Nspin(1))
end if
end subroutine over_lap_BCS_dc
    

!calculate the determinate of a matrix
subroutine caldet_dc(ovp,imp)
use param
use lattice_param
use model_param
use caldet_module
implicit none
complex(kind=8),intent(IN)::ovp(Ntot,Ntot)
complex(kind=8),intent(OUT)::imp
complex(kind=8)::imp1,imp2
if(dtype.EQ.'c') then
  call caldet(Ntot,ovp(1:Ntot,1:Ntot),imp)
else if(dtype.EQ.'d') then
  call caldet(Nspin(1),ovp(1:Nspin(1),1:Nspin(1)),imp1)
  call caldet(Nspin(2),ovp((Nspin(1)+1):Ntot,(Nspin(1)+1):Ntot),imp2)
  imp=imp1*imp2
end if
end subroutine caldet_dc



!get the Inverse of ovp 
subroutine inverse_dc(ovp)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(INOUT)::ovp(Ntot,Ntot)

if(dtype.EQ.'c') then
  call inverse(ovp(1:Ntot,1:Ntot),Ntot)
else if(dtype.EQ.'d') then
  call inverse(ovp(1:Nspin(1),1:Nspin(1)),Nspin(1))
  call inverse(ovp((Nspin(1)+1):Ntot,(Nspin(1)+1):Ntot),Nspin(2))
end if
end subroutine inverse_dc


!get the Inverse of ovp with the determinant 
subroutine inverse_d_dc(ovp,det)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(INOUT)::ovp(Ntot,Ntot)
complex(kind=8),intent(OUT)::det
complex(kind=8)::det1,det2

if(dtype.EQ.'c') then
  call inverse_d(ovp(1:Ntot,1:Ntot),Ntot,det)
else if(dtype.EQ.'d') then
  call inverse_d(ovp(1:Nspin(1),1:Nspin(1)),Nspin(1),det1)
  call inverse_d(ovp((Nspin(1)+1):Ntot,(Nspin(1)+1):Ntot),Nspin(2),det2)
  det=det1*det2
end if
end subroutine inverse_d_dc
    

!get the Inverse[<pht|ph>]=ovp 
subroutine over_lap_inv_dc(pht,ph,ovp)
use param
use lattice_param
use model_param

implicit none
complex(kind=8),intent(IN)::pht(2*Nsite,Ntot)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::ovp(Ntot,Ntot)

call over_lap_dc(pht,ph,ovp)
call inverse_dc(ovp)
end subroutine over_lap_inv_dc



!Get the important function
subroutine imp_fun_BCS_dc(ph,imp)
use param
use lattice_param
use model_param
use caldet_module

use sc_loop_param
use mpi_serial_param
implicit none
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::imp
complex(kind=8)::ovp(Nspin(1),Nspin(1))

call over_lap_BCS_dc(ph(1,1),ovp(1,1))
call caldet(sc_loop_Na,ovp(1:sc_loop_Na,1:sc_loop_Na),imp)

end subroutine imp_fun_BCS_dc


    
!Get the important function
subroutine imp_fun_dc(pht,ph,imp)
use param
use lattice_param
use model_param
use caldet_module
implicit none
complex(kind=8),intent(IN)::pht(2*Nsite,Ntot)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::imp
complex(kind=8)::ovp(Ntot,Ntot)

call over_lap_dc(pht(1,1),ph(1,1),ovp(1,1))
call caldet_dc(ovp,imp)
end subroutine imp_fun_dc
    

!Get the exp(K)|ph>=|phn>
subroutine k_to_ph_dc(expk,ph,phn)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::expk(2*Nsite,2*Nsite)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::phn(2*Nsite,Ntot)

if(dtype.EQ.'c') then
  call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,expk,2*Nsite,ph(1,1),2*Nsite,zero,phn(1,1),2*Nsite)
else if(dtype.EQ.'d') then
  call zgemm('N','N',Nsite,Nspin(1),Nsite,one,expk(1,1),2*Nsite,ph(1,1),2*Nsite,zero,phn(1,1),2*Nsite)
  call zgemm('N','N',Nsite,Nspin(2),Nsite,one,expk(Nsite+1,Nsite+1),2*Nsite,ph(Nsite+1,Nspin(1)+1),2*Nsite, &
           & zero,phn(Nsite+1,Nspin(1)+1),2*Nsite)
end if
end subroutine k_to_ph_dc


!Get the dagger{exp(K)}|ph>=|phn>
subroutine dk_to_ph_dc(expk,ph,phn)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::expk(2*Nsite,2*Nsite)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::phn(2*Nsite,Ntot)

if(dtype.EQ.'c') then
  call zgemm('C','N',2*Nsite,Ntot,2*Nsite,one,expk,2*Nsite,ph(1,1),2*Nsite,zero,phn(1,1),2*Nsite)
else if(dtype.EQ.'d') then
  call zgemm('C','N',Nsite,Nspin(1),Nsite,one,expk(1,1),2*Nsite,ph(1,1),2*Nsite,zero,phn(1,1),2*Nsite)
  call zgemm('C','N',Nsite,Nspin(2),Nsite,one,expk(Nsite+1,Nsite+1),2*Nsite,ph(Nsite+1,Nspin(1)+1),2*Nsite, &
           & zero,phn(Nsite+1,Nspin(1)+1),2*Nsite)
end if
end subroutine dk_to_ph_dc


!This subroutine calculate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r> matrix
subroutine cal_Amat_withovlpinv_dc(phi_l,phi_r,ovlpinv,Amat)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN):: phi_l(2*Nsite,Ntot),phi_r(2*Nsite,Ntot),ovlpinv(Ntot,Ntot)
complex(kind=8),intent(OUT):: Amat(2*Nsite,2*Nsite)

if(dtype.EQ.'c') then
  call cal_Amat_withovlpinv(2*Nsite,Ntot,phi_l,phi_r,ovlpinv,Amat) 
else if(dtype.EQ.'d') then
  !test
  !Amat=zero
  !end test
  call cal_Amat_withovlpinv(Nsite,Nspin(1),phi_l(1:Nsite,1:Nspin(1)),phi_r(1:Nsite,1:Nspin(1)), &
       & ovlpinv(1:Nspin(1),1:Nspin(1)),Amat(1:Nsite,1:Nsite))
  call cal_Amat_withovlpinv(Nsite,Nspin(2),phi_l((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot), &
       & phi_r((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot), &
       & ovlpinv((Nspin(1)+1):Ntot,(Nspin(1)+1):Ntot),Amat((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite)))
end if



end subroutine cal_Amat_withovlpinv_dc   
 
    
    
subroutine cal_BCS_Amat_withovlpinv_dc(phR,inver,GF_up,GF_dn)
use param
use model_param
use phiT_param

use lattice_param
use phi_param

use sc_loop_param
implicit none
complex(kind=8),intent(IN)::inver(Nspin(1),Nspin(1)),phR(2*Nsite,Ntot)
complex(kind=8),external::zdotu
complex(kind=8)::GF_up(Nsite,Nsite),GF_dn(Nsite,Nsite)
complex(kind=8)::temp_BCS_sc1(Nsite,Nspin(1)),temp_BCS_sc2(Nspin(1),Nsite)
integer::i,j
if(dtype.EQ.'c') then
  write(*,*)'no couple in BCS_sc now'
  call mystop
else if(dtype.EQ.'d') then 
     !--------------------------------
     temp_BCS_sc1=0
     temp_BCS_sc2=0 
     call zgemm('N','N',Nsite,Nspin(1),Nsite,one,BCS_sc,Nsite,&
phR(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,temp_BCS_sc1,Nsite)
     call zgemm('N','T',Nspin(1),Nsite,Nspin(1),one,inver(1:Nspin(1),1:Nspin(1)),Nspin(1),&
phR(1:Nsite,1:Nspin(1)),Nsite,zero,temp_BCS_sc2,Nspin(1))
     call zgemm('N','N',Nsite,Nsite,Nspin(1),one,temp_BCS_sc1,Nsite,temp_BCS_sc2,Nspin(1),zero,GF_up,Nsite)

     !may have some problem!!!!!!!!!!
     temp_BCS_sc1=0
     temp_BCS_sc2=0 
     call zgemm('N','N',Nsite,Nspin(1),Nsite,one,transpose(BCS_sc),Nsite,&
phR(1:Nsite,1:Nspin(1)),Nsite,zero,temp_BCS_sc1,Nsite)
     call zgemm('N','T',Nspin(1),Nsite,Nspin(1),one,transpose(inver(1:Nspin(1),1:Nspin(1))),Nspin(1),&
phR(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,temp_BCS_sc2,Nspin(1))
     call zgemm('N','N',Nsite,Nsite,Nspin(1),one,temp_BCS_sc1,Nsite,temp_BCS_sc2,Nspin(1),zero,GF_dn,Nsite)
     !-------------------------------
end if
end subroutine cal_BCS_Amat_withovlpinv_dc      
    
    
!calculate the greens function of <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r> 
subroutine cal_cidcj_withovlpinv_dc(phi_l,phi_r,ovlpinv,Amat,i,j)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::phi_l(2*Nsite,Ntot),phi_r(2*Nsite,Ntot),ovlpinv(Ntot,Ntot)
complex(kind=8),intent(OUT):: Amat
integer,intent(IN)::i,j

if(dtype.EQ.'c') then
  call cal_cidcj_withovlpinv(2*Nsite,Ntot,phi_l,phi_r,ovlpinv,Amat,i,j)
else if(dtype.EQ.'d') then
  if(i.LE.Nsite.AND.j.LE.Nsite) then
    call cal_cidcj_withovlpinv(Nsite,Nspin(1),phi_l(1:Nsite,1:Nspin(1)),phi_r(1:Nsite,1:Nspin(1)), &
           &  ovlpinv(1:Nspin(1),1:Nspin(1)),Amat,i,j)
  else if(i.GT.Nsite.AND.j.GT.Nsite) then
    call cal_cidcj_withovlpinv(Nsite,Nspin(2),phi_l((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot), &
    & phi_r((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot),ovlpinv((Nspin(1)+1):Ntot,(Nspin(1)+1):Ntot),Amat,(i-Nsite),(j-Nsite))
  else
    Amat=zero
  end if
end if

end subroutine cal_cidcj_withovlpinv_dc
