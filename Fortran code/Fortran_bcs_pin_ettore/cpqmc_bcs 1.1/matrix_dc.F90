subroutine copy_G_dc(A,B)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::A(2*Nsite,2*Nsite)
complex(kind=8),intent(OUT)::B(2*Nsite,2*Nsite)

if(dtype.EQ.'c') then
  call zcopy(2*Nsite*2*Nsite,A(1,1),1,B(1,1),1)
elseif(dtype.EQ.'d') then
  B(1:Nsite,1:Nsite)=A(1:Nsite,1:Nsite)
  B((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=A((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))
endif

end subroutine copy_G_dc


subroutine copy_B_dc(A,B)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::A(2*Nsite,2*Nsite)
complex(kind=8),intent(OUT)::B(2*Nsite,2*Nsite)

  B(1:Nsite,1:Nsite)=A(1:Nsite,1:Nsite)
  B((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=A((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))

end subroutine copy_B_dc



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


!Get the important function bcs unpaired case
subroutine unp_bcs_over_lap_dc(F,D,ph,ovp)
use param
use lattice_param
use model_param
use caldet_module
implicit none
complex(kind=8),intent(IN)::F(Nsite,Nsite),D(Nsite,Nzeta)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::ovp(Nspin(1),Nspin(1))
complex(kind=8)::phu(Nsite,Nspin(1)),phd(Nsite,Nspin(2))
complex(kind=8)::Fstarphd(Nsite,Nspin(2)),DstarFstarphd(Nsite,Nspin(1))
integer::i_deb,j_deb


phu(:,:)=ph(1:Nsite,1:Nspin(1))
phd(:,:)=ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot)

call zgemm('N','N',Nsite,Nspin(2),Nsite,one,conjg(F(:,:)),Nsite,phd,Nsite,zero,Fstarphd,Nsite)
DstarFstarphd(1:Nsite,1:Nzeta)=conjg(D(1:Nsite,1:Nzeta))
DstarFstarphd(1:Nsite,Nzeta+1:Nspin(1))=Fstarphd(1:Nsite,1:Nspin(2))
call zgemm('T','N',Nspin(1),Nspin(1),Nsite,one,phu,Nsite,DstarFstarphd,Nsite,zero,ovp,Nspin(1))


!DEBUG
!write(*,*)
!write(*,*)'In unp_bcs_over_lap_dc '
!write(*,*)
!write(*,*)'phu '
!do i_deb=1,Nsite
!  write(*,*)(dble(phu(i_deb,j_deb)),j_deb=1,Nspin(1))
!enddo
!write(*,*)
!write(*,*)'phd '
!do i_deb=1,Nsite
!  write(*,*)(dble(phd(i_deb,j_deb)),j_deb=1,Nspin(2))
!enddo
!write(*,*)
!write(*,*)'Dstar '
!do i_deb=1,Nsite
!  write(*,*)(dble(conjg(D(i_deb,j_deb))),j_deb=1,Nzeta)
!enddo
!write(*,*)
!write(*,*)'Fstar Phidn '
!do i_deb=1,Nsite
!  write(*,*)(dble(Fstarphd(i_deb,j_deb)),j_deb=1,Nspin(2))
!enddo
!write(*,*)
!write(*,*)
!write(*,*)'DstarFstarphd '
!do i_deb=1,Nsite
!  write(*,*)(dble(DstarFstarphd(i_deb,j_deb)),j_deb=1,Nspin(1))
!enddo
!write(*,*)

!write(*,*)'ovp '
!do i_deb=1,Nspin(1)
!  write(*,*)(dble(ovp(i_deb,j_deb)),j_deb=1,Nspin(1))
!enddo
!write(*,*)
!stop 'deb'



end subroutine unp_bcs_over_lap_dc



!Get the <BCS|ph>=ovp
subroutine bcs_over_lap_dc(F,ph,ovp)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::F(Nsite,Nsite)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::ovp(Nspin(1),Nspin(1))
complex(kind=8)::dummy(Nspin(1),Nsite)
!dummy = phi_up^ T * F^conjg
! ovp = dummy * phi_dn
dummy=0
ovp=0
call zgemm('T','N',Nspin(1),Nsite,Nsite,one,ph(1:Nsite,1:Nspin(1)),Nsite,conjg(F(:,:)),Nsite,zero,dummy,Nspin(1))
call zgemm('N','N',Nspin(1),Nspin(1),Nsite,one,dummy,Nspin(1),ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,ovp,Nspin(1))
end subroutine bcs_over_lap_dc



!get the <pht|ph>=ovp 
subroutine over_lap_dc(pht,ph,ovp)
use param
use lattice_param
use model_param
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




!Get the important function bcs unpaired case
subroutine unp_bcs_imp_fun_dc(F,D,ph,imp)
use param
use lattice_param
use model_param
use caldet_module
implicit none
complex(kind=8),intent(IN)::F(Nsite,Nsite),D(Nsite,Nzeta)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::imp
complex(kind=8)::ovp(Nspin(1),Nspin(1))
complex(kind=8)::phu(Nsite,Nspin(1)),phd(Nsite,Nspin(2))
complex(kind=8)::Fstarphd(Nsite,Nspin(2)),DstarFstarphd(Nsite,Nspin(1))
integer::i_deb,j_deb


phu(:,:)=ph(1:Nsite,1:Nspin(1))
phd(:,:)=ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot)

call zgemm('N','N',Nsite,Nspin(2),Nsite,one,conjg(F(:,:)),Nsite,phd,Nsite,zero,Fstarphd,Nsite)
DstarFstarphd(1:Nsite,1:Nzeta)=conjg(D(1:Nsite,1:Nzeta))
DstarFstarphd(1:Nsite,Nzeta+1:Nspin(1))=Fstarphd(1:Nsite,1:Nspin(2))
call zgemm('T','N',Nspin(1),Nspin(1),Nsite,one,phu,Nsite,DstarFstarphd,Nsite,zero,ovp,Nspin(1))
call caldet(Nspin(1),ovp(1:Nspin(1),1:Nspin(1)),imp)

end subroutine unp_bcs_imp_fun_dc



!Get the important function bcs case
subroutine bcs_imp_fun_dc(F,ph,imp)
use param
use lattice_param
use model_param
use caldet_module
implicit none
complex(kind=8),intent(IN)::F(Nsite,Nsite)
complex(kind=8),intent(IN)::ph(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::imp
complex(kind=8)::dummy(Nspin(1),Nsite),ovp(Nspin(1),Nspin(1))
!complex(kind=8)::phu(Nsite,Nspin(1)),phd(Nsite,Nspin(2))
integer::i_deb,j_deb
!dummy = phi_up^ T * F^conjg
! ovp = dummy * phi_dn

!phu(:,:)=ph(1:Nsite,1:Nspin(1))
!phd(:,:)=ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot)


!DEB
!write(*,*)'phi up '
!do i_deb=1,Nsite
!  write(*,*)(phu(i_deb,j_deb),j_deb=1,Nspin(1))
!enddo
!write(*,*)
!write(*,*)'phi dn '
!do i_deb=1,Nsite
!  write(*,*)(ph(i_deb,j_deb),j_deb=1,Nspin(2))
!enddo
!write(*,*)
!write(*,*)'F '
!do i_deb=1,Nsite
!  write(*,*)(F(i_deb,j_deb),j_deb=1,Nsite)
!enddo

!write(*,*)'m,n,k,lda,ldb,ldc'
!write(*,*)Nspin(1),Nsite,Nsite,Nsite,Nsite,Nsite
!write(*,*)'alpha,beta'
!write(*,*)one,zero

call zgemm('T','N',Nspin(1),Nsite,Nsite,one,ph(1:Nsite,1:Nspin(1)),Nsite,conjg(F(:,:)),Nsite,zero,dummy,Nspin(1))
call zgemm('N','N',Nspin(1),Nspin(1),Nsite,one,dummy,Nspin(1),ph(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,ovp,Nspin(1))


!write(*,*)'  phi_up^ T * F^conjg '
!do i_deb=1,Nspin(1)
!  write(*,*)(dummy(i_deb,j_deb),j_deb=1,Nsite)
!enddo
!write(*,*)' A matrix '
!do i_deb=1,Nspin(1)
!  write(*,*)(ovp(i_deb,j_deb),j_deb=1,Nspin(1))
!enddo


call caldet(Nspin(1),ovp(1:Nspin(1),1:Nspin(1)),imp)  !calculate the determinant of a matrix


!write(*,*)'Importance function = ',imp

end subroutine bcs_imp_fun_dc


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

!This subroutine calculate <phi_l|ci cj^+|phi_r>/<phi_l|phi_r> matrix
subroutine cal_Amat_withovlpinv_dc2(phi_l,phi_r,ovlpinv,Amat)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN):: phi_l(2*Nsite,Ntot),phi_r(2*Nsite,Ntot),ovlpinv(Ntot,Ntot)
complex(kind=8),intent(OUT):: Amat(2*Nsite,2*Nsite)

if(dtype.EQ.'c') then
  call cal_Amat_withovlpinv2(2*Nsite,Ntot,phi_l,phi_r,ovlpinv,Amat)
else if(dtype.EQ.'d') then
  !test
  !Amat=zero
  !end test
  call cal_Amat_withovlpinv2(Nsite,Nspin(1),phi_l(1:Nsite,1:Nspin(1)),phi_r(1:Nsite,1:Nspin(1)),&
       & ovlpinv(1:Nspin(1),1:Nspin(1)),Amat(1:Nsite,1:Nsite))
  call cal_Amat_withovlpinv2(Nsite,Nspin(2),phi_l((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot), &
       & phi_r((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot), &
       & ovlpinv((Nspin(1)+1):Ntot,(Nspin(1)+1):Ntot),Amat((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite)))
end if



end subroutine cal_Amat_withovlpinv_dc2


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


!This subroutine calculate <BCS|ci^+ cj|phi_r>/<BCS|phi_r> matrix
subroutine bcs_cal_Amat_withovlpinv_dc(F,phi_r,ovlpinv,Amat,Bmat)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::F(Nsite,Nsite),phi_r(2*Nsite,Ntot),ovlpinv(Ntot,Ntot)
complex(kind=8),intent(OUT):: Amat(2*Nsite,2*Nsite),Bmat(2*Nsite,2*Nsite)
complex(kind=8)::FT(Nsite,Nsite)
complex(kind=8)::dummy1(Nsite,Nspin(1)),dummy2(Nspin(1),Nsite)
complex(kind=8)::dummy3(Nsite,Nspin(1)),dummy4(Nsite,Nspin(1))
integer::i,j

!dummy1 = F* * phi_dn
!dummy2 = A^{-1} * phi_up^ T 
!Amat_upup = dummy1 * dummy2
call zgemm('N','N',Nsite,   Nspin(1),Nsite,   one,conjg(F(:,:)),Nsite                        &
     &      ,phi_r(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,dummy1,Nsite)
call zgemm('N','T',Nspin(1),Nsite,   Nspin(1),one,ovlpinv(1:Nspin(1),1:Nspin(1)),Nspin(1)    &
     &      ,phi_r(1:Nsite,1:Nspin(1)),Nsite,zero,dummy2,Nspin(1))
call zgemm('N','N',Nsite   ,Nsite,   Nspin(1),one,dummy1,Nsite                               &
     &      ,dummy2,Nspin(1),zero,Amat(1:Nsite,1:Nsite),Nsite)

!dummy3 = F^{dagger} * phi_up
!dummy4 = dummy3 * A^{-1}
!Amat_dndn = dummy4 * phi_dn^T
FT = transpose(F)
call zgemm('N','N',Nsite, Nspin(1),Nsite,   one,conjg(FT(:,:)),Nsite                         &
     &      ,phi_r(1:Nsite,1:Nspin(1)),Nsite,zero,dummy3,Nsite)
call zgemm('N','T',Nsite, Nspin(1),Nspin(1),one,dummy3,Nsite                                 &
     &      ,ovlpinv(1:Nspin(1),1:Nspin(1)),Nspin(1),zero,dummy4,Nsite)
call zgemm('N','T',Nsite, Nsite,   Nspin(1),one,dummy4,Nsite                                 &
     &      ,phi_r(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,Amat(Nsite+1:2*Nsite,Nsite+1:2*Nsite),Nsite)

!dummy2^T =  phi_up * (A^{-1})^T
!Bmat_upup = dummy2^T * phi_dn^T (it is the anomalous <cc>)
call zgemm('T','T',Nsite, Nsite,Nspin(2),one,dummy2,Nspin(1)                                   &
     &      ,phi_r(Nsite+1:2*Nsite,Nspin(1)+1:Ntot),Nsite,zero,Bmat(1:Nsite,1:Nsite),Nsite)

!dummy1 = F* * phi_dn
!dummy3^T = (F^{dagger} * phi_up)^T
!(new) dummy2 = A^{-1} * (F^{dagger} * phi_up)^T
!Bmat_dndn = dummy1 * dummy2 - F   (it is the anomalous <c+c+>)
call zgemm('N','T',Nspin(1),Nsite,Nspin(1),one,ovlpinv(1:Nspin(1),1:Nspin(1)),Nspin(1)       &
     &      ,dummy3,Nsite,zero,dummy2,Nspin(1))
call zgemm('N','N',Nsite,Nsite,Nspin(1),one,dummy1,Nsite                                     &
     &      ,dummy2,Nspin(1),zero,Bmat(Nsite+1:2*Nsite,Nsite+1:2*Nsite),Nsite)

do i=1,Nsite
  do j=1,Nsite
    Bmat(Nsite+j,Nsite+i)=Bmat(Nsite+j,Nsite+i)-conjg(F(j,i))
  enddo
enddo

end subroutine bcs_cal_Amat_withovlpinv_dc


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
