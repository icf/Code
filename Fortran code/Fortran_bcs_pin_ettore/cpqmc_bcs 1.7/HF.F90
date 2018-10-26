!HF model from "Spin- and charge-density waves in the Hartree-Fock ground state of the &
!& two-dimensional   Hubbard model"
subroutine initial_HF_n(n_new_in)
use lattice_param
use model_param
implicit none
complex(kind=8),intent(OUT)::n_new_in(2*Nsite)
complex(kind=8)::nn
integer::iy,jx
integer::Lx,Ly

Lx=Nl(1)
Ly=Nl(2)
nn=(Nspin(1)+Nspin(2))/Nsite
n_new_in=0
do iy=1,Ly
    do jx=1,Lx
        if (mod(jx+iy-1,2) .EQ. 1)then
           n_new_in(jx+(iy-1)*Lx)=0
           n_new_in(jx+(iy-1)*Lx+Nsite)=nn
        else
           n_new_in(jx+(iy-1)*Lx+Nsite)=0
           n_new_in(jx+(iy-1)*Lx)=nn
        endif
    enddo
enddo

end subroutine initial_HF_n
    
    
subroutine get_phi_from_HF(n,phi_HF)
use lattice_param
use model_param
implicit none
real(kind=8)::U
complex(kind=8)::U2_diag
complex(kind=8)::U1_up(Nsite,Nsite),U1_dn(Nsite,Nsite)
complex(kind=8)::H_HF_up(Nsite,Nsite),H_HF_dn(Nsite,Nsite)
complex(kind=8)::D_up(Nsite),D_dn(Nsite)
complex(kind=8),intent(in)::n(2*Nsite)
complex(kind=8),intent(OUT)::phi_HF(2*Nsite,Ntot)

integer::i

U=onsitU

U2_diag=0
do i=1,Nsite
    U2_diag=U2_diag-0.5*U*n(i)*n(i+Nsite);
enddo

U1_up=0
U1_dn=0
do i=1,Nsite
    U1_up(i,i)=U*n(i+Nsite)+U2_diag/(Ntot);
    U1_dn(i,i)=U*n(i)+U2_diag/(Ntot);
enddo

H_HF_up=Hzero(1:Nsite,1:Nsite)+U1_up;
H_HF_dn=Hzero(Nsite+1:2*Nsite,Nsite+1:2*Nsite)+U1_dn;

D_up=0
D_dn=0
call eigen(H_HF_up,Nsite,D_up)
call eigen(H_HF_dn,Nsite,D_dn)

!E_nonint_up=diag(D_up);
!E_nonint_dn=diag(D_dn);
!E=sum(E_nonint_up(1:N_up))+sum(E_nonint_dn(1:N_dn));
!-------------------------------
!check if choose those min terms
!-------------------------------
phi_HF=0
phi_HF(1:Nsite,1:Nspin(1))=H_HF_up(1:Nsite,1:Nspin(1))
phi_HF(Nsite+1:2*Nsite,Nspin(1)+1:Ntot)=H_HF_dn(1:Nsite,1:Nspin(2)) 

end subroutine get_phi_from_HF
    
  
!-----------------------
subroutine get_n_from_phi(n,phi_HF)
use lattice_param
use model_param
implicit none
complex(kind=8)::inv_ovp(Ntot,Ntot)
complex(kind=8)::Green_HF(2*Nsite,2*Nsite),n_up(Nsite),n_dn(Nsite)
complex(kind=8),intent(IN)::phi_HF(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::n(2*Nsite)

integer::i

inv_ovp=0
Green_HF=0
call over_lap_inv_dc(phi_HF,phi_HF,inv_ovp)
call cal_Amat_withovlpinv_dc(phi_HF,phi_HF,inv_ovp,Green_HF)

do i=1,Nsite
   n_up(i)=abs(Green_HF(i,i))
   n_dn(i)=abs(Green_HF(Nsite+i,Nsite+i))
enddo
n(1:Nsite)=n_up(1:Nsite)
n(Nsite+1:2*Nsite)=n_dn(1:Nsite)

end subroutine get_n_from_phi
    
    
subroutine get_new_n(a,n_old_in,n_old_out,n_new_in)
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::n_old_in(2*Nsite),n_old_out(2*Nsite),a
complex(kind=8),intent(OUT)::n_new_in(2*Nsite)

n_new_in=(1-a)*n_old_in+a*n_old_out;

end subroutine get_new_n
