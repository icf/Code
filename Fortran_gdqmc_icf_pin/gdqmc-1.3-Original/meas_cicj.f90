!--------------------------------
!measure the cicj_local from Amat
!--------------------------------
subroutine measure_cicj(Amat,cicj_local)
use param
use lattice_param
use phiT_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::cicj_local(2*Nsite,Dtot)
integer,external::latt_label
integer,external::bound
integer::i,j,k,m,n
integer::cc(1:Dimen),ctmp

cicj_local=zero
do k=1,Dtot,1

   !average cicj to cicj_local
   do i=1,Nsite,1
      do j=1,Nsite,1
         
         do m=1,Dimen,1
            !We calculate (Ci^+)Cj==>C(i-j+1)^+ C1
            !so we need to focus on i-j+1
            ctmp=coor(i,m)-coor(j,m)+1
            cc(m)=bound(ctmp,Nl(m))
         end do
         n=latt_label(cc(1:Dimen))
         cicj_local(n,k)=cicj_local(n,k)+Amat(i,j,k)
         cicj_local(n+Nsite,k)=cicj_local(n+Nsite,k)+Amat(i+Nsite,j+Nsite,k)
      end do
   end do
   do n=1,2*Nsite,1
      cicj_local(n,k)=cicj_local(n,k)/dcmplx(Nsite)
   end do

end do
end subroutine measure_cicj




!----------------------------------------------------------------
!fourier cicj_l(1:Nsample,1:Nsite,i) to ck_l(1:Nsample,1:Nsite,i)
!----------------------------------------------------------------
subroutine fourier_cij(i)
use param
use lattice_param
use mc_loop_param
use meas_param
implicit none
integer,intent(IN)::i
integer::j,k,m,n
integer::kl(Dimen)
complex(kind=8)::kkup,kkdn
real(kind=8)::ph



do j=1,Nsamples,1

   do k=1,Nsite,1
      !The number of kl for k
      do n=1,Dimen,1
         kl(n)=coor(k,n)-1
      end do

      !fourier of the momentum k
      kkup=zero;kkdn=zero
      do m=1,Nsite,1
         ph=0.d0
         do n=1,Dimen,1
            ph=ph+2.d0*Pi*dble(kl(n)*(coor(m,n)-1))/dble(Nl(n))
         end do
         kkup=kkup+exp(Xi*ph)*cicj_l(j,m,i)
         kkdn=kkdn+exp(Xi*ph)*cicj_l(j,m+Nsite,i)
        !write(*,*) exp(Xi*ph);pause
      end do
      ck_l(j,k,i)=dble(kkup)
      ck_l(j,k+Nsite,i)=dble(kkdn)
   end do

end do
end subroutine fourier_cij
