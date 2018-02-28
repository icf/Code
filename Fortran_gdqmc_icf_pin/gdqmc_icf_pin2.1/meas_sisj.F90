!--------------------------------
!measure the sisj_local from Amat
!--------------------------------
subroutine measure_sisj(Amat,sisj_local)
use param
use lattice_param
use phiT_param
use model_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::sisj_local(Nsite,Dtot)
integer,external::latt_label
integer,external::bound
complex(kind=8)::sisj(Nsite,Nsite)
integer::i,j,k,m,n
integer::cc(1:Dimen),ctmp

sisj_local=zero
do k=1,Dtot,1
   !Get sisj for k
   do i=1,Nsite,1
      sisj_local(i,k)=sisj_local(i,k)+(Amat(i,i,k)-Amat(i+Nsite,i+Nsite,k))/2.d0
   end do

end do
end subroutine measure_sisj




!----------------------------------------------------------------
!fourier sisj_l(1:Nsample,1:Nsite,i) to sk_l(1:Nsample,1:Nsite,i)
!----------------------------------------------------------------
subroutine fourier_sij(i)
use param
use lattice_param
use mc_loop_param
use meas_param
implicit none
integer,intent(IN)::i
integer::j,k,m,n
integer::kl(Dimen)
complex(kind=8)::kk
real(kind=8)::ph



do j=1,Nsamples,1

   do k=1,Nsite,1
      !The number of kl for k
      do n=1,Dimen,1
         kl(n)=coor(k,n)-1
      end do

      !fourier of the momentum k
      kk=zero
      do m=1,Nsite,1
         ph=0.d0
         do n=1,Dimen,1
            ph=ph+2.d0*Pi*dble(kl(n)*(coor(m,n)-1))/dble(Nl(n))
         end do
         kk=kk+exp(Xi*ph)*sisj_l(j,m,i)
        !write(*,*) exp(Xi*ph);pause
      end do
      sk_l(j,k,i)=dble(kk)
   end do

end do
end subroutine fourier_sij
