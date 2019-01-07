subroutine  bcs_measure_sisj(Amat,sisj_local)
use param
use lattice_param
use phiT_param
use model_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite)
complex(kind=8),intent(OUT)::sisj_local(Nbravais,Dtot)
integer,external::latt_label
integer,external::bound
complex(kind=8)::sisj(Nbravais,Nbravais)
integer::i,j,k,m,n
integer::cc(1:Dimen),ctmp

sisj_local=zero

k=1

sisj=zero

do i=1,Nbravais,1
  do j=1,Nbravais,1

        if(i.eq.j)then
             sisj(i,j)=sisj(i,j)+( Amat(i,i) + Amat(i+Nsite,i+Nsite)           &
     &                     -2.d0*(Amat(i,i)*Amat(i+Nsite,i+Nsite))             &
                           +2.d0*(Amat(i,i+Nsite)*Amat(i+Nsite,i)) )*3.d0
        else
             !up up up up
             sisj(i,j)=sisj(i,j)+Amat(i,i)*Amat(j,j)-Amat(i,j)*Amat(j,i)
             !dn dn dn dn
             sisj(i,j)=sisj(i,j)+Amat(i+Nsite,i+Nsite)*Amat(j+Nsite,j+Nsite)   &
     &                -Amat(i+Nsite,j+Nsite)*Amat(j+Nsite,i+Nsite)
             !other term R,S and S,R
             sisj(i,j)=sisj(i,j)+(-2.d0)*(Amat(i,j)*Amat(j+Nsite,i+Nsite)      &
     &                                   +Amat(i+Nsite,j+Nsite)*Amat(j,i))
             !other term R,R and S,S
             sisj(i,j)=sisj(i,j)-Amat(i,i)*Amat(j+Nsite,j+Nsite)               &
     &                          -Amat(i+Nsite,i+Nsite)*Amat(j,j)

!anomalous terms
             sisj(i,j)=sisj(i,j)-2.d0*(  Amat(i,j+Nsite) * (-Amat(j+Nsite,i))       &
     &                             +  (-Amat(j,i+Nsite))*   Amat(i+Nsite,j) )

              sisj(i,j)=sisj(i,j)+     (  Amat(i,j+Nsite) *   Amat(i+Nsite,j)         &
     &                              +  (-Amat(j,i+Nsite))* (-Amat(j+Nsite,i)))

        endif

        sisj(i,j)=sisj(i,j)/4.d0

  end do
end do


!average sisj to sisj_local
do i=1,Nbravais,1
  do j=1,Nbravais,1

     do m=1,Dimen,1
            !We calculate SiSj==>S1S(j-i+1)
            !so we need to focus on j-i+1
        ctmp=coor(j,m)-coor(i,m)+1
        cc(m)=bound(ctmp,Nl(m))
     end do
     n=latt_label(cc(1:Dimen))
     sisj_local(n,k)=sisj_local(n,k)+sisj(i,j)
  end do
end do
do n=1,Nbravais,1
   sisj_local(n,k)=sisj_local(n,k)/dcmplx(Nbravais)
end do

end subroutine  bcs_measure_sisj






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
complex(kind=8),intent(OUT)::sisj_local(Nbravais,Dtot)
integer,external::latt_label
integer,external::bound
complex(kind=8)::sisj(Nbravais,Nbravais)
integer::i,j,k,m,n
integer::cc(1:Dimen),ctmp

sisj_local=zero
do k=1,Dtot,1
   !Get sisj for k
   sisj=zero

   if(dtype.EQ.'c') then

     do i=1,Nbravais,1
        do j=1,Nbravais,1

 !          if(i.eq.j)then
 !            sisj(i,j)=sisj(i,j)+( Amat(i,i,k) + Amat(i+Nsite,i+Nsite,k)            &
 !    &                           -2.d0*(Amat(i,i,k)*Amat(i+Nsite,i+Nsite,k))        &
 !    &                           +2.d0*(Amat(i,i+Nsite,k)*Amat(i+Nsite,i,k)))*3.d0
 !          else
 !            !up up up up
 !            sisj(i,j)=sisj(i,j)+Amat(i,i,k)*Amat(j,j,k)-Amat(i,j,k)*Amat(j,i,k)
 !            !dn dn dn dn
 !            sisj(i,j)=sisj(i,j)+Amat(i+Nsite,i+Nsite,k)*Amat(j+Nsite,j+Nsite,k)   &
 !    &                -Amat(i+Nsite,j+Nsite,k)*Amat(j+Nsite,i+Nsite,k)      
 !            !other term R,S and S,R
 !            sisj(i,j)=sisj(i,j)+(-2.d0)*(Amat(i,j,k)*Amat(j+Nsite,i+Nsite,k)      &
 !    &                                   +Amat(i+Nsite,j+Nsite,k)*Amat(j,i,k))     &
 !    &                          +Amat(i,j+Nsite,k)*Amat(j+Nsite,i,k)               &
 !    &                          +Amat(i+Nsite,j,k)*Amat(j,i+Nsite,k)
 !            !other term R,R and S,S
 !            sisj(i,j)=sisj(i,j)-Amat(i,i,k)*Amat(j+Nsite,j+Nsite,k)               &
 !    &                          -Amat(i+Nsite,i+Nsite,k)*Amat(j,j,k)               &
 !    &                          +2.d0*Amat(i,i+Nsite,k)*Amat(j+Nsite,j,k)          &
 !    &                          +2.d0*Amat(i+Nsite,i,k)*Amat(j,j+Nsite,k)
!
!           endif
!         
!           sisj(i,j)=sisj(i,j)/4.d0


           !term 1
           if(i.eq.j) sisj(i,j)=sisj(i,j)+0.25d0*Amat(i,i,k)
           sisj(i,j)=sisj(i,j)+0.25d0*Amat(i,i,k)*Amat(j,j,k)-0.25d0*Amat(i,j,k)*Amat(j,i,k)
           !term 2
           if(i.eq.j) sisj(i,j)=sisj(i,j)+0.25d0*Amat(i+Nsite,i+Nsite,k)
           sisj(i,j)=sisj(i,j)+0.25d0*Amat(i+Nsite,i+Nsite,k)*Amat(j+Nsite,j+Nsite,k)- &
                            &  0.25d0*Amat(i+Nsite,j+Nsite,k)*Amat(j+Nsite,i+Nsite,k)
           !term 3
           sisj(i,j)=sisj(i,j)-0.25d0*Amat(i+Nsite,i+Nsite,k)*Amat(j,j,k)+ &
                            &  0.25d0*Amat(i+Nsite,j,k)*Amat(j,i+Nsite,k)
           !term 4
           sisj(i,j)=sisj(i,j)-0.25d0*Amat(i,i,k)*Amat(j+Nsite,j+Nsite,k)+ &
                            &  0.25d0*Amat(i,j+Nsite,k)*Amat(j+Nsite,i,k)
           !term 5
           if(i.eq.j) sisj(i,j)=sisj(i,j)+0.5d0*Amat(i,i,k)
           sisj(i,j)=sisj(i,j)+0.5d0*Amat(i,i+Nsite,k)*Amat(j+Nsite,j,k)- &
                            &  0.5d0*Amat(i,j,k)*Amat(j+Nsite,i+Nsite,k)
           !term 6
           if(i.eq.j) sisj(i,j)=sisj(i,j)+0.5d0*Amat(i+Nsite,i+Nsite,k)
           sisj(i,j)=sisj(i,j)+0.5d0*Amat(i+Nsite,i,k)*Amat(j,j+Nsite,k)- &
                            &  0.5d0*Amat(i+Nsite,j+Nsite,k)*Amat(j,i,k)
        end do
     end do

   else if(dtype.EQ.'d') then

     do i=1,Nbravais,1
        do j=1,Nbravais,1

           if(i.eq.j)then
             sisj(i,j)=sisj(i,j)+( Amat(i,i,k) + Amat(i+Nsite,i+Nsite,k)             &
     &                           -2.d0*(Amat(i,i,k)*Amat(i+Nsite,i+Nsite,k)))*3.d0
           else
             !up up up up
             sisj(i,j)=sisj(i,j)+Amat(i,i,k)*Amat(j,j,k)-Amat(i,j,k)*Amat(j,i,k)
             !dn dn dn dn
             sisj(i,j)=sisj(i,j)+Amat(i+Nsite,i+Nsite,k)*Amat(j+Nsite,j+Nsite,k)     &
     &                -Amat(i+Nsite,j+Nsite,k)*Amat(j+Nsite,i+Nsite,k)
             !other term R,S and S,R
             sisj(i,j)=sisj(i,j)+(-2.d0)*(Amat(i,j,k)*Amat(j+Nsite,i+Nsite,k)        &
     &                                   +Amat(i+Nsite,j+Nsite,k)*Amat(j,i,k))
             !other term R,R and S,S
             sisj(i,j)=sisj(i,j)-Amat(i,i,k)*Amat(j+Nsite,j+Nsite,k)                 &
     &                          -Amat(i+Nsite,i+Nsite,k)*Amat(j,j,k)
           endif

           sisj(i,j)=sisj(i,j)/4.d0

!           !term 1
!           if(i.eq.j) sisj(i,j)=sisj(i,j)+0.25d0*Amat(i,i,k)
!           sisj(i,j)=sisj(i,j)+0.25d0*Amat(i,i,k)*Amat(j,j,k)-0.25d0*Amat(i,j,k)*Amat(j,i,k)
!           !term 2
!           if(i.eq.j) sisj(i,j)=sisj(i,j)+0.25d0*Amat(i+Nsite,i+Nsite,k)
!           sisj(i,j)=sisj(i,j)+0.25d0*Amat(i+Nsite,i+Nsite,k)*Amat(j+Nsite,j+Nsite,k)- &
!                            & 0.25d0*Amat(i+Nsite,j+Nsite,k)*Amat(j+Nsite,i+Nsite,k)
!           !term 3
!           sisj(i,j)=sisj(i,j)-0.25d0*Amat(i+Nsite,i+Nsite,k)*Amat(j,j,k)
!
!           !term 4
!!           sisj(i,j)=sisj(i,j)-0.25d0*Amat(i,i,k)*Amat(j+Nsite,j+Nsite,k)
!
!           !term 5
!           if(i.eq.j) sisj(i,j)=sisj(i,j)+0.5d0*Amat(i,i,k)
!           sisj(i,j)=sisj(i,j)-0.5d0*Amat(i,j,k)*Amat(j+Nsite,i+Nsite,k)

!           !term 6
!           if(i.eq.j) sisj(i,j)=sisj(i,j)+0.5d0*Amat(i+Nsite,i+Nsite,k)
!           sisj(i,j)=sisj(i,j)-0.5d0*Amat(i+Nsite,j+Nsite,k)*Amat(j,i,k)

        end do
     end do

   end if

   !average sisj to sisj_local
   do i=1,Nbravais,1
      do j=1,Nbravais,1
         
         do m=1,Dimen,1
            !We calculate SiSj==>S1S(j-i+1)
            !so we need to focus on j-i+1
            ctmp=coor(j,m)-coor(i,m)+1
            cc(m)=bound(ctmp,Nl(m))
         end do
         n=latt_label(cc(1:Dimen))
         sisj_local(n,k)=sisj_local(n,k)+sisj(i,j)
      end do
   end do
   do n=1,Nbravais,1
      sisj_local(n,k)=sisj_local(n,k)/dcmplx(Nbravais)
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

   do k=1,Nbravais,1
      !The number of kl for k
      do n=1,Dimen,1
         kl(n)=coor(k,n)-1
      end do

      !fourier of the momentum k
      kk=zero
      do m=1,Nbravais,1
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
