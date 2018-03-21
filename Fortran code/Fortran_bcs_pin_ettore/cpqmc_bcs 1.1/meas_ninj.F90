subroutine bcs_measure_ninj(Amat,ninj_local)
use param
use lattice_param
use phiT_param
use model_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite)
complex(kind=8),intent(OUT)::ninj_local(Nbravais,Nbands,Nbands,Dtot)
integer,external::latt_label
integer,external::bound
complex(kind=8)::ninj(Nbravais,Nbands,Nbravais,Nbands)
integer::i,j,k,m,n,ib,jb,sitei,sitej,idn,jdn
integer::cc(1:Dimen),ctmp

ninj_local=zero
k=1

ninj=zero

do sitei=1,Nbravais,1
  do ib=1,Nbands
     do sitej=1,Nbravais,1
       do jb=1,Nbands

         i=sitei+(ib-1)*Nbravais
         idn=i+Nsite
         j=sitej+(jb-1)*Nbravais
         jdn=j+Nsite


         !term 1 <c+r_up cr_up c+oup coup>
         if(i.eq.j)ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+Amat(i,j)
         ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+Amat(i,i)*Amat(j,j)  &
     &                                                  -Amat(i,j)*Amat(j,i)

         !term 2 <c+r_up cr_up c+odn codn>
         ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+Amat(i,i)*Amat(jdn,jdn)

         !term 3 <c+r_dn cr_dn c+oup coup>
         ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+Amat(idn,idn)*Amat(j,j)

         !term 4 <c+r_dn cr_dn c+odn codn>
         if(i.eq.j)ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+Amat(idn,jdn)
         ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+Amat(idn,idn)*Amat(jdn,jdn) &
     &                                                  -Amat(idn,jdn)*Amat(jdn,idn)


!Anomalous terms
        ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)-Amat(i,j+Nsite)*Amat(i+Nsite,j)
        ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)-Amat(j,i+Nsite)*Amat(j+Nsite,i)

       enddo
     enddo
  enddo
enddo



!average ninj to ninj_local
 do ib=1,Nbands
   do jb=1,Nbands
     do i=1,Nbravais,1
       do j=1,Nbravais,1

         do m=1,Dimen,1
              !We calculate NiNj==>NN(j-i+1)
              !so we need to focus on j-i+1
             ctmp=coor(j,m)-coor(i,m)+1
             cc(m)=bound(ctmp,Nl(m))
         end do

         n=latt_label(cc(1:Dimen))
         ninj_local(n,ib,jb,k)=ninj_local(n,ib,jb,k)+ninj(i,ib,j,jb)
       end do
     end do
   end do
 end do
 do ib=1,Nbands
   do jb=1,Nbands
     do n=1,Nbravais,1
       ninj_local(n,ib,jb,k)=ninj_local(n,ib,jb,k)/dcmplx(Nbravais)
     end do
   end do
 end do




end subroutine bcs_measure_ninj


!--------------------------------
!measure the ninj_local from Amat
!--------------------------------
subroutine measure_ninj(Amat,ninj_local)
use param
use lattice_param
use phiT_param
use model_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::ninj_local(Nbravais,Nbands,Nbands,Dtot)
integer,external::latt_label
integer,external::bound
complex(kind=8)::ninj(Nbravais,Nbands,Nbravais,Nbands)
integer::i,j,k,m,n,ib,jb,sitei,sitej,idn,jdn
integer::cc(1:Dimen),ctmp

ninj_local=zero
do k=1,Dtot,1
   !Get ninj for k
   ninj=zero

   if(dtype.EQ.'c') then

     do sitei=1,Nbravais,1
       do ib=1,Nbands
         do sitej=1,Nbravais,1
           do jb=1,Nbands
         
             i=sitei+(ib-1)*Nbravais
             idn=i+Nsite
             j=sitej+(jb-1)*Nbravais
             jdn=j+Nsite
             
 
             !term 1 <c+r_up cr_up c+oup coup>
             if(i.eq.j)ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(i,j,k)
             ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(i,i,k)*Amat(j,j,k)  &
     &                                                      - Amat(i,j,k)*Amat(j,i,k)

             !term 2 <c+r_up cr_up c+odn codn>
             ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(i,i,k)*Amat(jdn,jdn,k)  &
     &                                                      - Amat(i,jdn,k)*Amat(jdn,i,k)

             !term 3 <c+r_dn cr_dn c+oup coup>
             ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(idn,idn,k)*Amat(j,j,k)  &
     &                                                      - Amat(idn,j,k)*Amat(j,idn,k)

             !term 4 <c+r_dn cr_dn c+odn codn>
             if(i.eq.j)ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(idn,jdn,k)
             ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(idn,idn,k)*Amat(jdn,jdn,k)  &
     &                                                      - Amat(idn,jdn,k)*Amat(jdn,idn,k)
            
           end do
         end do    
       end do
     end do

   else if(dtype.EQ.'d') then
  
     do sitei=1,Nbravais,1
       do ib=1,Nbands
         do sitej=1,Nbravais,1
           do jb=1,Nbands

             i=sitei+(ib-1)*Nbravais
             idn=i+Nsite
             j=sitej+(jb-1)*Nbravais
             jdn=j+Nsite

             !term 1 <c+r_up cr_up c+oup coup>
             if(i.eq.j)ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(i,j,k)
             ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(i,i,k)*Amat(j,j,k)  &
     &                                                      - Amat(i,j,k)*Amat(j,i,k)

             !term 2 <c+r_up cr_up c+odn codn>
             ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(i,i,k)*Amat(jdn,jdn,k)   

             !term 3 <c+r_dn cr_dn c+oup coup>
             ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(idn,idn,k)*Amat(j,j,k)

             !term 4 <c+r_dn cr_dn c+odn codn>
             if(i.eq.j)ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(idn,jdn,k)
             ninj(sitei,ib,sitej,jb)=ninj(sitei,ib,sitej,jb)+ Amat(idn,idn,k)*Amat(jdn,jdn,k) &
     &                                                      - Amat(idn,jdn,k)*Amat(jdn,idn,k)
           end do
         end do
       end do
     end do

   end if

   !average ninj to ninj_local
   do ib=1,Nbands
     do jb=1,Nbands
       do i=1,Nbravais,1
         do j=1,Nbravais,1
         
           do m=1,Dimen,1
              !We calculate NiNj==>NN(j-i+1)
              !so we need to focus on j-i+1
             ctmp=coor(j,m)-coor(i,m)+1
             cc(m)=bound(ctmp,Nl(m))
           end do

           n=latt_label(cc(1:Dimen))
           ninj_local(n,ib,jb,k)=ninj_local(n,ib,jb,k)+ninj(i,ib,j,jb)
         end do
       end do      
     end do
   end do
   do ib=1,Nbands
     do jb=1,Nbands
       do n=1,Nbravais,1
         ninj_local(n,ib,jb,k)=ninj_local(n,ib,jb,k)/dcmplx(Nbravais)
       end do
     end do 
   end do

end do
end subroutine measure_ninj




!----------------------------------------------------------------
!fourier sisj_l(1:Nsample,1:Nsite,i) to sk_l(1:Nsample,1:Nsite,i)
!----------------------------------------------------------------
subroutine fourier_nij(i)
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
end subroutine fourier_nij


subroutine bcs_measure_szsz(Amat,szsz_local)
use param
use lattice_param
use phiT_param
use model_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite)
complex(kind=8),intent(OUT)::szsz_local(Nbravais,Nbands,Nbands,Dtot)
integer,external::latt_label
integer,external::bound
complex(kind=8)::szsz(Nbravais,Nbands,Nbravais,Nbands)
integer::i,j,k,m,n,ib,jb,sitei,sitej,idn,jdn
integer::cc(1:Dimen),ctmp

szsz_local=zero
k=1

szsz=zero

do sitei=1,Nbravais,1
  do ib=1,Nbands
     do sitej=1,Nbravais,1
       do jb=1,Nbands

         i=sitei+(ib-1)*Nbravais
         idn=i+Nsite
         j=sitej+(jb-1)*Nbravais
         jdn=j+Nsite


         !term 1 <c+r_up cr_up c+oup coup>
         if(i.eq.j)szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(i,j)
         szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(i,i)*Amat(j,j)  &
     &                                                  -Amat(i,j)*Amat(j,i)

         !term 2 <c+r_up cr_up c+odn codn>
         szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)-Amat(i,i)*Amat(jdn,jdn)

         !term 3 <c+r_dn cr_dn c+oup coup>
         szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)-Amat(idn,idn)*Amat(j,j)

         !term 4 <c+r_dn cr_dn c+odn codn>
         if(i.eq.j)szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(idn,jdn)
         szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(idn,idn)*Amat(jdn,jdn)  &
     &                                                  -Amat(idn,jdn)*Amat(jdn,idn)


!Anomalous terms
        szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(i,j+Nsite)*Amat(i+Nsite,j)
        szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(j,i+Nsite)*Amat(j+Nsite,i)

       enddo
     enddo
  enddo
enddo



!average szsz to szsz_local
 do ib=1,Nbands
   do jb=1,Nbands
     do i=1,Nbravais,1
       do j=1,Nbravais,1

         do m=1,Dimen,1
              !We calculate NiNj==>NN(j-i+1)
              !so we need to focus on j-i+1
             ctmp=coor(j,m)-coor(i,m)+1
             cc(m)=bound(ctmp,Nl(m))
         end do

         n=latt_label(cc(1:Dimen))
         szsz_local(n,ib,jb,k)=szsz_local(n,ib,jb,k)+szsz(i,ib,j,jb)
       end do
     end do
   end do
 end do
 do ib=1,Nbands
   do jb=1,Nbands
     do n=1,Nbravais,1
       szsz_local(n,ib,jb,k)=(szsz_local(n,ib,jb,k)/dcmplx(Nbravais))*0.25d0
     end do
   end do
 end do




end subroutine bcs_measure_szsz






!--------------------------------
!measure the szsz_local from Amat
!--------------------------------
subroutine measure_szsz(Amat,szsz_local)
use param
use lattice_param
use phiT_param
use model_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::szsz_local(Nbravais,Nbands,Nbands,Dtot)
integer,external::latt_label
integer,external::bound
complex(kind=8)::szsz(Nbravais,Nbands,Nbravais,Nbands)
integer::i,j,k,m,n,ib,jb,sitei,sitej,idn,jdn
integer::cc(1:Dimen),ctmp


szsz_local=zero
do k=1,Dtot,1
   !Get szsz for k
   szsz=zero
   
   if(dtype.EQ.'c') then
   
     do sitei=1,Nbravais,1
       do ib=1,Nbands
         do sitej=1,Nbravais,1
           do jb=1,Nbands

             i=sitei+(ib-1)*Nbravais
             idn=i+Nsite
             j=sitej+(jb-1)*Nbravais
             jdn=j+Nsite

             !term 1 <c+r_up cr_up c+oup coup>
             if(i.eq.j)szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(i,j,k)
             szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(i,i,k)*Amat(j,j,k)  &
     &                                                      -Amat(i,j,k)*Amat(j,i,k)

            !term 2 <c+r_up cr_up c+odn codn>
            szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)-Amat(i,i,k)*Amat(jdn,jdn,k) &
     &                                                     +Amat(idn,j,k)*Amat(j,idn,k)

            !term 3 <c+r_dn cr_dn c+oup coup>
            szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)-Amat(idn,idn,k)*Amat(j,j,k)  &
     &                                                     +Amat(idn,j,k)*Amat(j,idn,k)

            !term 4 <c+r_dn cr_dn c+odn codn>
            if(i.eq.j)szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(idn,jdn,k)
            szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(idn,idn,k)*Amat(jdn,jdn,k) &
     &                                                     -Amat(idn,jdn,k)*Amat(jdn,idn,k)

          end do
         end do
       end do
     end do

   else if(dtype.EQ.'d') then

     do sitei=1,Nbravais,1
       do ib=1,Nbands
         do sitej=1,Nbravais,1
           do jb=1,Nbands

             i=sitei+(ib-1)*Nbravais
             idn=i+Nsite
             j=sitej+(jb-1)*Nbravais
             jdn=j+Nsite

             !term 1 <c+r_up cr_up c+oup coup>
             if(i.eq.j)szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(i,j,k)
             szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(i,i,k)*Amat(j,j,k)  &
     &                                                      -Amat(i,j,k)*Amat(j,i,k)

            !term 2 <c+r_up cr_up c+odn codn>
            szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)-Amat(i,i,k)*Amat(jdn,jdn,k)

            !term 3 <c+r_dn cr_dn c+oup coup>
            szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)-Amat(idn,idn,k)*Amat(j,j,k)

            !term 4 <c+r_dn cr_dn c+odn codn>
            if(i.eq.j)szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(idn,jdn,k)
            szsz(sitei,ib,sitej,jb)=szsz(sitei,ib,sitej,jb)+Amat(idn,idn,k)*Amat(jdn,jdn,k)  &
     &                                                     -Amat(idn,jdn,k)*Amat(jdn,idn,k)

          end do
         end do
       end do
     end do     
 
   endif

   !average szsz to szsz_local
   do ib=1,Nbands
     do jb=1,Nbands
       do i=1,Nbravais,1
         do j=1,Nbravais,1

           do m=1,Dimen,1
              !We calculate NiNj==>NN(j-i+1)
              !so we need to focus on j-i+1
             ctmp=coor(j,m)-coor(i,m)+1
             cc(m)=bound(ctmp,Nl(m))
           end do

           n=latt_label(cc(1:Dimen))
           szsz_local(n,ib,jb,k)=szsz_local(n,ib,jb,k)+szsz(i,ib,j,jb)
         end do
       end do
     end do
   end do
   do ib=1,Nbands
     do jb=1,Nbands
       do n=1,Nbravais,1
         szsz_local(n,ib,jb,k)=(szsz_local(n,ib,jb,k)/dcmplx(Nbravais))*0.25d0
       end do
     end do
   end do
    
end do
end subroutine measure_szsz

    
