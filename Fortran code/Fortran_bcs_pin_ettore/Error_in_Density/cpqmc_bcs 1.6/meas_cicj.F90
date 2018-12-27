subroutine measure_obdm(Amat,obdm_local)
use param
use lattice_param
use phiT_param
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::obdm_local(2*Nsite,2*Nsite,Dtot)

!write(101,*)2*Nsite*2*Nsite*Dtot
!write(101,*)'Amat(1,1,1) = ',Amat(1,1,1)
!write(101,*)
!flush(101)


call zcopy(2*Nsite*2*Nsite*Dtot,Amat(1,1,1),1,obdm_local(1,1,1),1)

!write(101,*)'obdm_local(1,1,1) ',obdm_local(1,1,1)
!write(101,*)obdm_local-Amat
!write(101,*)
!flush(101)
!stop 'deb'


end subroutine measure_obdm



!--------------------------------
!measure the cacb_local from Amat
!--------------------------------
subroutine measure_cacb(Amat,cacb_local)
use param
use lattice_param
use phiT_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::cacb_local(Nbravais,Nbands,Nbands,Dtot)
integer,external::latt_label
integer,external::bound
integer::i,j,k,m,n,sitei,sitej,ib1,ib2,ii,jj,ib,jb,iprime
integer::cc(1:Dimen),ctmp

cacb_local=zero
do k=1,Dtot,1

 if(Nbands.eq.1)then

   !average cicj to cicj_local
   do i=1,Nbravais,1
      do j=1,Nbravais,1

         do m=1,Dimen,1
            !We calculate (Ci^+)Cj==>C(i-j+1)^+ C1
            !so we need to focus on i-j+1
            ctmp=coor(i,m)-coor(j,m)+1
            cc(m)=bound(ctmp,Nl(m))
         end do
         n=latt_label(cc(1:Dimen))
         cacb_local(n,1,1,k)=cacb_local(n,1,1,k)+0.5*(Amat(i,j,k)+Amat(i+Nsite,j+Nsite,k))

      end do
   end do
   do n=1,Nbravais,1
      cacb_local(n,1,1,k)=cacb_local(n,1,1,k)/dcmplx(Nsite)
   end do

 elseif(Nbands.eq.3)then


!   write(*,*)
!   write(*,*)'**********************************'
!   write(*,*)'Amat  '
!   write(*,*)
!   do sitei=1,Nsite
!     do sitej=1,Nsite
!       write(*,*)sitei,sitej,Amat(sitei,sitej,1)
!     enddo
!   enddo
!   write(*,*)
!   write(*,*)

!   stop


   do i=1,Nsite,1

      call unit_cell(i,ii,ib)

      do j=1,Nsite,1

        call unit_cell(j,jj,jb)

        do m=1,Dimen,1
          ctmp=coor(jj,m)-coor(ii,m)+1
          cc(m)=bound(ctmp,Nl(m))
        enddo
        n=latt_label(cc(1:Dimen))

        cacb_local(n,ib,jb,k)=cacb_local(n,ib,jb,k)+0.5*(Amat(i,j,k)+Amat(i+Nsite,j+Nsite,k))        
      enddo
   enddo
   do ib=1,Nbands
     do jb=1,Nbands
       do n=1,Nbravais
         cacb_local(n,jb,ib,k)=cacb_local(n,jb,ib,k)/dcmplx(Nbravais)
       enddo
     enddo
   end do

  ! write(*,*)
  ! write(*,*)'**********************************'
  ! write(*,*)'Amat mediata  '
  ! write(*,*)
  ! do sitei=1,Nsite
  !     write(*,*)1,sitei,cicj_local(sitei,k)
  ! enddo
  ! write(*,*)
  ! write(*,*)



 endif

end do
end subroutine measure_cacb



!--------------------------------
!measure the density <c+(R,alpha)c(R,beta)> from Amat
!--------------------------------
subroutine measure_nofr(Amat,nofr_local)
use param
use lattice_param
use phiT_param
use model_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::nofr_local(Nsite,2,2,Dtot)
integer::i,k

do k=1,Dtot,1

  if(dtype.EQ.'c') then

    do i=1,Nsite,1
      nofr_local(i,1,1,k)=Amat(i,i,k)
      nofr_local(i,1,2,k)=Amat(i,i+Nsite,k)
      nofr_local(i,2,1,k)=Amat(i+Nsite,i,k)
      nofr_local(i,2,2,k)=Amat(i+Nsite,i+Nsite,k)
      !write(*,*)i,Amat(i,i+Nsite,k),Amat(i+Nsite,i,k)
    enddo
    !stop

  else if(dtype.EQ.'d') then

    do i=1,Nsite,1
      nofr_local(i,1,1,k)=Amat(i,i,k)
      nofr_local(i,1,2,k)=zero
      nofr_local(i,2,1,k)=zero
      nofr_local(i,2,2,k)=Amat(i+Nsite,i+Nsite,k)
    enddo

  endif

end do
end subroutine measure_nofr



subroutine bcs_measure_nofr(Amat,nofr_local)
use param
use lattice_param
use phiT_param
use model_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite)
complex(kind=8),intent(OUT)::nofr_local(Nsite,2,2,Dtot)
integer::i,k

k=1

do i=1,Nsite,1
  nofr_local(i,1,1,k)=Amat(i,i)
  nofr_local(i,1,2,k)=zero
  nofr_local(i,2,1,k)=zero
  nofr_local(i,2,2,k)=Amat(i+Nsite,i+Nsite)
enddo

end subroutine bcs_measure_nofr




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
integer::i,j,k,m,n,sitei,sitej,ib1,ib2,ii,jj,ib,jb,iprime
integer::cc(1:Dimen),ctmp

cicj_local=zero
do k=1,Dtot,1

 if(Nbands.eq.1)then

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

 elseif(Nbands.eq.3)then


  ! write(*,*)
  ! write(*,*)'**********************************'
  ! write(*,*)'Amat  '
  ! write(*,*)
  ! do sitei=1,Nsite
  !   do sitej=1,Nsite
  !     write(*,*)sitei,sitej,Amat(sitei,sitej,1)
  !   enddo
  ! enddo
  ! write(*,*)
  ! write(*,*)


   do i=1,Nsite,1

      call unit_cell(i,ii,ib)
  
      do j=1,Nsite,1

        call unit_cell(j,jj,jb)
 
        do m=1,Dimen,1
          ctmp=coor(jj,m)-coor(ii,m)+1
          cc(m)=bound(ctmp,Nl(m))
        enddo
        n=latt_label(cc(1:Dimen))

        if(ib.eq.jb)then

          iprime=(ib-1)*Nbravais+n

          cicj_local(iprime,k)=cicj_local(iprime,k)+Amat(i,j,k)
          cicj_local(iprime+Nsite,k)=cicj_local(iprime+Nsite,k)+Amat(i,j,k)
        
        endif   
 
     enddo
   enddo
   do n=1,2*Nsite,1
      cicj_local(n,k)=cicj_local(n,k)/dcmplx(Nbravais)
   end do

  ! write(*,*)
  ! write(*,*)'**********************************'
  ! write(*,*)'Amat mediata  '
  ! write(*,*)
  ! do sitei=1,Nsite
  !     write(*,*)1,sitei,cicj_local(sitei,k)
  ! enddo
  ! write(*,*)
  ! write(*,*)



 endif

end do
end subroutine measure_cicj




!----------------------------------------------------------------
!fourier cicj_l(1:Nsample,1:Nsite,i) to ck_l(1:Nsample,1:Nsite,i)
!----------------------------------------------------------------
!subroutine fourier_cij(i)
!use param
!use lattice_param
!use mc_loop_param
!use meas_param
!implicit none
!integer,intent(IN)::i
!integer::j,k,m,n,ii,ib
!integer::kl(Dimen)
!complex(kind=8)::kkup,kkdn,kkupd,kkdnd,kkupx,kkdnx,kkupy,kkdny
!real(kind=8)::ph


!if(Nbands.eq.1)then

!  do j=1,Nsamples,1
!
!     do k=1,Nsite,1
      !The number of kl for k
!        do n=1,Dimen,1
!           kl(n)=coor(k,n)-1
!        end do

      !fourier of the momentum k
!        kkup=zero;kkdn=zero
!        do m=1,Nsite,1
!           ph=0.d0
!           do n=1,Dimen,1
!              ph=ph+2.d0*Pi*dble(kl(n)*(coor(m,n)-1))/dble(Nl(n))
!           end do
!           kkup=kkup+exp(Xi*ph)*cicj_l(j,m,i)
!           kkdn=kkdn+exp(Xi*ph)*cicj_l(j,m+Nsite,i)
!          !write(*,*) exp(Xi*ph);pause
!        end do
!        ck_l(j,k,i)=dble(kkup)
!        ck_l(j,k+Nsite,i)=dble(kkdn)
!     end do

!  end do

!elseif(Nbands.eq.3)then
! 
!  do j=1,Nsamples,1
!!
!     do k=1,Nbravais,1
!      !The number of kl for k

!        do n=1,Dimen,1
!           kl(n)=coor(k,n)-1
!        end do

      !fourier of the momentum k
!        kkup=zero;kkdn=zero
!        do m=1,Nsite,1
!          call unit_cell(m,ii,ib)

!          ph=0.d0
!          do n=1,Dimen,1
!            ph=ph+2.d0*Pi*dble(kl(n)*(coor(ii,n)-1))/dble(Nl(n))
!          end do
!          kkup=kkup+exp(Xi*ph)*cicj_l(j,m,i)
!          kkdn=kkdn+exp(Xi*ph)*cicj_l(j,m+Nsite,i)

!        enddo
!        ck_l(j,k,i)=dble(kkup)
!        ck_l(j,k+Nsite,i)=dble(kkdn)
      
!     end do

!  end do

  ! write(*,*)
  ! write(*,*)'**********************************'
  ! write(*,*)'Trasformata di Fourier  '
  ! write(*,*)
  ! do k=1,Nsite
  !     write(*,*)k,ck_l(1,k,i)
  ! enddo
  ! write(*,*)
  ! write(*,*)



!endif
!end subroutine fourier_cij
