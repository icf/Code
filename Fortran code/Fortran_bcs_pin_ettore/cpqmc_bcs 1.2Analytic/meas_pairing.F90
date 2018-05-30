subroutine bcs_measure_pairing(Amat,didj_local)
use param
use lattice_param
use phiT_param
use model_param
use method_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite)
complex(kind=8),intent(OUT)::didj_local(Nbravais,Npair,Dtot)
integer,external::latt_label
integer,external::bound
real(kind=8)::sign2d
complex(kind=8)::tmpcorr
integer::i,j,k,m,n,ii,ib,jj,jb,delta,den,dir,dirp,rho,kk,iprime,jprime,ipair
integer::cc(1:Dimen),ctmp

if(Npair.eq.0)return
didj_local=zero
k=1

if(ipinn.eq.0)then
 do i=1,Nsite,1
  call unit_cell(i,ii,ib)
  do j=1,Nsite,1
    call unit_cell(j,jj,jb)

    do m=1,Dimen,1
      ctmp=coor(jj,m)-coor(ii,m)+1
      cc(m)=bound(ctmp,Nl(m))
    end do
    n=latt_label(cc(1:Dimen))

    tmpcorr=Amat(j+Nsite,i+Nsite)*Amat(j,i)-Amat(j,j+Nsite)*Amat(i+Nsite,i)

    didj_local(n,1,k)=didj_local(n,1,k)+tmpcorr               !onsite s wave pairing

  enddo
 enddo
else

  i=Nbravais/2
  ii=Nbravais/2                     !reference point in the middle of the lattice 
  ib=1

  do j=1,Nsite,1
    call unit_cell(j,jj,jb)

    do m=1,Dimen,1
      ctmp=coor(jj,m)-coor(ii,m)+1
      cc(m)=bound(ctmp,Nl(m))
    end do
    n=latt_label(cc(1:Dimen))

    tmpcorr=Amat(j+Nsite,i+Nsite)*Amat(j,i)-Amat(j,j+Nsite)*Amat(i+Nsite,i)

    didj_local(n,1,k)=didj_local(n,1,k)+tmpcorr               !onsite s wave pairing

  enddo
endif



if(Npair.gt.1)then
 if(ipinn.eq.0)then
  do i=1,Nsite,1
  call unit_cell(i,ii,ib)
  do j=1,Nsite,1
    call unit_cell(j,jj,jb)

    do m=1,Dimen,1
      ctmp=coor(jj,m)-coor(ii,m)+1
      cc(m)=bound(ctmp,Nl(m))
    end do
    n=latt_label(cc(1:Dimen))

    do delta=1,Dimen
      do dir=0,1
        cc(:)=coor(jj,:)
        ctmp=cc(delta)+(-1)**dir
        cc(delta)=bound(ctmp,Nl(delta))
        jprime=latt_label(cc(1:Dimen))+(jb-1)*Nbravais
        do rho=1,Dimen
          if(delta.eq.rho)then
            sign2d=1.d0
          else
            sign2d=-1.d0
          endif
          do dirp=0,1
            cc(:)=coor(ii,:)
            ctmp=cc(rho)+(-1)**dirp
            cc(rho)=bound(ctmp,Nl(rho))
            iprime=latt_label(cc(1:Dimen))+(ib-1)*Nbravais

            tmpcorr=Amat(jprime+Nsite,iprime+Nsite)*Amat(j,i)+             &
                    (-Amat(j,jprime+Nsite))*Amat(i+Nsite,iprime)+          &
                    Amat(jprime+Nsite,i+Nsite)*Amat(j,iprime)+             &
                    ( Amat(j,jprime+Nsite))*(-Amat(iprime+Nsite,i))+   &
                    Amat(jprime,i)*Amat(j+Nsite,iprime+Nsite)+             &
                    (-Amat(jprime,j+Nsite))*Amat(i+Nsite,iprime)+          &
                    Amat(jprime,iprime)*Amat(j+Nsite,i+Nsite)+             &
                    ( Amat(jprime,j+Nsite))*(-Amat(iprime+Nsite,i))   
                    

            didj_local(n,2,k)=didj_local(n,2,k)+tmpcorr                      !extended nn s wave pairing
            if(Npair.gt.2)didj_local(n,3,k)=didj_local(n,3,k)+tmpcorr*sign2d !nn dwave pairing
          enddo
        enddo
      enddo
    enddo
  enddo
  enddo
 else

  i=Nbravais/2
  ii=Nbravais/2
  ib=1

  do j=1,Nsite,1
    call unit_cell(j,jj,jb)

    do m=1,Dimen,1
      ctmp=coor(jj,m)-coor(ii,m)+1
      cc(m)=bound(ctmp,Nl(m))
    end do
    n=latt_label(cc(1:Dimen))

    do delta=1,Dimen
      do dir=0,1
        cc(:)=coor(jj,:)
        ctmp=cc(delta)+(-1)**dir
        cc(delta)=bound(ctmp,Nl(delta))
        jprime=latt_label(cc(1:Dimen))+(jb-1)*Nbravais
        do rho=1,Dimen
          if(delta.eq.rho)then
            sign2d=1.d0
          else
            sign2d=-1.d0
          endif
          do dirp=0,1
            cc(:)=coor(ii,:)
            ctmp=cc(rho)+(-1)**dirp
            cc(rho)=bound(ctmp,Nl(rho))
            iprime=latt_label(cc(1:Dimen))+(ib-1)*Nbravais

            tmpcorr=Amat(jprime+Nsite,iprime+Nsite)*Amat(j,i)+             &
                    (-Amat(j,jprime+Nsite))*Amat(i+Nsite,iprime)+          &
                    Amat(jprime+Nsite,i+Nsite)*Amat(j,iprime)+             &
                    ( Amat(j,jprime+Nsite))*(-Amat(iprime+Nsite,i))+   &
                    Amat(jprime,i)*Amat(j+Nsite,iprime+Nsite)+             &
                    (-Amat(jprime,j+Nsite))*Amat(i+Nsite,iprime)+          &
                    Amat(jprime,iprime)*Amat(j+Nsite,i+Nsite)+             &
                    ( Amat(jprime,j+Nsite))*(-Amat(iprime+Nsite,i))


            didj_local(n,2,k)=didj_local(n,2,k)+tmpcorr                      !extended s wave pairing
            if(Npair.gt.2)didj_local(n,3,k)=didj_local(n,3,k)+tmpcorr*sign2d !dwave pairing
          enddo
        enddo
      enddo
    enddo
  enddo


 endif


endif


if(ipinn.eq.0)then
 do ipair=1,Npair
  do n=1,Nbravais,1
    didj_local(n,ipair,k)=didj_local(n,ipair,k)/dcmplx(Nbravais)
  enddo
 enddo
endif
    

end subroutine bcs_measure_pairing



!--------------------------------
!measure the pairing correlation function from Amat
!--------------------------------
subroutine measure_pairing(Amat,didj_local)
use param
use lattice_param
use phiT_param
use model_param
use method_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::didj_local(Nbravais,Npair,Dtot)
integer,external::latt_label
integer,external::bound
real(kind=8)::sign2d
complex(kind=8)::tmpcorr
integer::i,j,k,m,n,ii,ib,jj,jb,delta,den,dir,dirp,rho,kk,iprime,jprime,ipair
integer::cc(1:Dimen),ctmp

if(Npair.eq.0)return
didj_local=zero
do k=1,Dtot,1

  do i=1,Nsite,1
    call unit_cell(i,ii,ib)
    do j=1,Nsite,1
      call unit_cell(j,jj,jb)

      do m=1,Dimen,1
        ctmp=coor(jj,m)-coor(ii,m)+1
        cc(m)=bound(ctmp,Nl(m))
      end do
      n=latt_label(cc(1:Dimen))

      tmpcorr=Amat(j+Nsite,i+Nsite,k)*Amat(j,i,k)

      didj_local(n,1,k)=didj_local(n,1,k)+tmpcorr               !onsite s wave pairing

    enddo
  enddo


  if(Npair.gt.1)then

   !get didj_local
   do i=1,Nsite,1

      call unit_cell(i,ii,ib)

      do j=1,Nsite,1

         call unit_cell(j,jj,jb)

!         write(*,*)'Element i ',i
!         write(*,*)'Site ',ii
!         write(*,*)'Band ',ib
 
!         write(*,*)
!         write(*,*)

!         write(*,*)'Element j ',j
!         write(*,*)'Site ',jj
!         write(*,*)'Band ',jb

         !stop'DEB'

     
         do m=1,Dimen,1
            !We calculate SiSj==>S1S(j-i+1)
            !so we need to focus on j-i+1
            ctmp=coor(jj,m)-coor(ii,m)+1
            cc(m)=bound(ctmp,Nl(m))
         end do
         n=latt_label(cc(1:Dimen))

!         write(*,*)
!         write(*,*)'Punto R ',n
!         write(*,*)

         do delta=1,Dimen

           do dir=0,1

             cc(:)=coor(jj,:)
             ctmp=cc(delta)+(-1)**dir
             cc(delta)=bound(ctmp,Nl(delta))
             jprime=latt_label(cc(1:Dimen))+(jb-1)*Nbravais


!             write(*,*)
         !    write(*,*)'delta ',delta    
          !   write(*,*)
           !  write(*,*)
         !    write(*,*)'dir ',dir
         !    write(*,*) 
         !    write(*,*)
!             write(*,*)'j,jprime ',j,jprime
!             write(*,*)
             !stop'DEB'

             do rho=1,Dimen

               if(delta.eq.rho)then
                 sign2d=1.d0
               else
                 sign2d=-1.d0
               endif

               do dirp=0,1

                 cc(:)=coor(ii,:)
                 ctmp=cc(rho)+(-1)**dirp
                 cc(rho)=bound(ctmp,Nl(rho))
                 iprime=latt_label(cc(1:Dimen))+(ib-1)*Nbravais

!                 write(*,*)
         !        write(*,*)'rho ',rho
         !        write(*,*)
         !        write(*,*)
         !        write(*,*)'dirp ',dirp
         !        write(*,*)
         !        write(*,*)
!                 write(*,*)'i,iprime ',i,iprime
!                 write(*,*)

                 if(dtype.EQ.'c') then
                   tmpcorr=(Amat(jprime+Nsite,iprime+Nsite,k)*Amat(j,i,k)-   &
                            Amat(jprime+Nsite,i,k)*Amat(j,iprime+Nsite,k))-  &
                           (Amat(jprime+Nsite,iprime,k)*Amat(j,i+Nsite,k)-   &
                            Amat(jprime+Nsite,i+Nsite,k)*Amat(j,iprime,k))-  &
                           (Amat(jprime,iprime+Nsite,k)*Amat(j+Nsite,i,k)-   &
                            Amat(jprime,i,k)*Amat(j+Nsite,iprime+Nsite,k))+  &
                           (Amat(jprime,iprime,k)*Amat(j+Nsite,i+Nsite,k)-   &
                            Amat(jprime,i+Nsite,k)*Amat(j+Nsite,iprime,k))   

                 elseif(dtype.EQ.'d') then

                    tmpcorr=Amat(jprime+Nsite,iprime+Nsite,k)*Amat(j,i,k)+   &
                            Amat(jprime+Nsite,i+Nsite,k)*Amat(j,iprime,k)+   &
                            Amat(jprime,i,k)*Amat(j+Nsite,iprime+Nsite,k)+   &
                            Amat(jprime,iprime,k)*Amat(j+Nsite,i+Nsite,k)   

!Shiwei    tmpcorr=tmpfac*
!     +                  ( onemG(i1,i2,1)*onemG(i1pr1,i2pr2,2) per me 1 riga
!     +                   +onemG(i1pr1,i2pr2,1)*onemG(i1,i2,2) per me 4 riga
!     +                   +onemG(i1pr1,i2,1)*onemG(i1,i2pr2,2) per me 3 riga
!     +                   +onemG(i1,i2pr2,1)*onemG(i1pr1,i2,2) per me 2 riga

                 endif

                 didj_local(n,2,k)=didj_local(n,2,k)+tmpcorr           !s wave pairing
                 if(Npair.gt.2)didj_local(n,3,k)=didj_local(n,3,k)+tmpcorr*sign2d    !d wave pairing
                 !write(*,*)'n, estim ',n,didj_local(n,k)

               enddo               
             enddo  
           enddo
         enddo
      end do
    end do


   endif
   do ipair=1,Npair
    do n=1,Nbravais,1
      didj_local(n,ipair,k)=didj_local(n,ipair,k)/dcmplx(Nbravais)
!      write(*,*)'MEAS_PAIR ipair,n, estim ',ipair,n,didj_local(n,ipair,k)      
    end do
!    write(*,*)
!    write(*,*)
   enddo
!   stop'DEB'

end do
end subroutine measure_pairing


subroutine unit_cell(i,site,band)
use param
use lattice_param
use mc_loop_param
use meas_param
integer,intent(IN)::i
integer,intent(OUT)::site,band

if(i.le.Nbravais)then
  site=i
  band=1
elseif(i.le.2*Nbravais)then
  site=i-Nbravais
  band=2
elseif(i.le.3*Nbravais)then
  site=i-2*Nbravais
  band=3
else
  write(*,*)'Problem in unit_cell ',i
  stop
endif

end subroutine unit_cell

