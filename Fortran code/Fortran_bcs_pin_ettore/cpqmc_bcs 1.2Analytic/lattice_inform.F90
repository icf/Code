!This subroutine get the Hamilton of the no interaction lattice
subroutine set_lattice()
use lattice_param
use model_param
use mpi_serial_param
implicit none
integer::mi,mj,i1,i2
real(kind=8)::eps,Vext,x,y

 eps=1.d-5

 !Set the number of lattice and Nhop
 if(set.NE.1) then
   Nbravais=1
   do mi=1,Dimen,1
    Nbravais=Nbravais*Nl(mi)
   end do
   Nsite=Nbands*Nbravais
   if(Nbands.eq.1)then
    Nhop=Nsite*2*Dimen
   elseif(Nbands.eq.3)then
    Nhop=16*Nbravais
    if(Dimen.ne.2)then
     write(*,*)'Three bands only in 2d'
     stop
    endif
   endif
 end if
 if(Ntot.GT.2*Nsite) then
   write(*,*) "Number of electrons larger than Nsite"
   call mystop
 end if

 call allocate_lattice_array()
 call set_lattice_hop()

 Hzero=(0.d0,0.d0)
 do mi=1,Nhop,1
    Hzero(sit(mi,1),sit(mi,2))=Hzero(sit(mi,1),sit(mi,2))+hopt(mi)
    Hzero(sit(mi,1)+Nsite,sit(mi,2)+Nsite)=Hzero(sit(mi,1)+Nsite,sit(mi,2)+Nsite)+hopt(mi)
 end do

! write(*,*)'********',Hzero(3,5)
 if(Nbands.eq.3)then
  do mi=1,Nbravais
   Hzero(mi,mi)=Hzero(mi,mi)+epsd
   Hzero(mi+Nsite,mi+Nsite)=Hzero(mi+Nsite,mi+Nsite)+epsd
  enddo
  do mi=Nbravais+1,Nbands*Nbravais
   Hzero(mi,mi)=Hzero(mi,mi)+epsp
   Hzero(mi+Nsite,mi+Nsite)=Hzero(mi+Nsite,mi+Nsite)+epsp
  enddo
 endif

!  write(*,*)'********',Hzero(3,5)
 if(pinndir.eq.'z')then
   do mi=1,Nl(1)  !Pinning field on z axis
    Hzero(mi,mi)=Hzero(mi,mi)+Hpinn*((-1)**mi)
    Hzero(mi+Nsite,mi+Nsite)=Hzero(mi+Nsite,mi+Nsite)-Hpinn*((-1)**mi)
   enddo
!Both sides are pinned
   do mi=Nbravais-Nl(1)+1,Nbravais
     Hzero(mi,mi)=Hzero(mi,mi)+Hpinn*((-1)**(1+mi-(Nbravais-Nl(1))))
     Hzero(mi+Nsite,mi+Nsite)=Hzero(mi+Nsite,mi+Nsite)-Hpinn*((-1)**(1+mi-(Nbravais-Nl(1))))
   enddo
 elseif(pinndir.eq.'x')then
   do mi=1,Nl(1)
    Hzero(mi,mi+Nsite)=Hzero(mi,mi+Nsite)+Hpinn*((-1)**mi)
    Hzero(mi+Nsite,mi)=Hzero(mi+Nsite,mi)+Hpinn*((-1)**mi)
   enddo
 endif

 if(I_Vext.eq.1)then
     if(rank.eq.0)then
       write(*,*)
       write(*,*)'Read external field from file '
     endif
     open(2,file='ExternalScalarField',status='old')
       do mi=1,Nsite
         read(2,*)Vext
         Hzero(mi,mi)=Hzero(mi,mi)+Vext
         Hzero(mi+Nsite,mi+Nsite)=Hzero(mi+Nsite,mi+Nsite)+Vext
       enddo
     close(2)     
 endif

 if(I_ReadH0.eq.1)then
   open(3,file='hzero_rspace.out',status='old')
   do mi=1,2*Nsite
     do mj=1,2*Nsite
       read(3,*)i1,i2,x,y
       Hzero(mi,mj)=dcmplx(x,y)
     enddo
   enddo
   close(3)
 endif

 if(rank.eq.0)then
   open(2,file='hzero.info',status='unknown')
   write(*,*)'**************************************'
   write(*,*)
   write(*,*)
   do mi=1,2*Nsite
    do mj=1,2*Nsite
     write(2,*)mi,mj,dble(Hzero(mi,mj)),aimag(Hzero(mi,mj))
    enddo
  enddo
  close(2)
  open(2,file='hzero_up.info',status='unknown')
   write(*,*)'**************************************'
   write(*,*)
   write(*,*)
   do mi=1,Nsite
    do mj=1,Nsite
     write(2,*)mi,mj,dble(Hzero(mi,mj)),aimag(Hzero(mi,mj))
    enddo
   enddo
  close(2)
  open(2,file='hop.dat',status='unknown')
   do mi=1,2*Nsite
    do mj=1,2*Nsite
     if(abs(Hzero(mi,mj)).gt.eps)then
      write(2,*)mi,mj,Hzero(mi,mj)
     endif
    enddo
   enddo
  close(2)
 endif
 call check_Hermite_c(Hzero,2*Nsite)

! stop'DEBUG'


 call set_Tmatrix()

! stop'DEBUG'

 !------------------------------------------------
 !If (set.EQ.1): We need to pay attention that:
 !coor(:,:) and Tmatrix(:,:) will not be set here.
 !------------------------------------------------
 
end subroutine set_lattice



!-----------------------------------------------------------
!This subroutine get the sit(Nhop,2) and hopt(Nhop) and coor
!-----------------------------------------------------------
subroutine set_lattice_hop()
use lattice_param
implicit none
integer::mi,mj

if(set.EQ.1) then
  open(unit=10,file='hop',status='old')
    do mi=1,Nhop,1
       read(10,*) sit(mi,1),sit(mi,2),hopt(mi)
    end do
  close(10)
else if(set.EQ.2) then
  call setnumber()
  if(Nbands.eq.1)then
   if(boundary.eq.'c')then
     call sethopTBC()
   elseif(boundary.eq.'o')then
     call sethopTBC_open()
   elseif(boundary .EQ. 'a')then
     call sethopTBC_allopen()
   endif
  else
   if(boundary.eq.'o')then
     call sethopTBC_multi_open()
   elseif(boundary.eq.'a')then
     call sethopTBC_multi_allopen()
   else
     call sethopTBC_multi()
   endif
  endif
else
  write(*,*) "something is wrong with the set"
  call mystop
end if

end subroutine set_lattice_hop


!-----------------------------------------
!we label the lattice by 1~Nsite,set
!by (x,y,z) the coor(Nsite,Dimen) is
!dependent on the dimension of the lattice
!It can be think as the coordinate of the
!lattice point is from 1~Nl(:)
!-----------------------------------------
subroutine setnumber
use lattice_param
use mpi_serial_param
implicit none
integer::i,j,k
integer::ntemp,den
real(kind=8)::x,y

if(Nbands.eq.1)then
 do i=1,Nsite,1
    ntemp=i-1
    do j=Dimen,1,-1
       den=1
       do k=1,j-1,1
          den=den*Nl(k)
       end do !den is z coor Nl(1)*Nl(2)
              !       y coor Nl(1)
              !       x coor 1
       coor(i,j)=ntemp/den
       ntemp=ntemp-coor(i,j)*den
       coor(i,j)=coor(i,j)+1
    end do
 end do
else
 do i=1,Nbravais,1
    ntemp=i-1
    do j=Dimen,1,-1
       den=1
       do k=1,j-1,1
          den=den*Nl(k)
       end do !den is z coor Nl(1)*Nl(2)
              !       y coor Nl(1)
              !       x coor 1
       coor(i,j)=ntemp/den
       ntemp=ntemp-coor(i,j)*den
       coor(i,j)=coor(i,j)+1
    end do
 end do
 do i=1,Nbravais,1
   do j=Dimen,1,-1
     coor(i+Nbravais,j)=coor(i,j)
     coor(i+2*Nbravais,j)=coor(i,j)
   enddo
 enddo
endif

if(rank.eq.0)then
  open(2,file='lattice.info',status='unknown')
  do i=1,Nbravais,1
    !x=dble(coor(i,1))-dble(Nl(1)-1)*nint(dble(coor(i,1))/dble(Nl(1)-1))
    !y=dble(coor(i,2))-dble(Nl(2)-1)*nint(dble(coor(i,2))/dble(Nl(2)-1))
    write(2,*)coor(i,1)-1,coor(i,2)-1
  enddo
  close(2)
endif

end subroutine setnumber



subroutine sethopTBC_multi_allopen
use param
use model_param
use lattice_param
implicit none
integer::i,orb,j,k,ntemp,den
integer::Nh

Nh=0
do i=1,Nbravais,1

!inside triangle hopping
!
!      py
!      |  +
!      d----px
!



  Nh=Nh+1
  sit(Nh,1)=i
  sit(Nh,2)=i+Nbravais     !d --> px
  hopt(Nh)=-tpd


  Nh=Nh+1
  sit(Nh,1)=i+Nbravais
  sit(Nh,2)=i
  hopt(Nh)=-tpd            !d <-- px

  Nh=Nh+1
  sit(Nh,1)=i
  sit(Nh,2)=i+2*Nbravais   !d --> py
  hopt(Nh)=+tpd


  Nh=Nh+1
  sit(Nh,1)=i+2*Nbravais
  sit(Nh,2)=i
  hopt(Nh)=+tpd            !d <-- py

  Nh=Nh+1
  sit(Nh,1)=i+Nbravais
  sit(Nh,2)=i+2*Nbravais   !px --> py
  hopt(Nh)=+tpp


  Nh=Nh+1
  sit(Nh,1)=i+2*Nbravais
  sit(Nh,2)=i+Nbravais
  hopt(Nh)=+tpp            !px <-- py


!inter triangles hopping

!inter triangles hopping

  do j=1,Dimen,1

    den=1
    do k=1,j-1,1
     den=den*Nl(k)
    end do

    if(coor(i,j).EQ.Nl(j)) then
      ntemp=(1-Nl(j))*den+i
    else
      ntemp=i+den
    endif

    Nh=Nh+1
    sit(Nh,1)=i+j*Nbravais
    sit(Nh,2)=ntemp            !px or py ---> d
!    write(*,*)'px or py ---> d',Nh,sit(Nh,1),sit(Nh,2)
    hopt(Nh)=-tpd*(-1)**(j)*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
    if(j.eq.1.and.coor(i,j).EQ.Nl(j)) then
      hopt(Nh)=zero
    endif
    if(j.eq.2.and.coor(i,j).EQ.Nl(j)) then
      hopt(Nh)=zero
    endif

    if(j.eq.1)then
      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      sit(Nh,2)=ntemp+2*Nbravais   !px ----> py
!      write(*,*)'px ----> py',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=-tpp*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
      if(j.eq.1.and.coor(i,j).EQ.Nl(j)) then
        hopt(Nh)=zero
      endif
    else
      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais       !py ----> px
      sit(Nh,2)=ntemp+Nbravais
!      write(*,*)'py ----> px',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=-tpp*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
      if(j.eq.2.and.coor(i,j).EQ.Nl(j)) then
        hopt(Nh)=zero
      endif

      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais
      if(mod(ntemp-1,Nl(1)).eq.0)then
       sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
      else
       sit(Nh,2)=ntemp+Nbravais-1   !py ----> px
      endif
!      write(*,*)'py ----> px, strange',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=tpp*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))
      if(j.eq.2.and.coor(i,1).EQ.1) then
        hopt(Nh)=zero
      endif
      if(j.eq.2.and.coor(i,j).EQ.Nl(j)) then
        hopt(Nh)=zero
      endif
    endif

    if(coor(i,j).EQ.1) then
      ntemp=(Nl(j)-1)*den+i
    else
      ntemp=i-den
    endif

    Nh=Nh+1
    sit(Nh,1)=i
    sit(Nh,2)=ntemp+j*Nbravais   !dx ---> px or py 
    hopt(Nh)=-tpd*(-1)**(j)*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
    if(j.eq.1.and.coor(i,j).EQ.1) then
      hopt(Nh)=zero
    endif
    if(j.eq.2.and.coor(i,j).EQ.1) then
      hopt(Nh)=zero
    endif
!    write(*,*)'dx ---> px or py',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)

    if(j.eq.1)then
      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais
      sit(Nh,2)=ntemp+Nbravais   !py ----> px
      hopt(Nh)=-tpp*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
!      write(*,*)'py ----> px ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)
      if(j.eq.1.and.coor(i,j).EQ.1) then
        hopt(Nh)=zero
      endif
    else
      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      sit(Nh,2)=ntemp+2*Nbravais  !px ----> py
      hopt(Nh)=-tpp*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
!      write(*,*)'px ----> py ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)
      if(j.eq.2.and.coor(i,j).EQ.1) then
        hopt(Nh)=zero
      endif

      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      if(mod(ntemp,Nl(1)).eq.0)then
       sit(Nh,2)=ntemp-Nl(1)+1+2*Nbravais
      else
       sit(Nh,2)=ntemp+1+2*Nbravais
      endif
      hopt(Nh)=tpp*exp((0.d0,-1.d0)*(-kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))
!      write(*,*)'px ----> py -- strange ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)
      if(j.eq.2.and.coor(i,1).EQ.Nl(1)) then
        hopt(Nh)=zero
      endif
      if(j.eq.2.and.coor(i,j).EQ.1) then
        hopt(Nh)=zero
      endif
    endif



  enddo


enddo
if(Nh.ne.Nhop) then
  write(*,*) "Something is wrong with Nhop PBC"
  call mystop
end if


end subroutine sethopTBC_multi_allopen






subroutine sethopTBC_multi_open
use param
use model_param
use lattice_param
implicit none
integer::i,orb,j,k,ntemp,den
integer::Nh

Nh=0
do i=1,Nbravais,1

!inside triangle hopping
!
!      py
!      |  +
!      d----px
!



  Nh=Nh+1
  sit(Nh,1)=i
  sit(Nh,2)=i+Nbravais     !d --> px
  hopt(Nh)=-tpd


  Nh=Nh+1
  sit(Nh,1)=i+Nbravais
  sit(Nh,2)=i
  hopt(Nh)=-tpd            !d <-- px

  Nh=Nh+1
  sit(Nh,1)=i
  sit(Nh,2)=i+2*Nbravais   !d --> py
  hopt(Nh)=+tpd


  Nh=Nh+1
  sit(Nh,1)=i+2*Nbravais
  sit(Nh,2)=i
  hopt(Nh)=+tpd            !d <-- py

  Nh=Nh+1
  sit(Nh,1)=i+Nbravais
  sit(Nh,2)=i+2*Nbravais   !px --> py
  hopt(Nh)=+tpp


  Nh=Nh+1
  sit(Nh,1)=i+2*Nbravais
  sit(Nh,2)=i+Nbravais
  hopt(Nh)=+tpp            !px <-- py


!inter triangles hopping

  do j=1,Dimen,1

    den=1
    do k=1,j-1,1
     den=den*Nl(k)
    end do

    if(coor(i,j).EQ.Nl(j)) then
      ntemp=(1-Nl(j))*den+i
    else
      ntemp=i+den
    endif

    Nh=Nh+1
    sit(Nh,1)=i+j*Nbravais
    sit(Nh,2)=ntemp            !px or py ---> d
!    write(*,*)'px or py ---> d',Nh,sit(Nh,1),sit(Nh,2)
    hopt(Nh)=-tpd*(-1)**(j)*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
    if(j.eq.2.and.coor(i,j).EQ.Nl(j)) then
      hopt(Nh)=zero
    endif

    if(j.eq.1)then
      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      sit(Nh,2)=ntemp+2*Nbravais   !px ----> py
!      write(*,*)'px ----> py',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=-tpp*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
    else
      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais       !py ----> px
      sit(Nh,2)=ntemp+Nbravais
!      write(*,*)'py ----> px',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=-tpp*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
      if(j.eq.2.and.coor(i,j).EQ.Nl(j)) then
        hopt(Nh)=zero
      endif

      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais
      if(mod(ntemp-1,Nl(1)).eq.0)then
       sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
      else
       sit(Nh,2)=ntemp+Nbravais-1   !py ----> px
      endif
!      write(*,*)'py ----> px, strange',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=tpp*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))
      if(j.eq.2.and.coor(i,j).EQ.Nl(j)) then
        hopt(Nh)=zero
      endif
    endif

    if(coor(i,j).EQ.1) then
      ntemp=(Nl(j)-1)*den+i
    else
      ntemp=i-den
    endif

    Nh=Nh+1
    sit(Nh,1)=i
    sit(Nh,2)=ntemp+j*Nbravais   !dx ---> px or py 
    hopt(Nh)=-tpd*(-1)**(j)*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
    if(j.eq.2.and.coor(i,j).EQ.1) then
      hopt(Nh)=zero
    endif
!    write(*,*)'dx ---> px or py',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)

    if(j.eq.1)then
      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais
      sit(Nh,2)=ntemp+Nbravais   !py ----> px
      hopt(Nh)=-tpp*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
!      write(*,*)'py ----> px ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)
    else
      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      sit(Nh,2)=ntemp+2*Nbravais  !px ----> py
      hopt(Nh)=-tpp*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
!      write(*,*)'px ----> py ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)
      if(j.eq.2.and.coor(i,j).EQ.1) then
        hopt(Nh)=zero
      endif

      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      if(mod(ntemp,Nl(1)).eq.0)then
       sit(Nh,2)=ntemp-Nl(1)+1+2*Nbravais
      else
       sit(Nh,2)=ntemp+1+2*Nbravais
      endif
      hopt(Nh)=tpp*exp((0.d0,-1.d0)*(-kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))
!      write(*,*)'px ----> py -- strange ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)
      if(j.eq.2.and.coor(i,j).EQ.1) then
        hopt(Nh)=zero
      endif
    endif



  enddo


enddo
if(Nh.ne.Nhop) then 
  write(*,*) "Something is wrong with Nhop PBC"
  call mystop
end if


end subroutine sethopTBC_multi_open



!----------------------------------------------------------------------
!we get the hoping matrix by Twist boundary condition,site(:,2),hopt(:)
!----------------------------------------------------------------------
subroutine sethopTBC_multi
use param
use model_param
use lattice_param
implicit none
integer::i,orb,j,k,ntemp,den
integer::Nh

Nh=0
do i=1,Nbravais,1

!inside triangle hopping
!
!      py
!      |  +
!      d----px
!
  Nh=Nh+1
  sit(Nh,1)=i
  sit(Nh,2)=i+Nbravais     !d --> px
  hopt(Nh)=-tpd
 
  
  Nh=Nh+1
  sit(Nh,1)=i+Nbravais
  sit(Nh,2)=i
  hopt(Nh)=-tpd            !d <-- px

  Nh=Nh+1
  sit(Nh,1)=i
  sit(Nh,2)=i+2*Nbravais   !d --> py
  hopt(Nh)=+tpd


  Nh=Nh+1
  sit(Nh,1)=i+2*Nbravais
  sit(Nh,2)=i
  hopt(Nh)=+tpd            !d <-- py

  Nh=Nh+1
  sit(Nh,1)=i+Nbravais
  sit(Nh,2)=i+2*Nbravais   !px --> py   
  hopt(Nh)=+tpp

!  write(*,*)'**********',i,Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)


  Nh=Nh+1
  sit(Nh,1)=i+2*Nbravais
  sit(Nh,2)=i+Nbravais
  hopt(Nh)=+tpp            !px <-- py

!  write(*,*)'**********',i,Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)
  

!inter triangles hopping

    
  do j=1,Dimen,1

    den=1
    do k=1,j-1,1
     den=den*Nl(k)
    end do

    if(coor(i,j).EQ.Nl(j)) then
      ntemp=(1-Nl(j))*den+i
    else
      ntemp=i+den
    endif

    Nh=Nh+1
    sit(Nh,1)=i+j*Nbravais
    sit(Nh,2)=ntemp            !px or py ---> d
!    write(*,*)'px or py ---> d',Nh,sit(Nh,1),sit(Nh,2)
    hopt(Nh)=-tpd*(-1)**(j)*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

    if(j.eq.1)then
      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      sit(Nh,2)=ntemp+2*Nbravais   !px ----> py
!      write(*,*)'px ----> py',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=-tpp*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
    else
      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais       !py ----> px
      sit(Nh,2)=ntemp+Nbravais
!      write(*,*)'py ----> px',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=-tpp*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais
      if(mod(ntemp-1,Nl(1)).eq.0)then
       sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
      else
       sit(Nh,2)=ntemp+Nbravais-1   !py ----> px
      endif
!      write(*,*)'py ----> px, strange',Nh,sit(Nh,1),sit(Nh,2)
      hopt(Nh)=tpp*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))
    endif

    if(coor(i,j).EQ.1) then
      ntemp=(Nl(j)-1)*den+i
    else
      ntemp=i-den
    endif

    Nh=Nh+1
    sit(Nh,1)=i
    sit(Nh,2)=ntemp+j*Nbravais   !dx ---> px or py 
    hopt(Nh)=-tpd*(-1)**(j)*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
!    write(*,*)'dx ---> px or py',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)

    if(j.eq.1)then
      Nh=Nh+1
      sit(Nh,1)=i+2*Nbravais
      sit(Nh,2)=ntemp+Nbravais   !py ----> px
      hopt(Nh)=-tpp*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
!      write(*,*)'py ----> px ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)
    else
      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      sit(Nh,2)=ntemp+2*Nbravais  !px ----> py
      hopt(Nh)=-tpp*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
!      write(*,*)'px ----> py ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)

      Nh=Nh+1
      sit(Nh,1)=i+Nbravais
      if(mod(ntemp,Nl(1)).eq.0)then
       sit(Nh,2)=ntemp-Nl(1)+1+2*Nbravais
      else
       sit(Nh,2)=ntemp+1+2*Nbravais
      endif
      hopt(Nh)=tpp*exp((0.d0,-1.d0)*(-kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))
!      write(*,*)'px ----> py -- strange ',Nh,sit(Nh,1),sit(Nh,2),hopt(Nh)

    endif



  enddo
enddo
if(Nh.ne.Nhop) then
  write(*,*) "Something is wrong with Nhop PBC"
  call mystop
end if


end subroutine sethopTBC_multi


subroutine sethopTBC_open
use param
use model_param
use lattice_param
implicit none
integer::i,j,k,ntemp,den
integer::Nh


Nh=0
do i=1,Nsite,1
   do j=1,Dimen,1

      if(coor(i,j).EQ.1) then

          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=(Nl(j)-1)*den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
          if(j.eq.2) then
            hopt(Nh)=zero
          endif

          ntemp=den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))


      else if(coor(i,j).EQ.Nl(j)) then

          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=(1-Nl(j))*den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
          if(j.eq.2) then
            hopt(Nh)=zero
          endif

          ntemp=i-den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

      else

          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=i+den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

          ntemp=i-den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

      end if

   end do
end do
if(Nh.ne.Nhop) then
  write(*,*) "Something is wrong with Nhop PBC"
  call mystop
end if
        

end subroutine sethopTBC_open


subroutine sethopTBC_allopen
use param
use model_param
use lattice_param
implicit none
integer::i,j,k,ntemp,den
integer::Nh


Nh=0
do i=1,Nsite,1
   do j=1,Dimen,1

      if(coor(i,j).EQ.1) then

          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=(Nl(j)-1)*den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
          if(j.eq.2) then
            hopt(Nh)=zero
          endif
          if(j.eq.1) then
            hopt(Nh)=zero
          endif

          ntemp=den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))


      else if(coor(i,j).EQ.Nl(j)) then

          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=(1-Nl(j))*den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
          if(j.eq.2) then
            hopt(Nh)=zero
          endif
          if(j.eq.1) then
            hopt(Nh)=zero
          endif

          ntemp=i-den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

      else

          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=i+den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

          ntemp=i-den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

      end if

   end do
end do
if(Nh.ne.Nhop) then
  write(*,*) "Something is wrong with Nhop PBC"
  call mystop
end if
        

end subroutine sethopTBC_allopen




!----------------------------------------------------------------------
!we get the hoping matrix by Twist boundary condition,site(:,2),hopt(:)
!----------------------------------------------------------------------
subroutine sethopTBC
use param
use model_param
use lattice_param
implicit none
integer::i,j,k,ntemp,den
integer::Nh

Nh=0
do i=1,Nsite,1
   do j=1,Dimen,1

      if(coor(i,j).EQ.1) then

          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=(Nl(j)-1)*den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

          ntemp=den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))          
        
      else if(coor(i,j).EQ.Nl(j)) then
          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=(1-Nl(j))*den+i
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

          ntemp=i-den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

      else
          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do

          ntemp=i+den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

          ntemp=i-den
          Nh=Nh+1
          sit(Nh,1)=i
          sit(Nh,2)=ntemp
          hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))        
      end if

   end do
end do
if(Nh.ne.Nhop) then
  write(*,*) "Something is wrong with Nhop PBC"
  call mystop
end if
end subroutine sethopTBC


!-------------------------------------------
!get the shift matrix in different direction
!-------------------------------------------
subroutine set_Tmatrix
use lattice_param
implicit none
integer,external::bound
integer::i,j,k
integer::NT(Dimen)
integer::tmp,dem
do i=1,Nbravais,1
   do j=1,Dimen,1

      tmp=1;dem=1
      do k=1,Dimen,1
         NT(k)=coor(i,k)
         if(k.eq.j) then
           NT(k)=bound(NT(k)+1,Nl(k))
         end if

         tmp=tmp+(NT(k)-1)*dem
         dem=dem*Nl(k)
      end do

      Tmatrix(i,j)=tmp

      if(tmp.GT.Nsite.OR.tmp.LT.1) then
        write(*,*) "Something is wrong with the Tmatrix",tmp
        call mystop
      end if
   end do
end do
end subroutine set_Tmatrix


!--------------------------------------------------------
!It is a bound function use to fit the boundary condtion.
!--------------------------------------------------------
integer function bound(i,ni)
integer::i,ni
bound=i
if(i.GT.ni) bound=i-ni
if(i.LT.1) bound=i+ni
!DEB
!write(*,*)
!write(*,*)'I am in bound '
!write(*,*)'i,ni ',i,ni
!write(*,*)'bound ',bound
!write(*,*)
if(bound.GT.ni.OR.bound.LT.1) write(*,*) "Something is wrong is the boundary condition"
end function


!-------------------------------------------------------
!give in coord(Dimen), give out the number label 1~Nsite
!-------------------------------------------------------
integer function latt_label(coord)
use lattice_param
implicit none
integer,intent(IN)::coord(1:Dimen)
integer::den,i

latt_label=1
den=1
do i=1,Dimen,1
   latt_label=latt_label+(coord(i)-1)*den
   den=den*Nl(i)
end do
!DEB
!write(*,*)'latt_label '
!write(*,*)'coord ',coord
!write(*,*)'latt_label ',latt_label
!write(*,*)
if(latt_label.LT.1.OR.latt_label.GT.Nsite) then
  write(*,*) "Something is wrong with latt_label output:",latt_label
  call mystop
end if
end function


!test
!subroutine test()
!use lattice_param
!implicit none
!integer,external::latt_label
!integer,external::bound
!integer::i,j,m,n
!integer::cc(1:Dimen),ctmp
!do i=1,Nsite,1
!   do j=1,Nsite,1
!
!      do m=1,Dimen,1
!         ctmp=coor(j,m)-coor(i,m)+1
!         cc(m)=bound(ctmp,Nl(m))
!      end do
!      n=latt_label(cc(1:Dimen))
!      write(*,*) i,j,n;pause
!   end do
!end do
!write(*,*) coor(1,1:Dimen)
!write(*,*) latt_label(coor(1,1:Dimen))
!call mystop
!end subroutine test
!end test
