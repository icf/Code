!-------------------------------------------------
!This subroutine measure the correlation functions
!-------------------------------------------------
subroutine meas_corr()
use param
use io_module
implicit none
complex(kind=8),external::zdotc
complex(kind=8)::Kinm,Vonm

!--------------------------------------------------------------
!for measure (Si^z)(Sj^z) and (Di^+)(Dj)+DiDj^,Si^+Sj^-+Si^-Sj+
!--------------------------------------------------------------
integer::i,j,k
complex(kind=8)::corrSz,corrSd,corrDd

call getbasename()
call getcorrSzname()
call getcorrDdname()
call getcorrSdname()
call openUnit(corrSzname,66,'R')
call openUnit(corrDdname,77,'R')
call openUnit(corrSdname,88,'R')
!do j=1,Nlattice,1
   j=1
   do k=1,Nlattice,1
      call SjSk_towf(eigenstate(1,Nei),Phi(1,1),j,k)
      corrSz=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)
      call DjDk_towf(eigenstate(1,Nei),Phi(1,1),j,k)
      corrDd=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)
      call SjmSkd_towf(eigenstate(1,Nei),Phi(1,1),j,k)
      corrSd=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)
      write(66,'(2I3,f15.8)') j,k,dble(corrSz)
      write(77,'(2I3,f15.8)') j,k,dble(corrDd)
      write(88,'(2I3,f15.8)') j,k,dble(corrSd)
   end do
!end do
close(66)
close(77)
close(88)


!Make fourier transformation of Sz
call fourierSz()
end subroutine meas_corr


!---------------------------------
!Make fourier transformation of Sz
!---------------------------------
subroutine fourierSz()
use param
use io_module
implicit none
integer::i,j,k
real(kind=8)::sz(Nlattice)
complex(kind=8)::sk(Nlattice)
integer::kl(Dimen)
real(kind=8)::ph

call getbasename()
call getcorrSzname()
call getcorrSkname()

call openUnit(corrSzname,66,'O')
do i=1,Nlattice,1
   read(66,'(2I3,f15.8)') j,k,sz(i)
end do
close(66)

call openUnit(corrSkname,77,'R')

sk=zero
do k=1,Nlattice,1
   do i=1,Dimen,1
      kl(i)=coor(k,i)-1
   end do

   do i=1,Nlattice,1
      ph=0.d0
      do j=1,Dimen,1
         ph=ph+2.d0*Pi*dble(kl(j)*(coor(i,j)-1))/dble(Nl(j))
      end do
      sk(k)=sk(k)+exp(Xi*ph)*sz(i)
   end do

   write(77,*) kl
   write(77,*) sk(k)
end do


close(77)
end subroutine fourierSz



!-------------------------------------------------------------
!-------------------------------------------------------------
!S1zSjz to the old wave function and get the new wave function
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine SjSk_towf(coeo,coen,m,n)
use param
implicit none
complex(kind=8),intent(in)::coeo(Nhilbert)
complex(kind=8),intent(out)::coen(Nhilbert)
integer,intent(in)::m,n
integer::ci,cj,coe
integer(kind=8)::num,numn
integer::i,j,k
integer(kind=8)::cyci,cycj,cycm  !Use to define the next step in cycle.

coen=zero
cycm=Binomial(Nup,Nlattice)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(num,i,j,k,ci,cj,coe,numn,cyci,cycj)
!$OMP DO SCHEDULE(GUIDED,chunk)
do num=1,Nhilbert,1

   cycj=(num-1)/cycm+1
   cyci=num-(cycj-1)*cycm

   numn=num
   coe=(Hcoe_up(cyci,m,m)-Hcoe_dn(cycj,m,m))*(Hcoe_up(cyci,n,n)-Hcoe_dn(cycj,n,n)) 
   coen(num)=coen(num)+coeo(numn)*0.25d0*dcmplx(coe)
end do
!$OMP END DO
!$OMP END PARALLEL

!if(cyci.NE.Nhup.OR.(cycj).NE.Nhdn) then
!  write(*,*) "Something is wrong in Htowf"
!end if
end subroutine SjSk_towf


!------------------------------------------------------------------
!------------------------------------------------------------------
!D1dDj+D1Djd to the old wave function and get the new wave function
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine DjDk_towf(coeo,coen,m,n)
use param
implicit none
complex(kind=8),intent(in)::coeo(Nhilbert)
complex(kind=8),intent(out)::coen(Nhilbert)
integer,intent(in)::m,n
integer::ci,cj,coe
integer(kind=8)::num,numn
integer::i,j,k
integer(kind=8)::cyci,cycj,cycm  !Use to define the next step in cycle.

coen=zero
cycm=Binomial(Nup,Nlattice)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(num,i,j,k,ci,cj,coe,numn,cyci,cycj)
!$OMP DO SCHEDULE(GUIDED,chunk)
do num=1,Nhilbert,1

   cycj=(num-1)/cycm+1
   cyci=num-(cycj-1)*cycm

   !DmdDn
   numn=Hm_up(cyci,m,n)+(Hm_dn(cycj,m,n)-1)*cycm
   coe=Hcoe_up(cyci,m,n)*Hcoe_dn(cycj,m,n)
   coen(num)=coen(num)+coeo(numn)*dcmplx(coe)

   !DmDnd
   if(m.eq.n) then
     numn=num
     coe=(1-Hcoe_up(cyci,m,m))*(1-Hcoe_dn(cycj,m,m))
   else
     numn=Hm_up(cyci,n,m)+(Hm_dn(cycj,n,m)-1)*cycm
     coe=Hcoe_up(cyci,n,m)*Hcoe_dn(cycj,n,m)
   end if
   coen(num)=coen(num)+coeo(numn)*dcmplx(coe)

end do
!$OMP END DO
!$OMP END PARALLEL

!if(cyci.NE.Nhup.OR.(cycj).NE.Nhdn) then
!  write(*,*) "Something is wrong in Htowf"
!end if
end subroutine DjDk_towf

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!S1dSjm+S1mS1d to the old wave function and get the new wave function
!--------------------------------------------------------------------
!--------------------------------------------------------------------
subroutine SjmSkd_towf(coeo,coen,m,n)
use param
implicit none
complex(kind=8),intent(in)::coeo(Nhilbert)
complex(kind=8),intent(out)::coen(Nhilbert)
integer,intent(in)::m,n
integer::ci,cj,coe
integer(kind=8)::num,numn
integer::i,j,k
integer(kind=8)::cyci,cycj,cycm  !Use to define the next step in cycle.

coen=zero
cycm=Binomial(Nup,Nlattice)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(num,i,j,k,ci,cj,coe,numn,cyci,cycj)
!$OMP DO SCHEDULE(GUIDED,chunk)
do num=1,Nhilbert,1

   cycj=(num-1)/cycm+1
   cyci=num-(cycj-1)*cycm

   !Sm+Sn-
   if(m.eq.n) then
     numn=num
     coe=Hcoe_up(cyci,m,n)*(1-Hcoe_dn(cycj,n,m))
   else
     numn=Hm_up(cyci,m,n)+(Hm_dn(cycj,n,m)-1)*cycm
     coe=-Hcoe_up(cyci,m,n)*Hcoe_dn(cycj,n,m)
   end if
   coen(num)=coen(num)+coeo(numn)*dcmplx(coe)

  !Sm-Sn+
   if(m.eq.n) then
     numn=num
     coe=(1-Hcoe_up(cyci,n,m))*Hcoe_dn(cycj,m,n)
   else
     numn=Hm_up(cyci,n,m)+(Hm_dn(cycj,m,n)-1)*cycm
     coe=-Hcoe_up(cyci,n,m)*Hcoe_dn(cycj,m,n)
   end if
   coen(num)=coen(num)+coeo(numn)*dcmplx(coe)

end do
!$OMP END DO
!$OMP END PARALLEL

!if(cyci.NE.Nhup.OR.(cycj).NE.Nhdn) then
!  write(*,*) "Something is wrong in Htowf"
!end if
end subroutine SjmSkd_towf
