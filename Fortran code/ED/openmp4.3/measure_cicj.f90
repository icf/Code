!----------------------------------------------------------------------------
!This subroutine measure the ci^+cj functions and make fourier transformation
!----------------------------------------------------------------------------
subroutine meas_cicj()
use param
use io_module
implicit none
complex(kind=8),external::zdotc
integer::i,j,k
integer::kl(Dimen)
complex(kind=8)::cijup(Nlattice,Nlattice)
complex(kind=8)::cijdn(Nlattice,Nlattice)
complex(kind=8)::kkup(Nlattice),kkdn(Nlattice)
real(kind=8)::ph

write(*,*) "measure ci+cj"
call getbasename()
call getcoorcicjname()
call openUnit(coorcicjname,99,'R')
do i=1,Nlattice,1
   !do j=1,Nlattice,1 
      j=1
      call cicjup_towf(eigenstate(1,Nei),Phi(1,1),i,j)
      cijup(i,j)=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)
      call cicjdn_towf(eigenstate(1,Nei),Phi(1,1),i,j)
      cijdn(i,j)=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)
   !end do
end do
   do i=1,Nlattice,1
      write(99,*) cijup(i,j)
   end do
   do i=1,Nlattice,1
      write(99,*) cijdn(i,j)
   end do
close(99)

!Make fourier transformation
call getbasename()
call getcoornkname()
call openUnit(coornkname,99,'R')

kkup=zero
kkdn=zero
do k=1,Nlattice,1
   do i=1,Dimen,1
      kl(i)=coor(k,i)-1
   end do
 
   do i=1,Nlattice,1
      ph=0.d0
      do j=1,Dimen,1
         ph=ph+2.d0*Pi*dble(kl(j)*(coor(i,j)-1))/dble(Nl(j))
      end do
      kkup(k)=kkup(k)+exp(Xi*ph)*cijup(i,1)
      kkdn(k)=kkdn(k)+exp(Xi*ph)*cijdn(i,1)
      !write(*,*) exp(Xi*ph);pause
   end do

end do

do k=1,Nlattice,1
   do i=1,Dimen,1
      kl(i)=coor(k,i)-1
   end do
   write(99,*) kl
   write(99,*) kkup(k)
end do
do k=1,Nlattice,1
   do i=1,Dimen,1
      kl(i)=coor(k,i)-1
   end do
   write(99,*) kl
   write(99,*) kkdn(k)
end do

close(99)
!call mystop
end subroutine meas_cicj


!--------------------------------------------------------------
!--------------------------------------------------------------
!cicjup to the  old wave function and get the new wave function
!--------------------------------------------------------------
!--------------------------------------------------------------
subroutine cicjup_towf(coeo,coen,i,j)
use param
implicit none
complex(kind=8),intent(IN)::coeo(Nhilbert)
complex(kind=8),intent(OUT)::coen(Nhilbert)
integer,intent(in)::i,j
integer::ci,cj,coe
integer(kind=8)::num,numn
integer(kind=8)::cyci,cycj,cycm  !Use to define the next step in cycle.

coen=zero
cycm=Binomial(Nup,Nlattice)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(num,ci,cj,coe,numn,cyci,cycj)
!$OMP DO SCHEDULE(GUIDED,chunk)
do num=1,Nhilbert,1

   cycj=(num-1)/cycm+1
   cyci=num-(cycj-1)*cycm

   !We calculate (ci^{+}cj)^+ 
   ci=j;cj=i
   !For the spin up
    numn=Hm_up(cyci,ci,cj)+(cycj-1)*cycm
    coe=Hcoe_up(cyci,ci,cj)
    coen(num)=coen(num)+coeo(numn)*dcmplx(coe)
end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine cicjup_towf


!--------------------------------------------------------------
!--------------------------------------------------------------
!cicjdn to the  old wave function and get the new wave function
!--------------------------------------------------------------
!--------------------------------------------------------------
subroutine cicjdn_towf(coeo,coen,i,j)
use param
implicit none
complex(kind=8),intent(IN)::coeo(Nhilbert)
complex(kind=8),intent(OUT)::coen(Nhilbert)
integer,intent(in)::i,j
integer::ci,cj,coe
integer(kind=8)::num,numn
integer(kind=8)::cyci,cycj,cycm  !Use to define the next step in cycle.

coen=zero
cycm=Binomial(Nup,Nlattice)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(num,ci,cj,coe,numn,cyci,cycj)
!$OMP DO SCHEDULE(GUIDED,chunk)
do num=1,Nhilbert,1

   cycj=(num-1)/cycm+1
   cyci=num-(cycj-1)*cycm

   !We calculate (ci^{+}cj)^+ 
   ci=j;cj=i
   !For the spin dn
   numn=cyci+(Hm_dn(cycj,ci,cj)-1)*cycm
   coe=Hcoe_dn(cycj,ci,cj)
   coen(num)=coen(num)+coeo(numn)*dcmplx(coe)
end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine cicjdn_towf
