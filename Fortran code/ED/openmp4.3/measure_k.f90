!--------------------------------
!This subroutine measure momentum
!--------------------------------
subroutine meas_k()
use param
use io_module
implicit none
complex(kind=8),external::zdotc
complex(kind=8)::b,temp
integer(kind=8)::i
integer::j,k
do j=1,Dimen,1

   call k_towf(eigenstate(1,Nei),Phi(1,1),j)

   temp=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)
 
   k_int(j)=log(temp)/Xi
   k_int(j)=-1.d0*k_int(j)*dcmplx(Nl(j)/(2.d0*Pi))

   write(*,*) "k_int point:",j,k_int(j)

end do
end subroutine meas_k


!---------------------------------------------------------------
!---------------------------------------------------------------
!T matrix to the old wave function and get the new wave function
!---------------------------------------------------------------
!---------------------------------------------------------------
subroutine k_towf(coeo,coen,dd)
use param
implicit none
complex(kind=8),intent(in)::coeo(Nhilbert)
complex(kind=8),intent(out)::coen(Nhilbert)
integer,intent(in)::dd
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

   numn=Tm_up(cyci,dd)+(Tm_dn(cycj,dd)-1)*cycm
   coe=Tcoe_up(cyci,dd)*Tcoe_dn(cycj,dd)
   coen(num)=coen(num)+coeo(numn)*dcmplx(coe)

end do
!$OMP END DO
!$OMP END PARALLEL


!if(cyci.NE.Nhup.OR.(cycj).NE.Nhdn) then
!  write(*,*) "Something is wrong in Htowf"
!end if
end subroutine k_towf
