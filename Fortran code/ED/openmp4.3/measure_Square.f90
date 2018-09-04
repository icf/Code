!--------------------------------
!This subroutine measure S_Square
!--------------------------------
subroutine meas_Square()
use param
use io_module
implicit none
complex(kind=8),external::zdotc
complex(kind=8)::b
integer(kind=8)::i
integer::j,k

call Square_towf(eigenstate(1,Nei),Phi(1,1))
S_square=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)

write(*,*) "Square S:",S_square
end subroutine meas_Square


!-------------------------------------------------------------
!-------------------------------------------------------------
!S1zSjz to the old wave function and get the new wave function
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine Square_towf(coeo,coen)
use param
implicit none
complex(kind=8),intent(in)::coeo(Nhilbert)
complex(kind=8),intent(out)::coen(Nhilbert)
integer::ci,cj,coe
integer(kind=8)::num,numn
integer::i,j,k,m,n
integer(kind=8)::cyci,cycj,cycm  !Use to define the next step in cycle.

coen=zero
cycm=Binomial(Nup,Nlattice)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(num,i,j,k,m,n,ci,cj,coe,numn,cyci,cycj)
!$OMP DO SCHEDULE(GUIDED,chunk)
do num=1,Nhilbert,1

   cycj=(num-1)/cycm+1
   cyci=num-(cycj-1)*cycm

   do m=1,Nlattice,1
      do n=1,Nlattice,1

        !SmzSnz
         numn=num
         coe=(Hcoe_up(cyci,m,m)-Hcoe_dn(cycj,m,m))*(Hcoe_up(cyci,n,n)-Hcoe_dn(cycj,n,n))
         coen(num)=coen(num)+coeo(numn)*0.25d0*dcmplx(coe)

        !Sm+Sn-
         if(m.eq.n) then
           numn=num
           coe=Hcoe_up(cyci,m,n)*(1-Hcoe_dn(cycj,n,m))
         else
           numn=Hm_up(cyci,m,n)+(Hm_dn(cycj,n,m)-1)*cycm
           coe=-Hcoe_up(cyci,m,n)*Hcoe_dn(cycj,n,m)
         end if
         coen(num)=coen(num)+coeo(numn)*0.5d0*dcmplx(coe)

        !Sm-Sn+
         if(m.eq.n) then
           numn=num
           coe=(1-Hcoe_up(cyci,n,m))*Hcoe_dn(cycj,m,n)
         else
           numn=Hm_up(cyci,n,m)+(Hm_dn(cycj,m,n)-1)*cycm
           coe=-Hcoe_up(cyci,n,m)*Hcoe_dn(cycj,m,n)
         end if
         coen(num)=coen(num)+coeo(numn)*0.5d0*dcmplx(coe)

      end do
   end do

end do
!$OMP END DO
!$OMP END PARALLEL

!if(cyci.NE.Nhup.OR.(cycj).NE.Nhdn) then
!  write(*,*) "Something is wrong in Htowf"
!end if
end subroutine Square_towf
