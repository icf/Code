!This page contains the subroutine need to measure the physical quantity.

!1.meas: measure the energy,kinetic energy,potential energy and write them into the file.
!2.Kintowf:kinectic hamiltonian to old wave function and get the new wave function
!3.Vtowf:potential haliltonian to the old wave function and get the new wave function
!We measure:(Si^z)(Sj^z)
!           (Di^+)(Dj)+DiDj^
!           Si^+Sj^-+Si^-Sj+
!           K
!           V

!-------------------------------
!Measure all the physics we need
!-------------------------------
subroutine meas()
use param
use io_module
implicit none
call meas_Square()
call meas_k()
call meas_corr()
call meas_cicj()
call meas_energy()
end subroutine meas


!--------------------------------------------------------------------------------
!measure the energy,kinetic energy,potential energy and write them into the file.
!--------------------------------------------------------------------------------
subroutine meas_energy()
use param
use io_module
implicit none
complex(kind=8),external::zdotc
complex(kind=8)::Kinm,Vonm

!--------------------------------------------------------------
!for measure (Si^z)(Sj^z) and (Di^+)(Dj)+DiDj^,Si^+Sj^-+Si^-Sj+
!--------------------------------------------------------------
!integer::i,j
!complex(kind=8)::corrSz,corrSd,corrDd


call Kintowf(eigenstate(1,Nei),Phi(1,1))
Kinm=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)


call Vtowf(eigenstate(1,Nei),Phi(1,1))
Vonm=zdotc(Nhilbert,eigenstate(1,Nei),1,Phi(1,1),1)


write(30,*) dble(Kinm),dble(Vonm),Real(eigenvalue(Nei))

write(*,*) "Kinm",Nei,dble(Kinm)
write(*,*) "Vonm",Nei,dble(Vonm)
write(*,*) "Em",Real(eigenvalue(Nei))

end subroutine meas_energy

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!kinectic hamiltonian to old wave function and get the new wave function
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine Kintowf(coeo,coen)
use param
implicit none
complex(kind=8),intent(IN)::coeo(Nhilbert)
complex(kind=8),intent(OUT)::coen(Nhilbert)
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

   do i=1,Nhop,1
      ci=sit(i,1);cj=sit(i,2)
      !For the spin up
      numn=Hm_up(cyci,ci,cj)+(cycj-1)*cycm
      coe=Hcoe_up(cyci,ci,cj)
      coen(num)=coen(num)+coeo(numn)*conjg(hopt(i))*dcmplx(coe)
      !For the spin dn
      numn=cyci+(Hm_dn(cycj,ci,cj)-1)*cycm
      coe=Hcoe_dn(cycj,ci,cj)
      coen(num)=coen(num)+coeo(numn)*conjg(hopt(i))*dcmplx(coe)
   end do

end do
!$OMP END DO
!$OMP END PARALLEL

!if(cyci.NE.Nhup.OR.(cycj).NE.Nhdn) then
!  write(*,*) "Something is wrong in Htowf"
!end if
end subroutine Kintowf

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!potential haliltonian to the old wave function and get the new wave function
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine Vtowf(coeo,coen)
use param
implicit none
complex(kind=8),intent(IN)::coeo(Nhilbert)
complex(kind=8),intent(OUT)::coen(Nhilbert)
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

   do i=1,Nlattice,1
      numn=num
      coe=Hcoe_up(cyci,i,i)*Hcoe_dn(cycj,i,i)
      coen(num)=coen(num)+coeo(numn)*dcmplx(onsitU)*dcmplx(coe)
   end do

end do
!$OMP END DO
!$OMP END PARALLEL

!if(cyci.NE.Nhup.OR.(cycj).NE.Nhdn) then
!  write(*,*) "Something is wrong in Htowf"
!end if
end subroutine Vtowf


