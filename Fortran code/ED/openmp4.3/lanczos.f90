!This page contains the lanczos method of the model,it can be used to calculate
!different model we only need to change Htowf subroutine depend on different
!model.


!1.lanczos:the main lanczos subroutine,it has the lanczos cycle,measure,write,and stop parts.
!2.initial:initial the lanczos wave funtion Phi(1:Nhilbert,1)
!3.lanmatrix:lanmatrix() get the lanczos matrix,put the diagnal element in dig(LanM)
!            put the offdiagnal element in offdig(LanM)
!4lanmatrix_loop(j,i1,i2,i3) it is the loop in lanmatrix subrutine, it can be
!call by subroutine newphi to get the new phi(:,i2) in eigenstate(:,Nei)
!5.Htowf:calculate the Hamiltonian to wave funtion,from coeo to coen.
!6.schmidt:enforce orthogonoliaze the phi(1:Nhilbert,k) with the eigenstate(1:Nhilbet,1:Nei)
!          and also phi(1:Nhilbert,k) with phi(1:Nhilbert,1-k) only when NlanM
!          is larger than 100 it is needed;or we can use schmidt2.
!7.schmidt_ab:Orthogonalize the Phi_b with phi_a
!8.schmidt_eigs(phi_orth):enforce orthogonoliaze the phi_orth with the eigenstate(1:Nhilbet,1:Nei)
!9.diglan:digonal the lanczos matrix
!10.newPhi:get the new phi(1:Nhilbert,1) depend on the diglan subroutine.
!11.wsave:put the dig(1) to eigenvalue(Nei) and phi(:,1) to eigenstate(:,Nei)
!   write eigenstate(:,Nei) to the hard disk,eigenvalue is write in meas subroutine
!12.resetLanM:rest LanM to (LanMmin,LanMmax-1) by lanadd or by
!             random number generator when lanM is larger than lanMmax
!13.mintomax:get the LanM=[LanMmin,LanMmax-1] randomly
!------------------------------------------------------------------------------------------
!lanczos: the main lanczos subroutine,it has the lanczos cycle,measure,write,and stop part.
!------------------------------------------------------------------------------------------
subroutine lanczos
use param
use rand_num
use io_module
use timing_module
implicit none
!integer(kind=8)::i
!complex(kind=8)::coe

call getbasename()
call geteigvname()
call openUnit(eigvname,30,'R')

Nei=1
Nlanmatrix=1
LanMmem=LanM
varE=2.d0**31
waveconv=.false.
!---------------------------------------------
!initial the phi0 of the lanczos wave function
!---------------------------------------------
call initial()

do
  !--------------------------------------------------------------
  !get the lanczos matrix,put into dig(LanM) and offdiag(LanMoff)
  !--------------------------------------------------------------
  call lanmatrix()
  write(*,*) " "
  write(*,*) " "
  write(*,*) "get lanmatrix",Nlan
  write(*,*) "the number of dig lanmatrix",Nlanmatrix
  call PrintTiming()
  Nlanmatrix=Nlanmatrix+1

  if(Nlan.EQ.1) then
    write(*,*) " "
    write(*,*) " "
    write(*,*) "Get the ",Nei,"eigenstate------done"
    call PrintTiming()
    write(*,*) " "
    write(*,*) " "
    write(*,*) " "


    call wsave()
    call meas()
    Nei=Nei+1
    Nlanmatrix=1
    LanM=LanMmem
    varE=2.d0**31
    waveconv=.false.

    if(Nei.GT.Nexcit) then
      write(*,*) "Get the whole eigenvalue and eigenstate------done"
      call PrintTiming()

      close(30)
      !call mystop
      return
    end if

    call initial()
  else 
    !---------------------------
    !digonal the lanczos matrix
    !---------------------------
     call diglan()


    !---------------------
    !get the new Phi(:,i2)
    !---------------------
     call newPhi()

     write(*,*) " "
     write(*,*) " "
     write(*,*) "get the new phi(:,i2) in eigenstate(:,Nei)"
     call PrintTiming()

    !-------------------------------------------------------------------------
    !reset the LanM to help converge--Sometime fix LanM will cost too many cpu
    !time.---add by shihao 2010.11.11 
    !-------------------------------------------------------------------------
    if(Nlanmatrix.GT.Nlanmatrix_set) then
       call resetLanM()
       write(*,*) " "
       write(*,*) " "
       write(*,*) "Reset the LanM, the new LanM is:",LanM
       call PrintTiming()
    end if
  end if

end do

end subroutine lanczos



!--------------------------------
!initial the lanczos wave funtion
!--------------------------------
subroutine initial
use param
use rand_num
implicit none
integer(kind=8)::i,j
real(kind=8)::a,b

do i=1,Nhilbert,1
   a=rndm()-0.5d0
   b=rndm()-0.5d0
   eigenstate(i,Nei)=dcmplx(a,b)
end do

!-------------------------------------------------------------------
!If calculate the excited state, orthogonalize with low energy state
!-------------------------------------------------------------------
call schmidt_eigs(eigenstate(1:Nhilbert,Nei))


!----------------------------
!Project S^2 to wave function
!----------------------------
if(PS.eq.'Y') then
  call S_con(eigenstate(1:Nhilbert,Nei))
else if(PS.eq.'N') then
else
  write(*,*) "Something is wrong with PS inupt"
end if

!--------------------------
!Project T to wave function
!--------------------------
if(PK.eq.'Y') then
  call k_con(eigenstate(1:Nhilbert,Nei))
else if(PK.eq.'N') then
else
  write(*,*) "Something is wrong with PK inupt"
end if

call norm_wave(eigenstate(1:Nhilbert,Nei))

!write(*,*) Phi(1,1),Phi(5,1),Phi(Nhilbert,1);pause 

!sqrtnorm=sqrt(norm)
!do i=1,Nhilbert,1
!   temp=Phi(i,1)
!   Phi(i,1)=temp/sqrtnorm
!end do
!write(*,*) Phi(1,1),Phi(5,1),Phi(Nhilbert,1);pause 

!norm=(0.d0,0.d0)
!do i=1,Nhilbert,1
!   norm=norm+Abs(Phi(i,1))**2
!end do
!write(*,*) norm
!pause


!test
!do j=1,9,1
!   do i=1,Nhilbert,1
!      a=rndm()-0.5d0
!      b=rndm()-0.5d0
!      Phi(i,j)=dcmplx(a,b)
!   end do
!
!   if(j.NE.1) then 
!     call schmidt(j)
!   end if
!   
!   norm=(0.d0,0.d0)
!   do i=1,Nhilbert,1
!      norm=norm+Abs(Phi(i,j))**2
!   end do   
!   sqrtnorm=sqrt(norm) 
!
!   do i=1,Nhilbert,1
!      temp=Phi(i,j)
!      Phi(i,j)=temp/sqrtnorm
!   end do
!end do
!
!do j=1,6,1
!   norm=(0.d0,0.d0)
!   do i=1,Nhilbert,1
!      norm=norm+conjg(Phi(i,j))*Phi(i,6)
!   end do
!   write(*,*) norm,j,Nhilbert
!   pause
!end do
!call mystop
end subroutine initial





!---------------------------------------------------
!lanmatrix() get the lanczos matrix,put the diagnal
!element in dig(LanM),put the offdiagnal element in
!offdig(LanM)
!---------------------------------------------------
subroutine lanmatrix
use param
implicit none
integer::j
integer::i1,i2,i3

dig=0.d0;offdig=0.d0
!Set the phi(:,3) first
i1=1;i2=2;i3=3
call zcopy(Nhilbert,eigenstate(1,Nei),1,Phi(1,i2),1)

Nlan=LanM
do j=1,LanM,1
  call lanmatrix_loop(j,i1,i2,i3)
  if(j.eq.1) then
    write(*,*) ""
    write(*,*) "Get the variational energy:",dig(j)
    write(*,*) ""
    write(*,*) ""

    if(conv.eq.'E') then
      if(.not.waveconv) then
        if(dig(j).GT.varE) then
          write(*,*) "The variational energy is going up"
          write(*,*) "Energy converge has been to its numerical limit."
          write(*,*) "We need to stop the lanczos due to the converge condition"
          Nlan=j
        end if
      end if
    else if(conv.eq.'W') then
    else
      write(*,*) "Something is wrong with the conv input: E OR W"
    end if

    varE=dig(j)
  end if

  if(Nlan.eq.j) exit
  i1=i1+1;if(i1.GT.3) i1=1
  i2=i2+1;if(i2.GT.3) i2=1
  i3=i3+1;if(i3.GT.3) i3=1
end do

end subroutine lanmatrix


!---------------------------------------------------
!lanmatrix() get the lanczos matrix,put the diagnal
!element in dig(LanM),put the offdiagnal element in
!offdig(LanM), the loop of lanmatrix. 
!---------------------------------------------------
subroutine lanmatrix_loop(j,i1,i2,i3)
use param
implicit none
integer,intent(IN)::j,i1,i2,i3
complex(kind=8),external::zdotc
real(kind=8),external::dznrm2
complex(kind=8)::a,b


!H-->Phi1
call Htowf(Phi(1,i2),Phi(1,i3))

!get a1----dig(1)
b=zdotc(Nhilbert,Phi(1,i2),1,Phi(1,i3),1) 
if(AIMAG(b).GT.1.D-10) then
  write(*,*) "Something is wrong with the Lanczos a"
  call mystop
end if
dig(j)=Real(b)

!For the LanM, we only need dig(LanM)
if(j.EQ.LanM) then
  return
end if


! get b2*Phi2
b=(-1.d0,0.d0)*dcmplx(dig(j))
call zaxpy(Nhilbert,b,Phi(1,i2),1,Phi(1,i3),1)
if(j.GT.1) then
   b=(-1.d0,0.d0)*dcmplx(offdig(j))
   call zaxpy(Nhilbert,b,Phi(1,i1),1,Phi(1,i3),1)
end if
offdig(j+1)=dznrm2(Nhilbert,Phi(1,i3),1)


!-------------------------------------------------------- 
! make schmidt orthogonalize,the eigenstate orthogonalize
!-------------------------------------------------------- 
if(offdig(j+1).LT.lim.AND.offdig(j+1).GT.lit) then

   if(j.NE.1) then
     write(*,*) ""
     write(*,*) ""
     write(*,*) "It maybe full hilbert space at j+1 wavefunction:",j+1,offdig(j+1)
     write(*,*) "Reduce the size of LanM will be better for the algorithm."
     write(*,*) ""
     write(*,*) ""
     Nlan=j
     return
   else

     write(*,*) "Project for j+1:",j+1,offdig(j+1)

     !Normalize the Phi(:,i3):do this first to reduce the numerical error
     a=1.d0/dcmplx(offdig(j+1))
     call zscal(Nhilbert,a,Phi(1,i3),1)


     !GS Phi(:,i3) with Phi(:,i2)
     call schmidt_ab(Phi(1,i2),Phi(1,i3))

     !GS Phi(:,i3) with Eigenstate
     if(Nei.GT.1) then
       call schmidt_eigs(Phi(1:Nhilbert,i3)) 
     end if

     !Project Phi(:,i3) to S0
     if(PS.eq.'Y') then
       call S_con(Phi(1:Nhilbert,i3))
     else if(PS.eq.'N') then
     else
       write(*,*) "Something is wrong with PS inupt"
     end if

     !Project Phi(:,i3) to K0
     if(PK.eq.'Y') then
       call k_con(Phi(1:Nhilbert,i3))
     else if(PK.eq.'N') then
     else
       write(*,*) "Something is wrong with PK inupt"
     end if

     b=dcmplx(dznrm2(Nhilbert,Phi(1:Nhilbert,i3),1))
     offdig(j+1)=offdig(j+1)*dble(b)
     a=1.d0/b
     call zscal(Nhilbert,a,Phi(1,i3),1)

   end if

else if(offdig((j+1)).LT.lit) then
   Nlan=j
   if(j.eq.1) then
     write(*,*) "wave function has been converged."
     waveconv=.True.
   end if
   return
else
   !get Phi(:,i3)
   a=1.d0/dcmplx(offdig(j+1))
   call zscal(Nhilbert,a,Phi(1,i3),1)
end if


end subroutine lanmatrix_loop





!--------------------------------------------
!it calculate the Hamiltonian to wave funtion
!--------------------------------------------
subroutine Htowf(coeo,coen)
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
end subroutine Htowf


!-----------------------------------------------------------------
!Orthogonalize the Phi from 1 to k-1 to k also with the eigenstate
!-----------------------------------------------------------------
!subroutine schmidt(n)
!use param
!implicit none
!integer::n,j
!integer(kind=8)::i
!complex(kind=8)::coe(n-1)
!complex(kind=8)::temp
!complex(kind=8)::eigen(Nei-1)
!complex(kind=8),external::zdotc
!
!
!call zgemv('C',Nhilbert,n-1,(1.d0,0.d0),Phi(1,1),Nhilbert,Phi(1,n),1,(0.d0,0.d0),coe(1),1)
!!coe=(0.d0,0.d0)
!!do j=1,n-1,1
!!   coe(j)=zdotc(Nhilbert,Phi(1,j),1,Phi(1,n),1)
!!end do
!
!if(Nei.GT.1) then
!  call zgemv('C',Nhilbert,Nei-1,(1.d0,0.d0),eigenstate(1,1),Nhilbert,Phi(1,n),1,(0.d0,0.d0),eigen(1),1)
!end if
!!if(Nei.GT.1) then
!!  eigen=(0.d0,0.d0)
!!  do j=1,Nei-1,1
!!     eigen(j)=zdotc(Nhilbert,eigenstate(1,j),1,Phi(1,n),1)
!!  end do
!!  write(*,*) eigen;pause
!!end if
!
!do j=1,n-1,1
!   temp=-1.d0*coe(j)
!   call zaxpy(Nhilbert,temp,Phi(1,j),1,Phi(1,n),1)
!end do
!
!
!if(Nei.GT.1) then
!  do j=1,Nei-1,1
!     temp=-1.d0*eigen(j)
!     call zaxpy(Nhilbert,temp,eigenstate(1,j),1,Phi(1,n),1)
!  end do
!end if
!!normalize will be done in subroutine lanmatrix
!! temp=(0.d0,0.d0)
!! do i=1,Nhilbert,1
!!    temp=temp+Abs(Phi(i,n))**2
!! end do
!! offdig(n)=sqrt(Real(temp))
!!
!! do i=1,Nhilbert,1
!!    temp=Phi(i,n)
!!    Phi(i,n)=temp/dcmplx(offdig(n))
!! end do
!
!end subroutine schmidt


!----------------------------------
!Orthogonalize the Phi_b with phi_a
!----------------------------------
subroutine schmidt_ab(phi_a,phi_b)
use param
implicit none
complex(kind=8),intent(IN)::phi_a(Nhilbert)
complex(kind=8),intent(INOUT)::phi_b(Nhilbert)
complex(kind=8)::temp
complex(kind=8)::coe
complex(kind=8),external::zdotc

coe=zdotc(Nhilbert,Phi_a(1),1,Phi_b(1),1)
temp=-1.d0*coe
call zaxpy(Nhilbert,temp,phi_a(1),1,phi_b(1),1)
end subroutine schmidt_ab


!----------------------------------------
!Orthogonalize the Phi(n) with eigenstate
!----------------------------------------
subroutine schmidt_eigs(phi_orth)
use param
implicit none
complex(kind=8),intent(INOUT)::phi_orth(Nhilbert)
complex(kind=8)::temp
complex(kind=8)::eigen(Nei-1)
integer::i,j

if(Nei.GT.1) then
  call zgemv('C',Nhilbert,Nei-1,(1.d0,0.d0),eigenstate(1,1),Nhilbert,phi_orth(1),1,(0.d0,0.d0),eigen(1),1)
  do j=1,Nei-1,1
     temp=-1.d0*eigen(j)
     call zaxpy(Nhilbert,temp,eigenstate(1,j),1,phi_orth(1),1)
  end do
end if
end subroutine schmidt_eigs


!--------------------------
!digonal the lanczos matrix
!--------------------------
subroutine diglan
use param
implicit none
real(kind=8)::a(Nlan),b(Nlan-1)
integer::i,j,k
real(kind=8),allocatable::work(:)
integer::lwork
integer::info


a(1)=dig(1)
do i=2,Nlan,1
   j=i-1
   a(i)=dig(i)
   b(j)=offdig(i)
end do

lwork=max(1,2*Nlan-2)
allocate(work(lwork))


MatrixU=0.d0
call dsteqr('I',Nlan,a,b,MatrixU,LanM,work,info)

deallocate(work)

if(info.NE.0) then
  write(*,*) "info=",info
  write(*,*) "Dsteqr failed to get exact result"
  call mystop
end if

end subroutine diglan


!---------------------
!get the new Phi(:,i2)
!---------------------
subroutine newPhi()
use param
implicit none
integer::i,j,k
complex(kind=8)::Norm
complex(kind=8)::complexU
real(kind=8),external::dznrm2
integer::i1,i2,i3

i1=1;i2=2;i3=3
call zcopy(Nhilbert,eigenstate(1,Nei),1,Phi(1,i2),1)
complexU=dcmplx(MatrixU(1,1))
call zscal(Nhilbert,complexU,eigenstate(1,Nei),1)
Nlan=LanM
do j=1,LanM-1,1

  call lanmatrix_loop(j,i1,i2,i3)
  if(Nlan.eq.j) exit

  complexU=dcmplx(MatrixU(j+1,1))
  call zaxpy(Nhilbert,complexU,Phi(1,i3),1,eigenstate(1,Nei),1)

  i1=i1+1;if(i1.GT.3) i1=1
  i2=i2+1;if(i2.GT.3) i2=1
  i3=i3+1;if(i3.GT.3) i3=1
end do

norm=dcmplx(dznrm2(Nhilbert,eigenstate(1,Nei),1))
if(abs(norm-1.d0).GT.1d-8) then
  write(*,*) "Something is wrong with the new phi norm:",norm
  write(*,*) "Try to add rigour project condition in lanmatrix_loop if Nei>1,or PS PK is Y"
  write(*,*) "Otherwise try to reduce the lanM."
  write(*,*) "The lanczos matrix::"
  write(*,*) "1",dig(1),0.d0
  do i=2,Nlan,1
     write(*,*) i,dig(i),offdig(i)
  end do
  write(*,*) "The eigenvector will is:"
  do i=1,Nlan,1
     write(*,*) i,MatrixU(i,1)
  end do
  call mystop
end if

!-----------------------------------------------
!Orthogonalize with eigenstate for excited state
!-----------------------------------------------
call schmidt_eigs(eigenstate(1:Nhilbert,Nei))

!----------------------------
!Project S^2 to wave function
!----------------------------
if(PS.eq.'Y') then
  call S_con(eigenstate(1:Nhilbert,Nei))
else if(PS.eq.'N') then
else
  write(*,*) "Something is wrong with PS inupt"
end if

!--------------------------
!Project T to wave function
!--------------------------
if(PK.eq.'Y') then
  call k_con(eigenstate(1:Nhilbert,Nei))
else if(PK.eq.'N') then
else
  write(*,*) "Something is wrong with PK inupt"
end if

call norm_wave(eigenstate(1:Nhilbert,Nei))


end subroutine newPhi


!---------------------------------------
!save the eigenvalue and the eigenstates
!---------------------------------------
subroutine wsave()
use param
use io_module
implicit none
integer(kind=8)::i

eigenvalue(Nei)=dig(1)

if(wstate.EQ.1) then
  call geteigsname()
  call openUnit(eigsname,20,'R') 
  do i=1,Nhilbert,1
     write(20,*) Real(eigenstate(i,Nei)),Dimag(eigenstate(i,Nei))
  end do
  close(20)
else if(wstate.EQ.0) then
else
  write(*,*) "Something is wrong with the wstate input.1-->write down,0-->donot write"
  call mystop
end if


end subroutine wsave


!----------------------------------------------------
!add to rest lanM() contains LanM=(LanMmin,LanMmax-1)
!----------------------------------------------------
subroutine resetLanM() 
use param
implicit none

 LanM=LanM+lanadd
 if(LanM.GE.LanMmax) then
    call mintomax()
 end if
 Nlanmatrix=1
 deallocate(dig,offdig,MatrixU)
 allocate(dig(LanM),offdig(LanM+1),MatrixU(LanM,LanM)) 
end subroutine resetLanM
!-----------------------------------------
!get the LanM=[LanMmin,LanMmax-1] randomly
!-----------------------------------------
subroutine mintomax()
use param
use rand_num
implicit none
  LanM=(LanMmax-LanMmin)*rndm(); if(LanM>((LanMmax-LanMmin)-1)) LanM=(LanMmax-LanMmin)-1
  LanM=LanM+LanMmin
end subroutine mintomax
