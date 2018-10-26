!---------------------------------------------------------------
!Get the Hm_up,Hm_dn,Tm_up,Tm_dn,Hcoe_up,Hcoe_dn,Tcoe_up,Tcoe_dn
!---------------------------------------------------------------
subroutine get_Hm()
use param
implicit none
integer,allocatable::Nfi(:),dltNfi(:)
integer(kind=8)::num,numn
integer::ci,cj,coe,dd
integer::flag
integer::i,j,k

!----------------------
!Calculate spin up term
!----------------------
allocate(Nfi(Nup),dltNfi(Nup))
do i=1,Nup,1
   Nfi(i)=i
   dltNfi(i)=0
end do
dltNfi(Nup)=Nlattice-Nup

do num=1,Nhup,1
   do ci=1,Nlattice,1
      do cj=1,Nlattice,1
         numn=num
         call cicjBasis(ci,cj,coe,numn,Nfi,Nup)
         Hm_up(num,ci,cj)=numn
         Hcoe_up(num,ci,cj)=coe
      end do
   end do
   do dd=1,Dimen,1
      numn=num
      call TBasis(coe,numn,Nfi,Nup,dd)
      Tm_up(num,dd)=numn
      Tcoe_up(num,dd)=coe
   end do
   
   call next(Nfi,dltNfi,flag,Nup)
end do

if(flag.NE.0) then
  write(*,*) "Something is wrong with the next in Hm up subroutine."
  call mystop
end if
deallocate(Nfi,dltNfi)


!----------------------
!Calculate spin dn term
!----------------------
allocate(Nfi(Ndn),dltNfi(Ndn))
do i=1,Ndn,1
   Nfi(i)=i
   dltNfi(i)=0
end do
dltNfi(Ndn)=Nlattice-Ndn

do num=1,Nhdn,1
   do ci=1,Nlattice,1
      do cj=1,Nlattice,1
         numn=num
         call cicjBasis(ci,cj,coe,numn,Nfi,Ndn)
         Hm_dn(num,ci,cj)=numn
         Hcoe_dn(num,ci,cj)=coe
      end do
   end do
   do dd=1,Dimen,1
      numn=num
      call TBasis(coe,numn,Nfi,Ndn,dd)
      Tm_dn(num,dd)=numn
      Tcoe_dn(num,dd)=coe
   end do

   call next(Nfi,dltNfi,flag,Ndn)
end do

if(flag.NE.0) then
  write(*,*) "Something is wrong with the next in Hm dn subroutine."
  call mystop
end if
deallocate(Nfi,dltNfi)

end subroutine get_Hm


!------------------------------------------------
!subroutine ci+ cj --->Basis with particular spin
!------------------------------------------------
subroutine cicjBasis(ci,cj,coe,numn,Nfi,Nsp)
use param
implicit none
integer,intent(in)::ci,cj,Nsp
integer,intent(in)::Nfi(Nsp)
integer,intent(out)::coe
integer(kind=8),intent(inout)::numn
integer::bin(Nlattice)
integer::i,j,k
integer::si,sj,totals

if(ci.GT.Nlattice.OR.cj.GT.Nlattice) then
  write(*,*) "Both ci and cj should not be larger than the Nlattice."
  call mystop
end if

bin=0
do i=1,Nsp,1
   bin(Nfi(i))=1
end do

if(ci.ne.cj) then
  si=bin(ci);sj=bin(cj)
  if(si.EQ.0.AND.sj.EQ.1) then
    if(ci.GT.cj) then
      totals=0
      do i=cj+1,ci-1,1
         totals=totals+bin(i)
      end do
      coe=(-1)**totals
    else if(ci.LT.cj) then
      totals=0
      do i=ci+1,cj-1,1
         totals=totals+bin(i)
      end do
      coe=(-1)**totals
    else
       write(*,*) 'Something is wrong in ci+cj';call mystop
    end if
    bin(ci)=1;bin(cj)=0
    call get_num_sp(bin,Numn,Nsp)
  else
     coe=0
  end if
else
   coe=bin(ci)
end if
end subroutine cicjBasis


!------------------------------------------------------
!T operator to one particular basis for particular spin
!------------------------------------------------------
subroutine TBasis(coe,numn,Nfi,Nsp,dd)
use param
implicit none
integer,intent(in)::Nsp,dd
integer,intent(in)::Nfi(Nsp)
integer,intent(out)::coe
integer(kind=8),intent(inout)::numn
integer::bin(Nlattice)
integer::NNfi(Nsp)
integer::i,j,k
integer::totals

do i=1,Nsp,1
   NNfi(i)=Tmatrix(Nfi(i),dd)
end do

bin=0
do i=1,Nsp,1
   bin(NNfi(i))=1
end do

call get_num_sp(bin,Numn,Nsp)

totals=0
do i=Nsp-1,1,-1
   do j=i+1,Nsp,1
      if(NNfi(i).GT.NNfi(j)) then
         totals=totals+1
      end if
   end do
end do
coe=(-1)**totals

end subroutine TBasis
