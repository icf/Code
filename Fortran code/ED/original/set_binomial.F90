!This page contains the subroutine need to use in the new label routine.


!1.set_binomial:Set the binomial table we need to use in our label method and we
!             can use biomial(0:Nlattice,0:Nlattice) to the number of hilbert space
!2.next:we use it update the configuration, contains Nfi and dltNfi==>one spin
!3.printN: when one configuration nfi is given,it output the number of nfi
!4.get_num:give out the number base on the configuration bin(2*Nlattice)==>two
!          spin conditions
!5.get_nfi:When the number is given, it output the nfi==>one spin


!---------------------------------------------------------
!Set the binomial table we need to use in our label method.
!---------------------------------------------------------
subroutine set_binomial
use param
implicit none
integer::i,j,k
integer(kind=8)::binom(0:Nlattice),binomtemp(0:Nlattice)

Binomial=0 !Set Binomial
binom=0;binomtemp=0
binomtemp(0)=1;binomtemp(1)=1
Binomial(0,1)=1;Binomial(1,1)=1!Set Binomial
do i=2,Nlattice,1
   binom(0)=1
   do j=1,i-1,1
      binom(j)=binomtemp(j-1)+binomtemp(j)
   end do
   binom(i)=1
   Binomial(0:i,i)=binom(0:i) !Set Binomial
   binomtemp=binom
end do
end subroutine set_binomial


!-----------------------------------------------
!update the configuration
!if flag=1 means get the new configuration
!if flag=0 means there is no new configuration
!-----------------------------------------------
subroutine next(Nfi,dltNfi,flag,Ns)
implicit none
integer,intent(out)::flag
integer,intent(in)::Ns
integer,intent(inout)::Nfi(Ns),dltNfi(Ns)
integer::i,j
flag=0
do i=1,Ns,1
   if(dltNfi(i).GT.0) then
     Nfi(i)=Nfi(i)+1
     dltNfi(i)=dltNfi(i)-1
     do j=1,i-1,1
        Nfi(j)=j
        dltNfi(j)=0
     end do
     if(i.GT.1) dltNfi(i-1)=Nfi(i)-Nfi(i-1)-1
     flag=1
     exit
   end if
end do

end subroutine next


!--------------------------------------------------------------
!give out the number base on the configuration,based on one Nfi
!--------------------------------------------------------------
subroutine printN(Nfi,Ns)
use param
implicit none
!integer,intent(IN)::vec(Nlattice)
integer,intent(in)::Ns
integer,intent(in)::Nfi(Ns)
integer(kind=8)::Nklabel
integer::i,j
!j=0
!do i=1,Nlattice,1
!   if(vec(i).EQ.1) then
!     j=j+1
!     Nfi(j)=i     
!   end if
!end do
!if(j.NE.Ns) write(*,*) "Wrong in printN"

Nklabel=0
do i=Ns,1,-1
   Nklabel=Nklabel+Binomial(i,Nfi(i)-1)
end do
Nklabel=Nklabel+1

write(*,*) Nklabel
end subroutine printN


!-----------------------------------------------------------
!give out the number base on the configuration bin(Nlattice)
!-----------------------------------------------------------
subroutine get_num_sp(vec,Num,Nsp)
use param
implicit none
integer,intent(IN)::vec(Nlattice)
integer,intent(IN)::Nsp
integer(kind=8),intent(OUT)::Num
integer::Nfi(Nsp)
integer::i,j
j=0
do i=1,Nlattice,1
   if(vec(i).EQ.1) then
     j=j+1
     Nfi(j)=i
   end if
end do
if(j.NE.Nsp) then
  write(*,*) "Wrong in get_num_sp"
  call mystop
end if
Num=0
do i=Nsp,1,-1
   Num=Num+Binomial(i,Nfi(i)-1)
end do
Num=Num+1
end subroutine get_num_sp


!-------------------------------------------------------------
!give out the number base on the configuration bin(2*Nlattice)
!-------------------------------------------------------------
subroutine get_num(vec,Num)
use param
implicit none
integer,intent(IN)::vec(2*Nlattice)
integer(kind=8),intent(OUT)::Num
integer::gNfi_up(Nup),gNfi_dn(Ndn)
integer(kind=8)::Nkup,Nkdn
integer::i,j
j=0
do i=1,Nlattice,1
   if(vec(2*i-1).EQ.1) then
     j=j+1
     gNfi_up(j)=i 
   end if
end do
if(j.NE.Nup) then
  write(*,*) "Wrong in get_num up"
  call mystop
end if

j=0
do i=1,Nlattice,1
   if(vec(2*i).EQ.1) then
     j=j+1
     gNfi_dn(j)=i 
   end if
end do
if(j.NE.Ndn) then
  write(*,*) "Wrong in get_num down"
  call mystop
end if

Nkup=0
do i=Nup,1,-1
   Nkup=Nkup+Binomial(i,gNfi_up(i)-1)
end do
Nkup=Nkup+1

Nkdn=0
do i=Ndn,1,-1
   Nkdn=Nkdn+Binomial(i,gNfi_dn(i)-1)
end do
Nkdn=Nkdn+1

Num=Nkup+Binomial(Nup,Nlattice)*(Nkdn-1)
end subroutine get_num




!------------------------------------------------------------------------
!when you have the num,this subroutine will give out nfi and dltnfi 
!so that you can use at the first of next subroutine.(num is single spin)
!------------------------------------------------------------------------
subroutine get_nfi(num,nfi,dltnfi,ns)
use param
implicit none
integer(kind=8),intent(in)::num
integer,intent(in)::ns
integer,intent(out)::nfi(ns),dltnfi(ns)
integer(kind=8)::ntmp
integer::i,j

integer::ni
ntmp=num-1

!if(ntmp.eq.0) return

if(ns.eq.0) return

ni=0
do i=Nlattice,1,-1
   if(Binomial(Ns,i-1).LE.ntmp) then
     nfi(ns)=i
     dltnfi(ns)=Nlattice-i
     ntmp=ntmp-Binomial(Ns,i-1)
     ni=1
     exit
   end if
end do

!if(ntmp.eq.0) return

do j=ns-1,1,-1
   do i=nfi(j+1)-1,1,-1
      if(Binomial(j,i-1).LE.ntmp) then
         nfi(j)=i
         dltnfi(j)=nfi(j+1)-nfi(j)-1
         ntmp=ntmp-Binomial(j,i-1)
         ni=ni+1
         exit
      end if
   end do
end do

if(ntmp.ne.0) then
  write(*,*) "Something is wrong in get_nfi ntmp"
  call mystop
end if

if(ni.NE.ns) then
  write(*,*) "Something is wrong in get_nfi ni"
  call mystop
end if

end subroutine get_nfi
