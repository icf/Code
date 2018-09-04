!This page contains the hopping information of the hubbard model.

!1.set_lattice_hop:set hop information to sit and hopt arrays depend on set conditions.
!2.setnumber:label the lattice by 1~Nlattice,set by (x,y,z) into coor(Nlattice,Dimen)
!3.sethopPBC:set the sit and hopt by PBC,write into 'hop'
!4.sethopTBC:set the sit and hopt by TBC,write into 'hop'


!---------------------------------------------------------------------
!This subroutine set hop information to sit and hopt arrays depend on
!set conditions: by hand,PBC or TBC
!---------------------------------------------------------------------
subroutine set_lattice_hop()
use param
implicit none
integer::mi,mj,mk
if(set.EQ.1) then
  open(unit=10,file='hop',status='old')
    do mi=1,Nhop,1
       read(10,*) sit(mi,1),sit(mi,2),hopt(mi)
    end do
  close(10)
else
  call setnumber()
!  do mi=1,Nlattice,1
!     write(*,*) "mi=",mi
!     do mj=1,Dimen,1
!        write(*,*) coor(mi,mj)
!     end do
!     write(*,*) "---------------"
!     pause
!  end do
   if(set.EQ.2) then
     call sethopPBC()
   else if(set.EQ.3) then
     call sethopTBC()
   else
     write(*,*) "something is wrong with the set"
     call mystop
   end if
   open(unit=10,file='hop',status='old')
     do mi=1,Nhop,1
        read(10,*) sit(mi,1),sit(mi,2),hopt(mi)
     end do
   close(10)
end if
end subroutine set_lattice_hop


!-----------------------------------------
!we label the lattice by 1~Nlattice,set
!by (x,y,z) the coor(Nlattice,Dimen) is
!dependent on the dimension of the lattice
!-----------------------------------------
subroutine setnumber
use param
implicit none
integer::i,j,k
integer::ntemp,den
 do i=1,Nlattice,1
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
    !write(*,*) coor(i,1),coor(i,2)
 end do
 !call mystop
end subroutine setnumber


!--------------------------------------------------------------------
!we get the hoping matrix by period boundary condition,also test Nhop
!--------------------------------------------------------------------
subroutine sethopPBC
use param
implicit none
integer::i,j,k,ntemp,den
Nhop=0
open(unit=10,file="hop",status='replace')
do i=1,Nlattice,1
   do j=1,Dimen,1

      if(coor(i,j).EQ.1) then

        den=1
        do k=1,j-1,1
           den=den*Nl(k)
        end do
        ntemp=(Nl(j)-1)*den+i
        write(10,*) i,ntemp,kin
        Nhop=Nhop+1
        !Change the code
        ntemp=den+i
        write(10,*) i,ntemp,kin
        Nhop=Nhop+1
        !if(Nl(j).GT.2) then
        !  ntemp=den+i
        !  write(10,*) i,ntemp,kin
        !  Nhop=Nhop+1
        !end if
        
      else if(coor(i,j).EQ.Nl(j)) then
        den=1
        do k=1,j-1,1
           den=den*Nl(k)
        end do
        ntemp=(1-Nl(j))*den+i
        write(10,*) i,ntemp,kin   
        Nhop=Nhop+1
        !Change the code
        ntemp=i-den
        write(10,*) i,ntemp,kin
        Nhop=Nhop+1
        !if(Nl(j).GT.2) then
        !  ntemp=i-den
        !  write(10,*) i,ntemp,kin
        !  Nhop=Nhop+1
        !end if

      else
        den=1
        do k=1,j-1,1
           den=den*Nl(k)
        end do
        ntemp=i+den
        write(10,*) i,ntemp,kin
        Nhop=Nhop+1
        ntemp=i-den
        write(10,*) i,ntemp,kin
        Nhop=Nhop+1
        

      end if

   end do
end do
close(10)

if(Nhop.ne.Nlattice*2*Dimen) then
  write(*,*) "Something is wrong with Nhop PBC"
  call mystop
end if
end subroutine sethopPBC



!-------------------------------------------------------------------
!we get the hoping matrix by Twist boundary condition,also test Nhop
!-------------------------------------------------------------------
subroutine sethopTBC
use param
implicit none
integer::i,j,k,ntemp,den

Nhop=0
open(unit=10,file="hop",status='replace')
do i=1,Nlattice,1
   do j=1,Dimen,1

      if(coor(i,j).EQ.1) then

       ! if(Nl(j).GT.2) then
          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do
          ntemp=(Nl(j)-1)*den+i
          write(10,*) i,ntemp,kin*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
          Nhop=Nhop+1

          ntemp=den+i
          write(10,*) i,ntemp,kin*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
          Nhop=Nhop+1
       ! else

       !   den=1
       !   do k=1,j-1,1
       !      den=den*Nl(k)
       !   end do
       !   ntemp=(Nl(j)-1)*den+i
       !   write(10,*) i,ntemp,kin!*exp(-1.d0*(0.d0,1.d0)*kbound(j)*Pi) !we can decide whether
       !                          ! to use TBC when Nl(j)=2,it does not need to add twist
       !   Nhop=Nhop+1
       ! end if
        
      else if(coor(i,j).EQ.Nl(j)) then
       ! if(Nl(j).GT.2) then
          den=1
          do k=1,j-1,1
             den=den*Nl(k)
          end do
          ntemp=(1-Nl(j))*den+i
          write(10,*) i,ntemp,kin*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
          Nhop=Nhop+1

          ntemp=i-den
          write(10,*) i,ntemp,kin*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
          Nhop=Nhop+1
       ! else
       !
       !   den=1
       !   do k=1,j-1,1
       !      den=den*Nl(k)
       !   end do
       !   ntemp=(1-Nl(j))*den+i
       !   write(10,*) i,ntemp,kin!*exp(1.d0*kbound(j)*Pi) 
       !   Nhop=Nhop+1
       ! end if

      else
        den=1
        do k=1,j-1,1
           den=den*Nl(k)
        end do
        ntemp=i+den
        write(10,*) i,ntemp,kin*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
        Nhop=Nhop+1
        ntemp=i-den
        write(10,*) i,ntemp,kin*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
        Nhop=Nhop+1
        
      end if

   end do
end do
close(10)

!write(*,*) Nhop;call mystop
if(Nhop.ne.Nlattice*2*Dimen) then
  write(*,*) "Something is wrong with Nhop TBC"
  call mystop
end if
end subroutine sethopTBC

!-------------------------------------------
!get the shift matrix in different direction
!-------------------------------------------
subroutine set_Tmatrix
use param
implicit none
integer,external::bound
integer::i,j,k
integer::NT(Dimen)
integer::tmp,dem
do i=1,Nlattice,1
   do j=1,Dimen,1

      tmp=1;dem=1
      do k=1,Dimen,1
         NT(k)=coor(i,k)
         if(k.eq.j) then
           NT(k)=bound(NT(k)-1,Nl(k))  !Tmatrix is inverse Translate here.
         end if

         tmp=tmp+(NT(k)-1)*dem
         dem=dem*Nl(k)
      end do

      Tmatrix(i,j)=tmp

      if(tmp.GT.Nlattice.OR.tmp.LT.1) then
        write(*,*) "Something is wrong with the Tmatrix",tmp
        call mystop
      end if
      !if(rank.eq.0) write(*,*) i,j,Tmatrix(i,j)
      !pause
   end do
end do
      !call mystop
end subroutine set_Tmatrix

!--------------------------------------------------------
!It is a bound function use to fit the boundary condtion.
!--------------------------------------------------------
integer function bound(i,ni)
integer::i,ni
bound=i
if(i.GT.ni) bound=i-ni
if(i.LT.1) bound=i+ni
if(bound.GT.ni.OR.bound.LT.1) write(*,*) "Something is wrong is the boundary condition"
end function

