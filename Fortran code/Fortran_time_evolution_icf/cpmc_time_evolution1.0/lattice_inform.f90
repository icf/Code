!This subroutine get the Hamilton of the no interaction lattice
subroutine set_lattice()
use lattice_param
use model_param
use mpi_serial_param
implicit none
integer::mi,mj

 !Set the number of lattice and Nhop
 if(set.NE.1) then
   Nsite=1
   do mi=1,Dimen,1
    Nsite=Nsite*Nl(mi)
   end do
   Nhop=Nsite*2*Dimen
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
 call check_Hermite_c(Hzero,2*Nsite)

 call set_Tmatrix()

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
  call sethopTBC()
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
implicit none
integer::i,j,k
integer::ntemp,den
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
end subroutine setnumber



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
do i=1,Nsite,1
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
