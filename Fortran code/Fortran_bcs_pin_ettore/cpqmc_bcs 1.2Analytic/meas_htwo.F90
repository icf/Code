!----------------------------------------------
!Get the htwo_local(k)=<phR|H^2|phL>/ <phR|phL>
!the htwo_local should be in 1~Dtot
!----------------------------------------------
subroutine get_htwolocal(htwom_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(OUT)::htwom_local(Dtot)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)

complex(kind=8)::g_hopt(2*Nhop)
integer::g_sit(2*Nhop,2)

integer::i,j,k
integer::i1,i2,i3,i4
integer::i5,i6,i7,i8


complex(kind=8)::tmp

call get_ghop(g_hopt,g_sit)


htwom_local=zero

!<K.K>
call get_kklocal(g_hopt,g_sit,htwom_local,Amat)

!<K.V>+<V.K>
call get_kvlocal(g_hopt,g_sit,htwom_local,Amat)

!<V.V>
call get_vvlocal(g_hopt,g_sit,htwom_local,Amat)

end subroutine get_htwolocal




!get the general hopping terms 
subroutine get_ghop(g_hopt,g_sit)
use lattice_param
implicit none
complex(kind=8),intent(OUT)::g_hopt(2*Nhop)
integer,intent(OUT)::g_sit(2*Nhop,2)
integer::i,j
do i=1,Nhop,1
   g_hopt(i)=hopt(i)
   g_sit(i,1)=sit(i,1)
   g_sit(i,2)=sit(i,2)

   g_hopt(i+Nhop)=hopt(i)
   g_sit(i+Nhop,1)=sit(i,1)+Nsite
   g_sit(i+Nhop,2)=sit(i,2)+Nsite
end do
end subroutine get_ghop



!Get the <phR|KK|phL>/ <phR|phL>
subroutine get_kklocal(g_hopt,g_sit,htwom_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::g_hopt(2*Nhop)
integer,intent(IN)::g_sit(2*Nhop,2)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(INOUT)::htwom_local(Dtot)
integer::i,j,k
integer::i1,i2,i3,i4
integer::i5,i6,i7,i8
complex(kind=8)::tmp


if(dtype.EQ.'c') then !put if outside, avoid too many selections

 do k=1,Dtot,1
    do i=1,2*Nhop,1
       do j=1,2*Nhop,1
          i1=g_sit(i,1)
          i2=g_sit(i,2)
          i3=g_sit(j,1)
          i4=g_sit(j,2)
          call two_body(tmp,i1,i2,i3,i4,Amat(1:2*Nsite,1:2*Nsite,k))
          htwom_local(k)=htwom_local(k)+tmp*g_hopt(i)*g_hopt(j)
       end do
    end do
 end do

else if(dtype.EQ.'d') then

 do k=1,Dtot,1

    do i=1,Nhop,1
       do j=1,Nhop,1
          i1=g_sit(i,1)
          i2=g_sit(i,2)
          i3=g_sit(j,1)
          i4=g_sit(j,2)
          call two_body(tmp,i1,i2,i3,i4,Amat(1:2*Nsite,1:2*Nsite,k))
          htwom_local(k)=htwom_local(k)+tmp*g_hopt(i)*g_hopt(j)
       end do
    end do


    do i=1,Nhop,1
       do j=Nhop+1,2*Nhop,1
          i1=g_sit(i,1)
          i2=g_sit(i,2)
          i3=g_sit(j,1)
          i4=g_sit(j,2)
          htwom_local(k)=htwom_local(k)+g_hopt(i)*g_hopt(j)*Amat(i1,i2,k)*Amat(i3,i4,k)
       end do
    end do


    do i=Nhop+1,2*Nhop,1
       do j=1,Nhop,1
          i1=g_sit(i,1)
          i2=g_sit(i,2)
          i3=g_sit(j,1)
          i4=g_sit(j,2)
          htwom_local(k)=htwom_local(k)+g_hopt(i)*g_hopt(j)*Amat(i1,i2,k)*Amat(i3,i4,k)
       end do
    end do


    do i=Nhop+1,2*Nhop,1
       do j=Nhop+1,2*Nhop,1
          i1=g_sit(i,1)
          i2=g_sit(i,2)
          i3=g_sit(j,1)
          i4=g_sit(j,2)
          call two_body(tmp,i1,i2,i3,i4,Amat(1:2*Nsite,1:2*Nsite,k))
          htwom_local(k)=htwom_local(k)+tmp*g_hopt(i)*g_hopt(j)
       end do
    end do


 end do


end if

end subroutine get_kklocal




!Get the <phR|KV+VK|phL>/ <phR|phL>
subroutine get_kvlocal(g_hopt,g_sit,htwom_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::g_hopt(2*Nhop)
integer,intent(IN)::g_sit(2*Nhop,2)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(INOUT)::htwom_local(Dtot)
integer::i,j,k
integer::i1,i2,i3,i4
integer::i5,i6,i7,i8
complex(kind=8)::tmp


if(dtype.EQ.'c') then 

 do k=1,Dtot,1
    !kv
    do i=1,2*Nhop,1
       do j=1,Nsite,1
          i1=g_sit(i,1)
          i2=g_sit(i,2)

          i3=j
          i4=j

          i5=j+Nsite
          i6=j+Nsite

          call three_body(tmp,i1,i2,i3,i4,i5,i6,Amat(1:2*Nsite,1:2*Nsite,k))

          htwom_local(k)=htwom_local(k)+tmp*onsitU*g_hopt(i)
       end do
    end do
    !vk
    do j=1,Nsite,1
       do i=1,2*Nhop,1
          i1=j
          i2=j

          i3=j+Nsite
          i4=j+Nsite

          i5=g_sit(i,1)
          i6=g_sit(i,2)

          call three_body(tmp,i1,i2,i3,i4,i5,i6,Amat(1:2*Nsite,1:2*Nsite,k))

          htwom_local(k)=htwom_local(k)+tmp*onsitU*g_hopt(i)
       end do
    end do

 end do

else if(dtype.EQ.'d') then

 do k=1,Dtot,1

    !kv
    do i=1,Nhop,1
       do j=1,Nsite,1
          i1=g_sit(i,1)
          i2=g_sit(i,2)

          i3=j
          i4=j

          i5=j+Nsite
          i6=j+Nsite

          call two_body(tmp,i1,i2,i3,i4,Amat(1:2*Nsite,1:2*Nsite,k))

          htwom_local(k)=htwom_local(k)+tmp*Amat(i5,i6,k)*onsitU*g_hopt(i)
       end do
    end do

    do i=Nhop+1,2*Nhop,1
       do j=1,Nsite,1
          i1=g_sit(i,1)
          i2=g_sit(i,2)

          i3=j
          i4=j

          i5=j+Nsite
          i6=j+Nsite

          call two_body(tmp,i1,i2,i5,i6,Amat(1:2*Nsite,1:2*Nsite,k))

          htwom_local(k)=htwom_local(k)+tmp*Amat(i3,i4,k)*onsitU*g_hopt(i)
       end do
    end do

    !vk
    do j=1,Nsite,1
       do i=1,Nhop,1
          i1=j
          i2=j

          i3=j+Nsite
          i4=j+Nsite

          i5=g_sit(i,1)
          i6=g_sit(i,2)

          call two_body(tmp,i1,i2,i5,i6,Amat(1:2*Nsite,1:2*Nsite,k))

          htwom_local(k)=htwom_local(k)+tmp*Amat(i3,i4,k)*onsitU*g_hopt(i)
       end do
    end do


    do j=1,Nsite,1
       do i=Nhop+1,2*Nhop,1
          i1=j
          i2=j

          i3=j+Nsite
          i4=j+Nsite

          i5=g_sit(i,1)
          i6=g_sit(i,2)

          call two_body(tmp,i3,i4,i5,i6,Amat(1:2*Nsite,1:2*Nsite,k))

          htwom_local(k)=htwom_local(k)+tmp*Amat(i1,i2,k)*onsitU*g_hopt(i)
       end do
    end do

 end do

end if
end subroutine get_kvlocal



!Get the <phR|VV|phL>/ <phR|phL>
subroutine get_vvlocal(g_hopt,g_sit,htwom_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::g_hopt(2*Nhop)
integer,intent(IN)::g_sit(2*Nhop,2)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(INOUT)::htwom_local(Dtot)
integer::i,j,k
integer::i1,i2,i3,i4
integer::i5,i6,i7,i8
complex(kind=8)::tmp1,tmp2,tmp

if(dtype.EQ.'c') then !try to avoid too many selections

  do k=1,Dtot,1
     do i=1,Nsite,1
        do j=1,Nsite,1
           i1=i
           i2=i

           i3=i+Nsite
           i4=i+Nsite

           i5=j
           i6=j

           i7=j+Nsite
           i8=j+Nsite

           call four_body(tmp,i1,i2,i3,i4,i5,i6,i7,i8,Amat(1:2*Nsite,1:2*Nsite,k))

           htwom_local(k)=htwom_local(k)+tmp*onsitU*onsitU
        end do
     end do
  end do

else if(dtype.EQ.'d') then
  do k=1,Dtot,1
     do i=1,Nsite,1
        do j=1,Nsite,1
           i1=i
           i2=i

           i3=i+Nsite
           i4=i+Nsite

           i5=j
           i6=j

           i7=j+Nsite
           i8=j+Nsite

           call two_body(tmp1,i1,i2,i5,i6,Amat(1:2*Nsite,1:2*Nsite,k))
           call two_body(tmp2,i3,i4,i7,i8,Amat(1:2*Nsite,1:2*Nsite,k))
           tmp=tmp1*tmp2

           htwom_local(k)=htwom_local(k)+tmp*onsitU*onsitU
        end do
     end do
  end do
end if
end subroutine get_vvlocal



!get <c1+ c2 c3+ c4>
subroutine two_body(tmp,i1,i2,i3,i4,Amat)
use lattice_param
implicit none
complex(kind=8),intent(OUT)::tmp
integer,intent(IN)::i1,i2,i3,i4
complex(kind=8),intent(IN)::Amat(1:2*Nsite,1:2*Nsite)


tmp=Amat(i1,i2)*Amat(i3,i4)
tmp=tmp-Amat(i1,i4)*Amat(i3,i2)
if(i2.eq.i3) tmp=tmp+Amat(i1,i4)
end subroutine two_body



!get <c1+ c2 c3+ c4 c5+ c6>
subroutine three_body(tmp,i1,i2,i3,i4,i5,i6,Amat)
use lattice_param
implicit none
complex(kind=8),intent(OUT)::tmp
integer,intent(IN)::i1,i2,i3,i4,i5,i6
complex(kind=8),intent(IN)::Amat(1:2*Nsite,1:2*Nsite)
complex(kind=8)::tmp_two

!<c1+ c2> <c3+ c4 c5+ c6>
call two_body(tmp_two,i3,i4,i5,i6,Amat(1:2*Nsite,1:2*Nsite))
tmp=Amat(i1,i2)*tmp_two
!<c1+ c4> <c2 c3+ c5+ c6>
!=<c1+ c4>(d32<c5+ c6>-<c3+ c2 c5+ c6>)
call two_body(tmp_two,i3,i2,i5,i6,Amat(1:2*Nsite,1:2*Nsite))
tmp=tmp-Amat(i1,i4)*tmp_two
if(i3.eq.i2) tmp=tmp+Amat(i1,i4)*Amat(i5,i6)
!<c1+ c6> <c2 c3+ c4 c5+>
call two_body(tmp_two,i3,i2,i5,i4,Amat(1:2*Nsite,1:2*Nsite))
tmp=tmp+Amat(i1,i6)*tmp_two
if(i2.eq.i3) tmp=tmp-Amat(i1,i6)*Amat(i5,i4)
if(i5.eq.i4) tmp=tmp-Amat(i1,i6)*Amat(i3,i2)
if(i2.eq.i3.AND.i5.eq.i4) tmp=tmp+Amat(i1,i6)
end subroutine three_body


!get <c1+ c2 c3+ c4 c5+ c6 c7+ c8>
subroutine four_body(tmp,i1,i2,i3,i4,i5,i6,i7,i8,Amat)
use lattice_param
implicit none
complex(kind=8),intent(OUT)::tmp
integer,intent(IN)::i1,i2,i3,i4,i5,i6,i7,i8
complex(kind=8),intent(IN)::Amat(1:2*Nsite,1:2*Nsite)
complex(kind=8)::tmp_three
complex(kind=8)::tmp_two

!<c1+ c2> <c3+ c4 c5+ c6 c7+ c8>
call three_body(tmp_three,i3,i4,i5,i6,i7,i8,Amat(1:2*Nsite,1:2*Nsite))
tmp=Amat(i1,i2)*tmp_three
!<c1+ c4> <c2 c3+ c5+ c6 c7+ c8>
call three_body(tmp_three,i3,i2,i5,i6,i7,i8,Amat(1:2*Nsite,1:2*Nsite))
tmp=tmp-Amat(i1,i4)*tmp_three
if(i2.eq.i3) then
  call two_body(tmp_two,i5,i6,i7,i8,Amat(1:2*Nsite,1:2*Nsite))
  tmp=tmp+Amat(i1,i4)*tmp_two
end if
!<c1+ c6> <c2 c3+ c4 c5+ c7+ c8>
call three_body(tmp_three,i3,i2,i5,i4,i7,i8,Amat(1:2*Nsite,1:2*Nsite))
tmp=tmp+Amat(i1,i6)*tmp_three
if(i5.eq.i4) then
  call two_body(tmp_two,i3,i2,i7,i8,Amat(1:2*Nsite,1:2*Nsite))
  tmp=tmp-Amat(i1,i6)*tmp_two
end if
if(i3.eq.i2) then
  call two_body(tmp_two,i5,i4,i7,i8,Amat(1:2*Nsite,1:2*Nsite))
  tmp=tmp-Amat(i1,i6)*tmp_two
end if
if(i3.eq.i2.AND.i5.eq.i4) then
  tmp=tmp+Amat(i1,i6)*Amat(i7,i8)
end if
!<c1+ c8><c2 c3+ c4 c5+ c6 c7+>
if(i2.eq.i3.AND.i4.eq.i5.AND.i6.eq.i7) then
  tmp=tmp+Amat(i1,i8)
end if
if(i2.eq.i3.AND.i4.eq.i5) then
  tmp=tmp-Amat(i1,i8)*Amat(i7,i6)
end if
if(i2.eq.i3.AND.i6.eq.i7) then
  tmp=tmp-Amat(i1,i8)*Amat(i5,i4)
end if
if(i2.eq.i3) then
  call two_body(tmp_two,i5,i4,i7,i6,Amat(1:2*Nsite,1:2*Nsite))
  tmp=tmp+Amat(i1,i8)*tmp_two
end if
if(i6.eq.i7.AND.i4.eq.i5) then
  tmp=tmp-Amat(i1,i8)*Amat(i3,i2)
end if
if(i4.eq.i5) then
  call two_body(tmp_two,i3,i2,i7,i6,Amat(1:2*Nsite,1:2*Nsite))
  tmp=tmp+Amat(i1,i8)*tmp_two
end if
if(i6.eq.i7) then
  call two_body(tmp_two,i3,i2,i5,i4,Amat(1:2*Nsite,1:2*Nsite))
  tmp=tmp+Amat(i1,i8)*tmp_two
end if
call three_body(tmp_three,i3,i2,i5,i4,i7,i6,Amat(1:2*Nsite,1:2*Nsite))
tmp=tmp-Amat(i1,i8)*tmp_three
end subroutine four_body
