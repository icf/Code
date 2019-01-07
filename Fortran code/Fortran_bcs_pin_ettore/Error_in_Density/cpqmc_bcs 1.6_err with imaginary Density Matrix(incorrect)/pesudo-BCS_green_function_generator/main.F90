!main program of code
program main

implicit none
integer::Nsite,Nspin,i,j,total_num,fac_Nspin,fac_Nsite,fac_NN
complex(kind=8),allocatable::pBCS(:),mark(:),set(:),green(:)
complex(kind=8)::value,sum_green

!set lattice
Nsite=16
Nspin=7

call FAC(Nsite-1,fac_Nsite)
call FAC(Nspin-1,fac_Nspin)
call FAC(Nsite-Nspin,fac_NN)

total_num=fac_Nsite/(fac_Nspin*fac_NN)

allocate(mark(Nsite-1))
allocate(set(Nsite-1))
allocate(pBCS(Nsite))
allocate(green(Nsite))

!read pBCS
!open(unit=10,file='pesudo-BCS.inputdat',status='UNKNOWN')
!do i=1,Nsite,1
!        write(10,*)pBCS(i)
!end do
!close(10)
pBCS(1)=      2.0027441349543061E-002
pBCS(2)=    3.8298655004720572E-002
pBCS(3)=    4.5667145993058413E-002
pBCS(4)=    5.4371978169259441E-002
pBCS(5)=    7.0135519665785628E-002
pBCS(6)=   0.19177861125993251     
pBCS(7)=   0.21905259746086969     
pBCS(8)=   0.23423162495675626     
pBCS(9)=   0.30043604145437530     
pBCS(10)=   0.55237766742782823     
pBCS(11)=   0.65035292999084193     
pBCS(12)=   0.88629295354583792     
pBCS(13)=   0.90727229377100671     
pBCS(14)=   0.90727229377100671    
pBCS(15)=   0.90727229377100671     
pBCS(16)=   0.90727229377100671



do i=1,Nsite
   pBCS(i)=pBCS(i)/abs(1-pBCS(i))
   write(*,*)pBCS(i)
   !pBCS(i)=sqrt(pBCS(i))
enddo

green=0
value=0
!calculate green function
do i=1,Nsite

   call initial(pBCS,i,Nsite,Nspin,mark,set)
   !write(*,*)real(set)
   !write(*,*)real(mark);pause
   call get_mark_value(i,set,pBCS,mark,value,Nsite)
   green(i)=green(i)+value
   !write(*,*)real(value);pause
   !write(*,*)total_num
   do j=2,total_num
      call update_mark(mark,Nsite)
      !write(*,*)real(mark);pause
      call get_mark_value(i,set,pBCS,mark,value,Nsite)
      green(i)=green(i)+value
      !write(*,*)real(value);pause
   enddo

enddo

!output
sum_green=sum(green)
green=(green/sum_green)*Nspin

do i=1,Nsite
   write(*,*)real(green(i))
enddo
write(*,*)'sum',real(sum(green))

deallocate(set)
deallocate(mark)
deallocate(pBCS)
deallocate(green)

end program main
!--------------------------

subroutine initial(pBCS,i,Nsite,Nspin,mark,set)

implicit none
integer:: i,j,counter,Nsite,Nspin
complex(kind=8)::pBCS(Nsite)
complex(kind=8)::mark(Nsite-1),set(Nsite-1)

counter=0
do j=1,Nsite
   if(i .NE. j)then
      counter=counter+1
      set(counter)=(pBCS(j))
   endif
enddo

mark=0

do j=1,Nspin-1
   mark(j)=1
enddo

end subroutine initial
!-------------------------

subroutine get_mark_value(i,set,pBCS,mark,value,Nsite)

implicit none
integer::i,j,Nsite
complex(kind=8)::value
complex(kind=8)::pBCS(Nsite)
complex(kind=8)::mark(Nsite-1),set(Nsite-1)

value=pBCS(i)
do j=1,Nsite-1
   if(mark(j) .EQ. 1)value=value*set(j)
enddo

end subroutine get_mark_value
!------------------------

subroutine update_mark(mark,Nsite)

implicit none
integer::j,mark_down,Nsite,left_move
complex(kind=8)::mark(Nsite-1)


do j=1,Nsite-2
   if(mark(j) .EQ. 1 .AND. mark(j+1) .EQ. 0)then
      mark(j)=0
      mark(j+1)=1
      mark_down=j
      goto 10
   endif
enddo

10 left_move=0
do j=1,mark_down-1
   if (mark(j) .EQ. 0)left_move=left_move+1
   if (mark(j) .EQ. 1)goto 20
enddo

20 left_move=left_move
do j=1,mark_down-1
   if (mark(j) .EQ. 1 .AND. left_move .NE. 0)then
      mark(j-left_move)=1
      mark(j)=0
   endif
enddo

end subroutine update_mark
!-----------------------

subroutine FAC(N,f)

implicit none
integer::f,N,i

f=1
do i=1,N
   f=f*i
enddo

end subroutine FAC























