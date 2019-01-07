!main program of code
program main

implicit none
integer::Nsite,Nspin,i,j,total_num,fac_Nspin,fac_Nsite,fac_NN
complex(kind=8),allocatable::pBCS(:),mark(:),set(:),green(:)
complex(kind=8)::value,sum_green

!set lattice
Nsite=16
Nspin=5

call FNC(total_num,Nsite,Nspin)

write(*,*)"watch out!!! fac_Nsite may limited by int_size --> results for 4*I4"

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

!4477u=8 t0.01t0.02
pBCS(1)=       8.2944186404347406E-002
pBCS(2)=    0.13340081421483782 
pBCS(3)=    0.13344510870228987     
pBCS(4)=   0.13406325275718692 
pBCS(5)=    0.13425567229569418 
pBCS(6)=    0.16270273365080343    
pBCS(7)=  0.16380146332085119      
pBCS(8)=   0.44535609427839512    
pBCS(9)=   0.48102344479411840   
pBCS(10)=  0.51840202230960142    
pBCS(11)=  0.55101049225777377    
pBCS(12)=   0.79222764300357107     
pBCS(13)=    0.79428524007380519    
pBCS(14)=   0.79738870407163664    
pBCS(15)=   0.79839155969966680     
pBCS(16)=   0.87730109132826306

!4433u=8
pBCS(1)=   2.3042412541613433E-003
pBCS(2)=   4.9857982158276783E-003
pBCS(3)=   4.9857982165344584E-003
pBCS(4)=   4.9857982169589461E-003
pBCS(5)=   4.9857982176932346E-003
pBCS(6)=   1.0081847122723760E-002
pBCS(7)=   1.0081847125658964E-002
pBCS(8)=   1.0081847128341504E-002
pBCS(9)=   1.0081847129364723E-002
pBCS(10)=   1.0081847132195921E-002
pBCS(11)=   1.0081847134661363E-002
pBCS(12)=  0.49555768629773750     
pBCS(13)=  0.49555768676740131     
pBCS(14)= 0.49555768737821226     
pBCS(15)=  0.49555768784792076     
pBCS(16)=  0.93503073481453436 

!4433u=8 VMpBCS
pBCS(1)=   0
pBCS(2)=   0.00723835
pBCS(3)=   0.00723835
pBCS(4)=   0.00723835
pBCS(5)=   0.00723835
pBCS(6)=   0.0141306
pBCS(7)=   0.0141306
pBCS(8)=   0.0141306 
pBCS(9)=   0.0141306
pBCS(10)=  0.0141306
pBCS(11)=  0.0141306
pBCS(12)=  0.490692     
pBCS(13)=  0.490692     
pBCS(14)=  0.490692      
pBCS(15)=  0.490692    
pBCS(16)=  0.910919

!4433u=8 t0.01t0.02
pBCS(1)=   2.3047360190699140E-003;
pBCS(2)=   4.9687315076646678E-003;
pBCS(3)=   4.9798367561196847E-003;
pBCS(4)=   4.9973637749746768E-003;
pBCS(5)=   5.0037653317879125E-003;
pBCS(6)=   9.9778572858311744E-003;
pBCS(7)=   1.0048959477296589E-002;
pBCS(8)=   1.0082967971097020E-002;
pBCS(9)=   1.0104375649555835E-002;
pBCS(10)=   1.0120340369560648E-002;
pBCS(11)=  1.0191952546723332E-002;
pBCS(12)=  0.49530287841117793;     
pBCS(13)=  0.49542307337813962;     
pBCS(14)=  0.49567953412719967;     
pBCS(15)=  0.49581733122279181;     
pBCS(16)=  0.93499629617096369;

!4433u=8 t0.01t0.02 VMpBCS
pBCS(1)=   0
pBCS(2)=   0.00700254
pBCS(3)=   0.00701376
pBCS(4)=   0.00703146
pBCS(5)=   0.00703793
pBCS(6)=   0.0139749
pBCS(7)=   0.0140464
pBCS(8)=   0.0140806
pBCS(9)=   0.0141113
pBCS(10)=  0.0141553
pBCS(11)=  0.0142649
pBCS(12)=  0.490184     
pBCS(13)=  0.490303     
pBCS(14)=  0.490163     
pBCS(15)=  0.490176     
pBCS(16)=  0.91022

!4433u=8 t0.01t0.02 pinning 0.1
pBCS(1)=   2.3015934782632526E-003
pBCS(2)=   4.9593439516630016E-003
pBCS(3)=   4.9716478998392205E-003
pBCS(4)=   4.9892948426129508E-003
pBCS(5)=   4.9945004930035813E-003
pBCS(6)=    9.9199071856902361E-003
pBCS(7)=   9.9600103979486532E-003
pBCS(8)=   9.9739741341896507E-003
pBCS(9)=   1.0174015131442041E-002
pBCS(10)=   1.0178063723799209E-002
pBCS(11)=  1.0214687091310399E-002
pBCS(12)=  0.49530226356172502    
pBCS(13)=  0.49543869106856586      
pBCS(14)=  0.49569431240271261      
pBCS(15)=  0.49581553576199994     
pBCS(16)=  0.93511215887531818  

!4455u=8
pBCS(1)=   3.0158113924851477E-002
pBCS(2)=   3.0299109044342496E-002
pBCS(3)=   3.0542983660082741E-002
pBCS(4)=   3.0645025397496046E-002
pBCS(5)=   3.2516454282798000E-002
pBCS(6)=   7.8899421954550186E-002
pBCS(7)=   8.0644679069739073E-002
pBCS(8)=   8.1408948938182291E-002
pBCS(9)=   8.1591944889586099E-002
pBCS(10)=   8.2418332264636601E-002
pBCS(11)=   8.4204777156725127E-002
pBCS(12)=  0.85775916374131667     
pBCS(13)=  0.85888220980855245     
pBCS(14)=  0.86102579874921659     
pBCS(15)=  0.86206228442197930     
pBCS(16)=  0.91694075270317810   

!4455u=8 VMpBCS
pBCS(1)=0.0150479
pBCS(2)=0.0456552
pBCS(3)=0.0458477
pBCS(4)=0.0459439
pBCS(5)=0.0477725
pBCS(6)=0.108148
pBCS(7)=0.109822
pBCS(8)=0.110609
pBCS(9)=0.110806
pBCS(10)=0.111594
pBCS(11)=0.116523
pBCS(12)=0.809404
pBCS(13)=0.810511
pBCS(14)=0.812624
pBCS(15)=0.813645
pBCS(16)=0.886046



  
  


do i=1,Nsite
   pBCS(i)=pBCS(i)/abs(1-pBCS(i))
   pBCS(i)=sqrt(pBCS(i))
   write(*,*)pBCS(i)
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

value=pBCS(i)**2
do j=1,Nsite-1
   if(mark(j) .EQ. 1)value=value*set(j)**2
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

subroutine FNC(f, Nsite, Nspin)

implicit none
integer:: i,f,Nsite,Nspin
real(kind=8):: fr

   fr=1
   do i=1,Nsite-1
      fr=fr*real(i)
      if(i .LE. Nspin-1)fr=fr/real(i);
      if(i .LE. Nsite-Nspin)fr=fr/real(i);
   enddo

f=int(fr)

end subroutine FNC






















