!main program of code
program main

implicit none
integer:: Nsite,Nsitex,Nsitey
complex(kind=8),allocatable::cicj_sc_global(:,:),cicj_sc_global_error(:,:)

integer::sitex,sitey
integer::i,j

Nsitex=4
Nsitey=4
Nsite=16

allocate(cicj_sc_global(2*Nsite,2*Nsite))
allocate(cicj_sc_global_error(2*Nsite,2*Nsite))

!input
open(unit=10,file='cicj_sc_global_input.inputdat',status='OLD')
do i=1,2*Nsite,1
  do j=1,2*Nsite,1
        read(10,*)cicj_sc_global(i,j)
  end do
end do
close(10)

open(unit=20,file='cicj_sc_global_error_input.inputdat',status='OLD')
do i=1,2*Nsite,1
  do j=1,2*Nsite,1
        read(20,*)cicj_sc_global_error(i,j)
  end do
end do
close(20)





!output
open(unit=30,file='cicj_sc_global_input.inputdat',status='UNKNOWN')
do i=1,2*Nsite,1
  do j=1,2*Nsite,1

  end do
end do
close(30)

deallocate(cicj_sc_global)
deallocate(cicj_sc_global_error)
 
end program main
