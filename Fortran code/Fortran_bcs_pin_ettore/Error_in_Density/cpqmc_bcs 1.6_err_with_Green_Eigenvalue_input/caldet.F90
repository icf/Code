module caldet_module
implicit none

INTERFACE caldet
 MODULE PROCEDURE caldet_r, caldet_c
END INTERFACE

contains 

subroutine caldet_r(n,A,det)
implicit none
integer,intent(IN):: n
real(kind=8),intent(IN)::A(n,n)
real(kind=8),intent(OUT):: det

integer:: perm,i

       real(kind=8)::b(n,n)
       integer::ipiv(n),info

       b=a
       call dgetrf(n,n,b,n,ipiv,info)    
       if(info.lt.0) then
         write(6,*)'DGETRF error ,info=,n=',info,n
         print *, a
         stop
       else if(info.gt.0) then
         det=0.d0
         return
       endif

       call find_perm(n,ipiv,perm)

       det=1d0
       do i=1,n
         det=det*b(i,i)
       enddo 
         det=det *dble(perm)

end subroutine caldet_r


subroutine caldet_c(n,A,det)
implicit none
integer,intent(IN):: n
complex(kind=8),intent(IN)::A(n,n)
complex(kind=8),intent(OUT):: det

integer:: perm,i

       complex(kind=8)::b(n,n)
       integer::ipiv(n),info

       b=a
       call zgetrf(n,n,b,n,ipiv,info)
       if(info.lt.0) then
          write(6,*)'in caldet, ZGETRF error ,info=',info
          print *, 'n', n 
          print *, 'A', A
          stop
       else if(info.gt.0) then
          det=dcmplx(0.d0,0.d0)
          return
       endif

       call find_perm(n,ipiv,perm)

!       print *, ipiv, parity
!        pause

       det=1d0
       do i=1,n
         det=det*b(i,i)
       enddo 
         det=det *dble(perm)

end subroutine caldet_c

subroutine find_perm(n,ipiv,perm)
implicit none
integer,intent(IN):: n
integer,intent(IN):: ipiv(n)
integer,intent(OUT):: perm

integer:: i,j, num_inv

!find number of inversions
!num_inv=0
!  do i=1,n-1
!    do j=i+1,n
!        if (ipiv(i)>ipiv(j)) then 
!            num_inv=num_inv+1
!        endif
!    enddo 
!  enddo 

! changed 2010.8.29
num_inv=0
do i=1,n
 if (ipiv(i)/=i) then 
 num_inv=num_inv+1
 endif 
enddo 


     if (mod(num_inv,2)==0) then 
         perm=1
     else
         perm=-1  
     endif

end subroutine find_perm

end module caldet_module
