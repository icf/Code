subroutine caldet_c(n,A,det)
implicit none
integer,intent(IN):: n
complex(kind=8),intent(IN)::A(n,n)
complex(kind=8),intent(OUT):: det

integer:: parity,i

       complex(kind=8)::b(n,n)
       integer::ipiv(n),info

       b=a
       call zgetrf(n,n,b,n,ipiv,info)
       if(info/=0) then
        write(6,*)'ZGETRF error ,info=',info
        print *, 'incaldet a', a, n
        stop
       endif
 
       call find_parity(n,ipiv,parity)

       det=1d0
       do i=1,n
         det=det*b(i,i)
       enddo 
         det=det *dble(parity)

end subroutine caldet_c

subroutine find_parity(n,ipiv,parity)
implicit none
integer,intent(IN):: n
integer,intent(IN):: ipiv(n)
integer,intent(OUT):: parity

integer:: i,j, num_inv

!find number of inversions
num_inv=0
  do i=1,n-1
    do j=i+1,n
        if (ipiv(i)>ipiv(j)) then 
            num_inv=num_inv+1
        endif
    enddo 
  enddo 

     if (mod(num_inv,2)==0) then 
         parity=1
     else
         parity=-1  
     endif

end subroutine find_parity


subroutine cal_overlap(n,m,phi_l,phi_r,overlap)
implicit none

integer,intent(IN):: n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: overlap

complex(kind=8):: tempmat(m,m)

if (isnan(dble(phi_l(1,1)))) then 
print *, 'wrong mat to caldet'
stop
endif 

tempmat= matmul(transpose(conjg(phi_l)), phi_r)
call caldet_c(m,tempmat,overlap)

end subroutine 


