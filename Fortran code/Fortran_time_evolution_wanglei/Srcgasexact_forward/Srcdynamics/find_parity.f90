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
