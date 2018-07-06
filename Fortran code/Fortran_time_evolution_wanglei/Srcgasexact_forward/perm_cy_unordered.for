  !  
  !    Mar-11-08 from CY, modified by LW Mar-12-08
  !----------------------------
  !   a simple drive for illustration purpose is provided 
   !-------------------------------------------------------
   !   program main   
   !   implicit none
!	integer, parameter:: n=4
!	integer:: A(n),A0(n),  k
!
!	do k=1, n
!	A0(k)=k
!	enddo
!	
!	A=A0
!	do k=1, product(A0(:))
!	print *, A
!	call perm_cy(k,n, A0, A)
!	enddo
!	  
!
 !     end program

!Unordered generation
!For every number k, with 0 ¡Ü k < n!, the following algorithm generates a unique permutation of the initial sequence sj, j=1¡­n:
!
!Notation
!k / j denotes integer division of k by j, i.e. the integral quotient without any remainder, and
!k mod j is the remainder following integer division of k by j.
!s[n] denotes the nth element of sequence s.
      
      subroutine perm_cy(k, n, s_in, s_out)
	implicit none 
	integer, intent(IN)::k, n, s_in(n)
	integer, intent(OUT):: s_out(n)
	
      integer(kind=4):: factorial
	integer:: j
	factorial=1
      s_out(:)=s_in(:)
	do j=2, n
	 factorial=factorial*(j-1)
	call exchange_int(  s_out( j-  mod((k / factorial),j)),s_out(j))
      enddo
	return
	end subroutine
	
   !-------------------------
       subroutine exchange_int(a, b)  ! exchange integer numbers 
       implicit none

       integer, intent(INOUT):: a, b
        integer:: c

      c=a
      a=b
      b=c

      endsubroutine


