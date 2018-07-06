subroutine zqr(m,n,A,Q,R)
implicit none

integer, intent(IN):: m,n
complex(kind=8),intent(IN):: A(m,n)
complex(kind=8),intent(OUT):: Q(m,m), R(m,n)

integer:: lwork,info
complex(kind=8):: tau(min(n,n))

complex(kind=8),allocatable:: work(:)

complex(kind=8):: H(m,m), v(m)
complex(kind=8):: Temp_mn(m,n), Temp_mm(m,m)
integer::i,j,k


lwork=max(1,n)*10
allocate(work(max(1,lwork)))

!do i=1,m
!print *, i, A(i,:)
!enddo 


call ZGEQRF( M, N, A, M, TAU, WORK, LWORK, INFO )

!get R
R=0d0
do i=1,m
  do j=i,n
R(i,j)=A(i,j)
  enddo 
enddo 


!initialize Q
Q=0d0
do i=1,m
Q(i,i)=1d0
enddo 

do i=1,min(m,n)

v=0d0
v(i)=1d0
v(i+1:m)=A(i+1:m,i)

do j=1,m
  do k=1,m
      if (j==k) then 
        H(j,k)=1d0-tau(i)*v(j)*conjg(v(k))
      else
        H(j,k)=-tau(i)*v(j)*conjg(v(k))
      endif
   enddo
enddo

Q=matmul(Q,H)
enddo 

deallocate(work)

if (info/=0) then 
print *,'error happens when zqr'
stop
endif

end subroutine 
