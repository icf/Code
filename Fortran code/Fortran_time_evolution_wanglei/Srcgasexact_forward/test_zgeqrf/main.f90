program test_zgeqrf
implicit none


integer,parameter:: m=2
integer,parameter:: n=2

integer,parameter:: lwork=max(1,n)*10
integer:: info
complex(kind=8):: A(m,n), tau(min(n,n)), work(max(1,lwork))

complex(kind=8):: Q(m,m), R(m,n),H(m,m), v(m)
complex(kind=8):: Temp_mn(m,n), Temp_mm(m,m)
integer::i,j,k

!----- initializa A--------
A(:,1)=(/3,1/)
A(:,2)=(/2,2/)

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

!-------
print *, 'test QdagQ'
Temp_mm=matmul(conjg(transpose(Q)),Q)
do i=1,m
print *, i, Temp_mm(i,:)
enddo 


!----
print *, 'test QR'
Temp_mn=matmul(Q,R)
do i=1,m
print *, i, Temp_mn(i,:)
enddo 


end program 
