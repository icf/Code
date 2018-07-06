program test_zgeqrf
implicit none


integer,parameter:: m=2
integer,parameter:: n=1

complex(kind=8):: A(m,n)
complex(kind=8):: Q(m,m), R(m,n)
complex(kind=8):: Temp_mn(m,n), Temp_mm(m,m)
integer::i,j,k

!----- initializa A--------
A(:,1)=(/3,1/)
!A(:,2)=(/2,2/)

call zqr(m,n,A,Q,R)


!-------
print *, 'Q'
do i=1,m
print *, i, Q(i,:)
enddo 


!-------
print *, 'R'
do i=1,m
print *, i, R(i,:)
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
