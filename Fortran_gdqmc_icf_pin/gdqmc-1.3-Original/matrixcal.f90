!This page contains same general matrix operator subroutines. Get from Lei Wang
!1.eigen:Diagonalize hermitian Matrix
!2.general_eigen:Diagonalize of general matrix
!3.inverse:inverse a complex matrix.
!4.zdirectmul:direct multiply of two complex matrix
!5.reigen:Diagonalize real matrix
!6.rinverse:inverse a real matrix
!7.check_Hermite_c:check Hermite of a complex matrix


!This code is changed from Lei Wang's code, I add something more.
!Date: May-22-2013
!Author: Hao Shi
!Email: boruoshihao@gmail.com


!Diagonalization of One Hermitian Matrix, get the eigenvalue and eigenvectors
subroutine eigen(a,n,w)
implicit none
integer,intent(IN)::n
complex(kind=8),intent(INOUT)::a(n,n)
real(kind=8),intent(OUT)::w(n)
complex(kind=8),allocatable::work(:)
real(kind=8),allocatable::rwork(:)
complex(kind=8)::work_test(1)
integer::lwork,lrwork
integer::info

!allocate rwork
lrwork=max(1,3*n-2)
allocate(rwork(lrwork))


!get the best lwork
call zheev('V','U',n,a,n,w,work_test,-1,rwork,info)
lwork=work_test(1)!;lwork=10*n
allocate(work(lwork))


!digonalize the matrix
call zheev('V','U',n,a,n,w,work,lwork,rwork,info)


!deallocate the matrix
deallocate(work,rwork)


if(info/=0) then
   write(6,*)'ZHEEV error ,info=',info
   stop
endif

!call mystop
end subroutine eigen



!Lei's old code, need to be changed some day.  
subroutine general_eigen(a,n,w)
implicit none
integer::n,info
complex*16::a(n,n)
complex*16::w(n)
complex*16::work(10*n)
real*8::rwork(2*n)
complex*16:: vl(n,n),vr(n,n)
call zgeev('N','V',n,a,n,w,vl,n,vr,n,work,10*n,rwork,info)
if(info/=0) then
  write(6,*)'ZgEEV error ,info=',info
  stop
endif
  a=vr
end subroutine general_eigen

 
 

!Evaluate the Inverse of one complex Matrix
subroutine inverse(a,n)
implicit none
integer,intent(IN)::n
complex(kind=8),intent(INOUT)::a(n,n)
complex(kind=8),allocatable::work(:)
complex(kind=8)::work_test(1)
integer::lwork
integer::ipiv(n),info


call zgetrf(n,n,a,n,ipiv,info)
if (info < 0) then
  write(*,"(' zgetrf error info=',I4)") info
  write(*,"(' The',I4,'-th parameter is illegal.')") info
  stop
else if (info > 0) then
  write(*,"(' zgetrf error info=',I4)") info
  write(*,"(' Factorization is done. But the ',I4,'-th diagonal element is zero.')") info
  stop
! we can calculate determinant if everything is fine
end if


!get the best lwork
call zgetri(n,a,n,ipiv,work_test,-1,info)
lwork=work_test(1)!;lwork=10*n
allocate(work(lwork))


!get the inverse
call zgetri(n,a,n,ipiv,work,lwork,info)


!deallocate arrays
deallocate(work)


!Check the point
if (info < 0) then
  write(*,"(' zgetri error info=',I4)") info
  write(*,"(' The',I4,'-th parameter is illegal.')") info
  stop
else if (info > 0) then
  write(*,"(' zgetri error info=',I4)") info
  write(*,"(' The ',I4,'-th diagonal element is zero. Inversion cannot be performed.')") info
  stop
end if

end subroutine inverse



!Get the inverse of a 2*2 matrix
subroutine inv2(inv)
implicit none
complex(kind=8),intent(INOUT)::inv(2,2)
complex(kind=8)::a(2,2),b
a(1:2,1:2)=inv(1:2,1:2)
b=a(1,1)*a(2,2)-a(1,2)*a(2,1)
inv(1,1)=a(2,2)/b
inv(1,2)=-a(1,2)/b
inv(2,1)=-a(2,1)/b
inv(2,2)=a(1,1)/b
end subroutine inv2



!Evaluate the Inverse of One complex Matrix with the determinat
subroutine inverse_d(a,n,det)
implicit none
integer,intent(IN)::n
complex(kind=8),intent(INOUT)::a(n,n)
complex(kind=8),intent(OUT)::det
complex(kind=8),allocatable::work(:)
complex(kind=8)::work_test(1)
integer::lwork
integer::ipiv(n),info
integer::i

call zgetrf(n,n,a,n,ipiv,info)
if (info < 0) then
  write(*,"(' zgetrf error info=',I4)") info
  write(*,"(' The',I4,'-th parameter is illegal.')") info
  stop
else if (info > 0) then
  write(*,"(' zgetrf error info=',I4)") info
  write(*,"(' Factorization is done. But the ',I4,'-th diagonal element is zero.')") info
  stop
! we can calculate determinant if everything is fine
else
  det=dcmplx(1.d0,0.d0)
  do i=1,n,1
    if (ipiv(i).ne.i) then
      det=dcmplx(-1.d0,0.d0)*det*a(i,i)
    else
      det=det*a(i,i)
    end if
  end do
end if


!get the best lwork
call zgetri(n,a,n,ipiv,work_test,-1,info)
lwork=work_test(1)!;lwork=10*n
allocate(work(lwork))


!get the inverse
call zgetri(n,a,n,ipiv,work,lwork,info)


!deallocate arrays
deallocate(work)


!Check the point
if (info < 0) then
  write(*,"(' zgetri error info=',I4)") info
  write(*,"(' The',I4,'-th parameter is illegal.')") info
  stop
else if (info > 0) then
  write(*,"(' zgetri error info=',I4)") info
  write(*,"(' The ',I4,'-th diagonal element is zero. Inversion cannot be performed.')") info
  stop
end if
end subroutine inverse_d


 
!Direct Multiply of Two Matrices
subroutine zdirectmul(a,m,b,n,beta,c)
implicit none
integer::m,n,i,j,k,l
complex(kind=8):: beta
complex*16::a(m,m),b(n,n)
complex*16::c(m*n,m*n)
 
  forall(i=1:m,j=1:m,k=1:n,l=1:n) ! index of b is inner loop
  c((i-1)*n+k,(j-1)*n+l)=beta*c((i-1)*n+k,(j-1)*n+l)+a(i,j)*b(k,l)
  end forall
end subroutine zdirectmul



! Diagonalization of One Real Symmetry Matrix
subroutine reigen(a,n,w)
implicit none
integer::n,info
real*8::a(n,n)
real*8::w(n)
real*8::work(10*n)
real*8::rwork(10*n)

call dsyev('V','U',n,a,n,w,work,10*n,info)
if(info/=0) then
   write(6,*)'DSYEV error ,info=',info
   stop
endif
end subroutine reigen
  

! Evaluate the Inverse of one real Matrix
subroutine rinverse(a,n)
implicit none
integer::n
real*8::a(n,n)
real*8::b(n,n),work(10*n)
real*8::ab(3*n+1,n)
integer::ipiv(n),info,ii
       
if (n<1) then
   print *, 'n wrong in rinverse',n
endif 

b=a
call dgetrf(n,n,b,n,ipiv,info)    
if(info/=0) then
   write(6,*)'DGETRF error ,info=,n=',info,n
   print *, a
stop
endif

call dgetri(n,b,n,ipiv,work,10*n,info)
if(info/=0) then
   write(6,*)'DGETRI error ,info=',info
   stop
endif

a=b
end subroutine rinverse


!check the Hermition of the matrix H
subroutine check_Hermite_c(H,n)
implicit none
integer, intent(IN):: n
complex(kind=8), intent(IN):: H(n,n)
real(kind=8):: error
integer:: i, j

error=0d0
do i=1, n
   do j=i+1,n
      error=error+ abs(H(i,j)-conjg(H(j,i)))
   enddo 
enddo 

if (error>1d-8) then 
   print *, 'H not Hermite', error
   stop
endif 

end subroutine check_Hermite_c

