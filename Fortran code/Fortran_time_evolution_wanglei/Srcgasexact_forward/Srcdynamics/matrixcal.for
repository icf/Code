!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	   Diagonalization of One Hermitian Matrix
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eigen(a,n,w)
	 implicit none
       integer::n,info
	 complex*16::a(n,n)
	 real*8::w(n)
       complex*16::work(10*n)
       real*8::rwork(10*n)
	
	 call zheev('V','U',n,a,n,w,work,10*n,rwork,info)
       if(info/=0) then
        write(6,*)'ZHEEV error ,info=',info
        stop
       endif
      end subroutine eigen
  
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
  
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	   Evaluate the Inverse of One Matrix
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inverse(a,n)
       implicit none
       integer::n
       complex*16::a(n,n)
       complex*16::b(n,n),work(10*n)
       integer::ipiv(n),info
       b=a
       call zgetrf(n,n,b,n,ipiv,info)
       if(info/=0) then
        write(6,*)'ZGETRF error ,info=',info
        stop
       endif
       call zgetri(n,b,n,ipiv,work,10*n,info)
       if(info/=0) then
        write(6,*)'ZGETRI error ,info=',info
        stop
       endif
       a=b
      end subroutine inverse
 
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	   Direct Multiply of Two Matrices
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine directmul(a,m,b,n,c)
       implicit none
       integer::m,n,i,j,k,l
       complex*16::a(m,m),b(n,n)
       complex*16::c(m*n,m*n)
 
       forall(i=1:m,j=1:m,k=1:n,l=1:n)
         c((i-1)*n+k,(j-1)*n+l)=a(i,j)*b(k,l)
       end forall
      end subroutine directmul



!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	   Diagonalization of One Real Symmetry Matrix
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	   Evaluate the Inverse of One Matrix
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
