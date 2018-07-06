subroutine testn(Nold,Nnew,error,ave) 
implicit none

integer, intent(IN):: Nold, Nnew

real(kind=8):: w(Nold), phi_old(Nold), phi_new(Nnew)
real(kind=8):: deno, ave_before, ave_after, ave
real(kind=8):: error
integer:: i,table(Nnew)

real(kind=8),external:: rannyu

!real(kind=8),external:: r250

!call initialize_r250()

do i=1,Nold
w(i)=rannyu()
phi_old(i)=rannyu()
enddo 

ave_before=0d0
deno=0d0
do i=1,Nold
ave_before=ave_before+w(i)*phi_old(i)
deno=deno+w(i)
enddo 
ave_before=ave_before/deno


!call  reconfiguration(n,w,table)
call reconfiguration(Nold,Nnew,w,table)


do i=1, Nnew
phi_new(i)=phi_old(table(i))
print *, i,'-->', table(i)
pause
enddo 

ave_after=0d0
deno=0d0
do i=1,Nnew
ave_after=ave_after+phi_new(i)
deno=deno+1d0
enddo 
ave_after=ave_after/deno

error=abs(ave_before-ave_after)
ave= (ave_before+ave_after)*0.5d0

!do i=1,n
!print *, i, w(i) , table(i)
!enddo 

!print *, 'ave_before', ave_before
!print *, 'ave_after', ave_after

end subroutine 
