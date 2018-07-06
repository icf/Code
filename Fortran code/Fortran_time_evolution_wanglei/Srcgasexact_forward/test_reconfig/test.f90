program main
implicit none

integer,parameter:: n=50
real(kind=8):: w(n), phi(n), phi_temp(n)
real(kind=8):: deno, ave_before, ave_after
integer:: i,table(n)


real(kind=8),external:: r250

call initialize_r250()

do i=1,n
w(i)=r250()
phi(i)=r250()
enddo 

ave_before=0d0
deno=0d0
do i=1,n
ave_before=ave_before+w(i)*phi(i)
deno=deno+w(i)
enddo 
ave_before=ave_before/deno


call  reconfiguration(n,w,table)

phi_temp=phi
do i=1, n
phi(i)=phi_temp(table(i))
w(i)=1d0
enddo 

ave_after=0d0
deno=0d0
do i=1,n
ave_after=ave_after+w(i)*phi(i)
deno=deno+w(i)
enddo 
ave_after=ave_after/deno


!do i=1,n
!print *, i, w(i) , table(i)
!enddo 

print *, 'ave_before', ave_before
print *, 'ave_after', ave_after

end program
