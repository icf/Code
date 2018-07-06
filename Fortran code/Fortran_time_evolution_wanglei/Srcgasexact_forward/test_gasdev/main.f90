program main
implicit none

integer,parameter:: Nsample=1000000
integer,parameter:: Nslice=100

real(kind=8),parameter:: Xmin=-2d0
real(kind=8),parameter:: Xmax=2d0

real(kind=8):: x, yl,yr
integer:: i,j
integer:: hist(Nslice)

real(kind=8),external :: gasdev
real(kind=8),external :: rannyu

hist=0

do i=1,Nsample
x=gasdev()

!x=rannyu()

!print *, i, x
!pause

do j=1,Nslice 
 yl=Xmin+(Xmax-Xmin)/dble(Nslice)*dble(j-1)
 yr=Xmin+(Xmax-Xmin)/dble(Nslice)*dble(j)
 if (x>=yl .and. x<yr) then 
 hist(j)=hist(j)+1
 endif 
enddo 
enddo 


open(11,file='hist_gasdev.dat')
do j=1,Nslice
write(11,*) Xmin+(Xmax-Xmin)/dble(Nslice)*dble(j-1),  dble(hist(j))/dble(Nsample)
enddo 
close(11)

end program 
