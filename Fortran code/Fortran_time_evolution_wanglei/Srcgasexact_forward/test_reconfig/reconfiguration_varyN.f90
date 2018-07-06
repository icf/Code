subroutine reconfiguration(Nold,Nnew,w,table)
implicit none

integer,intent(IN):: Nold,Nnew
real(kind=8),intent(IN):: w(Nold)

integer,intent(OUT):: table(Nnew)

real(kind=8):: wsum,p(Nold),z
integer:: i,j

real(kind=8),external:: rannyu
!real(kind=8),external:: r250

wsum=0d0
do j=1,Nold
wsum=wsum+w(j)
enddo 

!---Prepare the probability
do j=1,Nold
p(j)=w(j)/wsum
!print *, 'j,p', j,p(j)
enddo 


!---initialize-----------
j=1
wsum=p(1)

do i=1,Nnew

z=dble((i-1)+rannyu())/dble(Nnew)
!z=dble((i-1)+0.1d0)/dble(n)

do while(z>wsum)
!print *, i,j,z,wsum
!pause
 j=j+1
 wsum=wsum+p(j) 
enddo 

table(i)=j
!print *, 'table' ,i,'=',j
enddo 

end subroutine reconfiguration
