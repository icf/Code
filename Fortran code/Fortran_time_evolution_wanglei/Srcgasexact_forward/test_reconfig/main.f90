program main
implicit none

integer:: n
real(kind=8):: error,ave
real(kind=8),external:: rannyu

!real(kind=8),external:: r250


!call initialize_r250()


open(11,file='n_error.dat')
do n=1000, 100000, 500

call testn(n,n,error,ave)

write(11,*) n, error, ave
enddo 
close(11)

end program
