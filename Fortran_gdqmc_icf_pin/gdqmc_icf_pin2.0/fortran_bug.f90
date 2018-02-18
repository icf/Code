module fortran_bug
contains

!use QR to do the modified GS
subroutine modGS(ph,L,N,det,Rmax)
implicit none
integer,intent(IN)::L,N
complex(kind=8),intent(INOUT)::ph(L,N)
real(kind=8),intent(OUT)::det
complex(kind=8),optional::Rmax(N,N)

complex(kind=8)::tau(N)
complex(kind=8),allocatable::work(:)
integer::lwork
integer::info

integer::i


allocate(work(1))
call zgeqrf(L,N,ph,L,tau,work,-1,info)
lwork=work(1)
deallocate(work)
if(info.NE.0) then
  write(*,*) "Something is wrong in QR:",info
end if


allocate(work(lwork))
call zgeqrf(L,N,ph,L,tau,work,lwork,info)
if(info.NE.0) then
  write(*,*) "Something is wrong in QR:",info
end if


det=1.d0
do i=1,N,1
   det=det/ph(i,i)
end do
call zungqr(L,N,N,ph,L,tau,work,lwork,info)
deallocate(work)


if(det.LT.0.d0) then
  det=-1.d0*det
  do i=1,L,1
     ph(i,1)=ph(i,1)*(-1.d0)
  end do
end if
end subroutine modGS


end module fortran_bug
