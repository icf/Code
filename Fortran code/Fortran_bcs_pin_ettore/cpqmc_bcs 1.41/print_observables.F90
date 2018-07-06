subroutine print_observables(quando)
use param
use mc_loop_param
use mpi_serial_param
use adET_param
use meas_param
use mpi_serial_param
use method_param
use meas_param
use lattice_param
use model_param
use adET_param
implicit none
integer, intent(IN)::quando

integer :: sitei,ipair,iunit,beta,alpha,jb,ib,i_beta,k,m,n
integer::kl(Dimen)
real(kind=8)::ph
real(kind=8) :: x,y
complex(kind=8)::kk
character*3 :: sfix
character*20 :: filename



end subroutine print_observables


subroutine xifs(sfix,i)
implicit none
character*3 sfix
integer i
if(i.lt.10)then
  write(sfix,'(i1)')i
elseif(i.lt.100)then
  write(sfix,'(i2)')i
elseif(i.lt.1000)then
  write(sfix,'(i3)')i
endif
return
end
