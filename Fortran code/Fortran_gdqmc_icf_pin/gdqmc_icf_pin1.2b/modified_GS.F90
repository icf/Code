!This subroutine do the Modified GS periodically.
subroutine Modified_GS()
use lattice_param
use model_param
use phi_param
use mpi_serial_param
use mc_loop_param
use method_param
implicit none
integer::i

do i=1,Nwalkers,1
   call modGS_i(i)
end do


if(crn.LT.0.d0)  then !cpmc call the imp
  call cal_imp_ovlap()
else
  !We want to check the abs_imp_avg, imp_avg, imp_err, so calculate 
  !imp all the time even when ccoe=0.d0
  !if(abs(ccoe).GT.1d-8) then !fpmc ccoe.ne.0 call the imp
    call cal_imp_ovlap()
  !end if
end if
end subroutine Modified_GS



subroutine modGS_i(i)
use lattice_param
use model_param
use phi_param
use mpi_serial_param
use mc_loop_param
use method_param

use fortran_bug
implicit none
integer,intent(IN)::i
real(kind=8)::anm,anm1,anm2
complex(kind=8)::Rmat(Ntot,Ntot)

if(dtype.EQ.'c') then
  call modGS(phi(1,1,i),2*Nsite,Ntot,anm,Rmat)
else if(dtype.EQ.'d') then
  call modGS(phi(1:Nsite,1:Nspin(1),i),Nsite,Nspin(1),anm1,Rmat(1:Nspin(1),1:Nspin(1)))
  call modGS(phi((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,i),Nsite,Nspin(2),anm2,Rmat(1:Nspin(2),1:Nspin(2)))
  anm=anm1*anm2
end if


!For the release or free projection 
if(crn.GT.0.d0) then
  if(anm.LT.0.d0) then
    write(*,*) "Something is wrong in Modified_GS."
    write(*,*) "anm:",anm
    call mystop
  end if
  if(weight(i).le.0.d0) return
  weight(i)=weight(i)/anm
  dlogw(i)=dlogw(i)-dlog(anm)
end if
end subroutine modGS_i



!Get the total impfunc and tot_imp
subroutine cal_imp_ovlap()
use param
use mc_loop_param
use phi_param
use mpi_serial_param

use sc_loop_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer::i
complex(kind=8)::abs_imp,imp,imp2,imp_all
abs_imp=zero;imp=zero;imp2=zero
do i=1,Nwalkers,1
   call get_imp(i)
   !tot_imp(i)=dcmplx(i+rank*Nwalkers,i+rank*Nwalkers-100)
   abs_imp=abs_imp+abs(tot_imp(i))
   imp=imp+tot_imp(i)
   imp2=imp2+tot_imp(i)*conjg(tot_imp(i))
end do

#ifdef MPI
call MPI_ALLREDUCE(abs_imp,imp_all,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
abs_imp=imp_all
call MPI_ALLREDUCE(imp,imp_all,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
imp=imp_all
call MPI_ALLREDUCE(imp2,imp_all,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
imp2=imp_all
#endif

abs_imp=abs_imp/dble(Nwalkers*Nsize)
imp=imp/dble(Nwalkers*Nsize)
imp2=imp2/dble(Nwalkers*Nsize)


!Complex defination of variance=E[(X-u)*(X-u)^{+}]=E(X*X^+)-u*u^+
abs_imp_avg=abs_imp
imp_avg=imp
imp_err=sqrt((imp2-imp*conjg(imp))/dble(Nwalkers*Nsize))

if(rank.eq.0) then
  write(*,*) "The average of abs important function after gs is:",abs_imp_avg
  write(*,*) "The average total important function with error bar after gs:",imp_avg,imp_err
end if

!write(*,*) tot_imp
!call mystop
end subroutine cal_imp_ovlap



!Get impfunc(k) and tot_imp(k)
subroutine get_imp(k)
use param
use phiT_param
use phi_param
use model_param
use mpi_serial_param
use lattice_param

use sc_loop_param
implicit none
integer,intent(IN)::k
complex(kind=8)::imp
integer::m
tot_imp(k)=zero
do m=1,Dtot,1
   call imp_fun_dc(phiT(1,1,m),phi(1,1,k),imp)
   !if(abs(imp).LT.1d-15) then
   !   write(*,*) "Small overlap phiT with phi",k,m,rank
   !   write(*,*) "imptant function is:",imp
   !  !call mystop
   !end if
   impfunc(m,k)=imp
   tot_imp(k)=tot_imp(k)+conjg(coe_multi(m))*impfunc(m,k)
end do
end subroutine get_imp


!Get impfunc(1:Dtot,k), tot_imp(k) and invop(1:Ntot,1:Ntot,1:Dtot)
subroutine get_imp_inv(k,invop)
use param
use phiT_param
use phi_param
use model_param
use mpi_serial_param
use lattice_param
implicit none
integer,intent(IN)::k
complex(kind=8),intent(OUT)::invop(Ntot,Ntot,Dtot)
complex(kind=8)::imp
integer::m
tot_imp(k)=zero
do m=1,Dtot,1
   !call imp_fun_dc(phiT(1,1,m),phi(1,1,k),imp)
   call over_lap_dc(phiT(1,1,m),phi(1,1,k),invop(1,1,m))
   call inverse_d_dc(invop(1,1,m),imp)
   !if(abs(imp).LT.1d-15) then
   !   write(*,*) "Small overlap phiT with phi",k,m,rank
   !   write(*,*) "imptant function is:",imp
   !  call mystop
   !end if
   impfunc(m,k)=imp
   tot_imp(k)=tot_imp(k)+conjg(coe_multi(m))*impfunc(m,k)
end do
end subroutine get_imp_inv




!use QR to do the modified GS
!subroutine modGS(ph,L,N,det,Rmax)
!in fortran_bug



!subroutine test()
!use param
!implicit none
!integer::i,j,k,l
!complex(kind=8)::t_up(Nup,Nup,Dup),t_dn(Ndn,Ndn,Ddn),tmp,tot
!complex(kind=8)::iup(Dup),idn(Ddn)
!do i=1,Nwalkers,1
!
!if(crn.LT.0.d0) then
!   t_up(1:Nup,1:Nup,1:Dup)=ovlap_up(1:Nup,1:Nup,1:Dup,i)
!   t_dn(1:Ndn,1:Ndn,1:Ddn)=ovlap_dn(1:Ndn,1:Ndn,1:Ddn,i)
!   iup(1:Dup)=impfunc0_up(1:Dup,i)
!   idn(1:Ddn)=impfunc0_dn(1:Ddn,i)
!   tot=tot_imp(i)
!
!
!   call cal_imp_ovlap()
!
!
!   tmp=zero
!   do l=1,Nup,1
!      do j=1,Nup,1
!         do k=1,Dup,1
!            tmp=tmp+t_up(l,j,k)-ovlap_up(l,j,k,i)
!         end do
!      end do
!   end do
!   if(abs(tmp).GT.1D-10) then
!     write(*,*) "wrong with ov up",tmp,i,rank
!     !write(*,*) "----------------------"
!     if(rank.eq.1) write(*,*) t_up(:,:,:)
!     if(rank.eq.1) write(*,*) "----------------------"
!     !if(rank.eq.1) write(*,*) ovlap_up(:,:,:,i)
!     !if(rank.eq.1) write(*,*) "----------------------"
!     call mystop
!   end if
!
!   tmp=zero
!   do l=1,Ndn,1
!      do j=1,Ndn,1
!         do k=1,Ddn,1
!            tmp=tmp+t_dn(l,j,k)-ovlap_dn(l,j,k,i)
!         end do
!      end do
!   end do
!   if(abs(tmp).GT.1D-10) then
!     write(*,*) "wrong with ov dn",tmp,i,rank
!     call mystop
!   end if
!
!   if(abs(tot-tot_imp(i)).GT.1d-10) then
!     write(*,*) "tot wrong",tot-tot_imp(i),i,rank
!     call mystop
!   end if
!
!
!   tmp=zero
!   do l=1,Dup,1
!      tmp=tmp+iup(l)-impfunc0_up(l,i)
!   end do
!   if(abs(tmp).GT.1D-10) then
!     write(*,*) "wrong with imp up",tmp,i,rank
!     call mystop
!   end if
!
!   tmp=zero
!   do l=1,Ddn,1
!      tmp=tmp+idn(l)-impfunc0_dn(l,i)
!   end do
!   if(abs(tmp).GT.1D-10) then
!     write(*,*) "wrong with imp dn",tmp,i,rank
!     call mystop
!   end if
!end if
!
!end do
!end subroutine test
