!-------------------------------------------------------------------
!Do the data manipulate of the step by step qmc, write into the file
!-------------------------------------------------------------------
subroutine data_step_mc()
use io_module
use meas_param
use mc_loop_param
use project_param
use lattice_param
use sc_loop_param
implicit none
integer::i,j,sitei,sitej
real(kind=8)::e_t,e_e,v_t,v_e,k_t,k_e
complex(kind=8)::var_t,var_e
real(kind=8)::nu_t,nu_e,nd_t,nd_e
real(kind=8)::mean,err
complex(kind=8)::err_c
real(kind=8)::pt
call get_filename()
call openUnit(EnergyName,10,'R')
call openUnit(ScorrName,20,'R')
call openUnit(SkName,30,'R')
call openUnit(ckName,40,'R')
call openUnit(numName,50,'R')
call openUnit(ccName,60,'R')
do i=1,max_local,1
   pt=dble(Thermblock*blockstep+i*meastep)*dt

   call err_anal(kin_l(1:Nsamples,i),Nsamples,k_t,k_e)
   call err_anal(v_l(1:Nsamples,i),Nsamples,v_t,v_e)
   call err_anal(e_l(1:Nsamples,i),Nsamples,e_t,e_e)
   call err_anal_c(var_l(1:Nsamples,i),Nsamples,var_t,var_e)
   write(10,'(10f15.8)') pt,k_t,k_e,v_t,v_e,e_t,e_e,abs(var_t),abs(var_e),abs(sig(i)/absig(i))
   
   call err_anal(nu_l(1:Nsamples,i),Nsamples,nu_t,nu_e)
   call err_anal(nd_l(1:Nsamples,i),Nsamples,nd_t,nd_e)
   write(50,'(5f15.8)') pt,nu_t,nu_e,nd_t,nd_e

   do j=1,Nsite,1
      call err_anal(sisj_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(20,'(I4,3f15.8)') j,pt,mean,err
   end do
   call fourier_sij(i)
   do j=1,Nsite,1
      call err_anal(sk_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(30,'(I4,3f15.8)') j,pt,mean,err
   end do
   call fourier_cij(i)
   do j=1,2*Nsite,1
      call err_anal(ck_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(40,'(I4,3f15.8)') j,pt,mean,err
   end do
   do sitei=1,2*Nsite,1
      do sitej=1,2*Nsite,1
         call err_anal_c(cicj_l_global(1:Nsamples,sitei,sitej,i),Nsamples,cicj_sc_global(sitei,sitej),err_c)
      enddo
      write(60,'(I4,3f15.8)') sitei,pt,abs(cicj_sc_global(sitei,sitei)),abs(err_c)  
   end do

end do
close(10)
close(20)
close(30)
close(40)
close(50)
close(60)
end subroutine data_step_mc



!-------------------------------------------------------------
!Get the mean value of the MC measurement of different i_local
!-------------------------------------------------------------
subroutine mean_meas()
use method_param
use lattice_param
use mc_loop_param
use meas_param
implicit none
integer::max_n
integer::sample,i,sitei,sitej
max_n=max(max_crn,0)
if(fw_bk.NE.'FW') max_n=0 !for release back propogation.
do sample=1,Nsamples,1
   do i=0,max_n,1
      kin_l(sample,i)=kin_l(sample,i)/dble(max_local)
      v_l(sample,i)=v_l(sample,i)/dble(max_local)
      e_l(sample,i)=e_l(sample,i)/dble(max_local)
      var_l(sample,i)=var_l(sample,i)/dble(max_local)
      nu_l(sample,i)=nu_l(sample,i)/dble(max_local)
      nd_l(sample,i)=nd_l(sample,i)/dble(max_local)
      do sitei=1,Nsite,1
         sisj_l(sample,sitei,i)=sisj_l(sample,sitei,i)/dble(max_local)
      end do
      do sitei=1,2*Nsite,1
         cicj_l(sample,sitei,i)=cicj_l(sample,sitei,i)/max_local
      end do
      do sitei=1,2*Nsite,1
         do sitej=1,2*Nsite,1
            cicj_l_global(sample,sitei,sitej,i)=cicj_l_global(sample,sitei,sitej,i)/max_local
         enddo
      end do
   end do
end do
!write(*,*) e_l
end subroutine mean_meas



!-----------------------------------------------------------------
!Do the data manipulate of the cpmc and rcpmc, write into the file
!-----------------------------------------------------------------
subroutine data_cpmc_rcpmc()
use io_module
use meas_param
use project_param
use mc_loop_param
use method_param
use lattice_param
use mpi_serial_param
use sc_loop_param
implicit none
integer::i,j,sitei,sitej
real(kind=8)::e_t,e_e,v_t,v_e,k_t,k_e
complex(kind=8)::var_t,var_e
real(kind=8)::nu_t,nu_e,nd_t,nd_e
real(kind=8)::mean,err
complex(kind=8)::err_c
real(kind=8)::pt
integer::max_n
call get_filename()
call openUnit(EnergyName,10,'R')
call openUnit(ScorrName,20,'R')
call openUnit(SkName,30,'R')
call openUnit(ckName,40,'R')
call openUnit(numName,50,'R')
call openUnit(ccName,60,'R')
call openUnit(SzName,70,'R')
max_n=max(max_crn,0)
if(fw_bk.NE.'FW') max_n=0 !for release back propogation.
do i=0,max_n,1
   pt=dble(i)*dt
   call err_anal(kin_l(1:Nsamples,i),Nsamples,k_t,k_e)
   call err_anal(v_l(1:Nsamples,i),Nsamples,v_t,v_e)
   call err_anal(e_l(1:Nsamples,i),Nsamples,e_t,e_e)
   call err_anal_c(var_l(1:Nsamples,i),Nsamples,var_t,var_e)
   
   write(10,'(10f15.8)') pt,k_t,k_e,v_t,v_e,e_t,e_e,abs(var_t),abs(var_e),abs(sig(i)/absig(i))

   call err_anal(nu_l(1:Nsamples,i),Nsamples,nu_t,nu_e)
   call err_anal(nd_l(1:Nsamples,i),Nsamples,nd_t,nd_e)
   write(50,'(5f15.8)') pt,nu_t,nu_e,nd_t,nd_e

   do j=1,Nsite,1
      call err_anal(sisj_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(20,'(I4,3f15.8)') j,pt,mean,err
   end do 
   call fourier_sij(i)
   do j=1,Nsite,1
      call err_anal(sk_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(30,'(I4,3f15.8)') j,pt,mean,err
   end do
   call fourier_cij(i)
   do j=1,2*Nsite,1
      call err_anal(ck_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(40,'(I4,3f15.8)') j,pt,mean,err
   end do
   do sitei=1,2*Nsite,1
      do sitej=1,2*Nsite,1
         call err_anal_c(cicj_l_global(1:Nsamples,sitei,sitej,i),Nsamples,cicj_sc_global(sitei,sitej),err_c)
      enddo
   end do

   do sitei=1,Nsite,1
      write(60,'(I4,3f15.8)') sitei,pt,1-abs(cicj_sc_global(sitei,sitei))-abs(cicj_sc_global(sitei+Nsite,sitei+Nsite))
      write(70,'(I4,3f15.8)') sitei,pt,((-1)**(coor(sitei,1)+coor(sitei,2)))* &
(abs(cicj_sc_global(sitei,sitei))-abs(cicj_sc_global(sitei+Nsite,sitei+Nsite)))
   enddo
end do
close(10)
close(20)
close(30)
close(40)
close(50)
close(60)
close(70)
end subroutine data_cpmc_rcpmc



!---------------------------
!Get the error bar of dat(N)
!---------------------------
subroutine err_anal(dat,N,m,er)
implicit none
integer,intent(IN)::N
real(kind=8),intent(IN)::dat(N)
real(kind=8),intent(OUT)::m,er
integer::i,j,k

if(N.LE.0) then
  write(*,*) "N should not be smaller than or EQ 0", N
  call mystop
else if(N.EQ.1) then
  !write(*,*) "N eq 1 warning", N
  m=dat(1)
  er=0.d0
  return
end if

!Get the mean
m=0.d0
do i=1,N,1
   m=m+dat(i)
end do
m=m/dble(N)

er=0.d0
do i=1,N,1
   er=er+dat(i)**2
end do
er=er/dble(N)

er=er-m**2

if(ABS(er).GT.1.d-10) then
  if(er.lt.0.d0) then
    write(*,*) "Something is wrong in err_anal",er
  end if
else
  er=0.d0
end if

er=sqrt(er/dble(N-1))
end subroutine err_anal



!------------------------------------
!Get the error bar of dat(N)--complex
!u=E[u]
!variance=E[(x-u)(x-u)*]
!error bar=sqrt[variance/(N-1)]
!------------------------------------
subroutine err_anal_c(dat,N,m,er)
implicit none
integer,intent(IN)::N
complex(kind=8),intent(IN)::dat(N)
complex(kind=8),intent(OUT)::m,er
integer::i,j,k

if(N.LE.0) then
  write(*,*) "N should not be smaller than or EQ 0", N
  call mystop
else if(N.EQ.1) then
  !write(*,*) "N eq 1 warning", N
  m=dat(1)
  er=dcmplx(0.d0,0.d0)
  return
end if

!Get the mean
m=dcmplx(0.d0,0.d0)
do i=1,N,1
   m=m+dat(i)
end do
m=m/dble(N)

er=dcmplx(0.d0,0.d0)
do i=1,N,1
   er=er+abs(dat(i)-m)**2
end do
er=er/dble(N)

if(ABS(er).GT.1.d-10) then
  !if(er.lt.0.d0) then
  !  write(*,*) "Something is wrong in err_anal",er
  !end if
else
  er=dcmplx(0.d0,0.d0)
end if

er=sqrt(er/dble(N-1))
end subroutine err_anal_c
