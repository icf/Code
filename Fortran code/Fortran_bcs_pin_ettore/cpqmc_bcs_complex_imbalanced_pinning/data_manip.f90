!-------------------------------------------------------------------
!Do the data manipulate of the step by step qmc, write into the file
!-------------------------------------------------------------------
subroutine data_step_mc()
use io_module
use meas_param
use mc_loop_param
use project_param
use lattice_param
use method_param
implicit none
integer::i,j,i_beta,ipair,ib,jb
real(kind=8)::e_t,e_e,v_t,v_e,k_t,k_e
complex(kind=8)::var_t,var_e,zmean
real(kind=8)::nu_t,nu_e,nd_t,nd_e
real(kind=8)::mean,err
real(kind=8)::pt
call get_filename()
call openUnit(EnergyName,10,'R')
call openUnit(ScorrName,20,'R')
call openUnit(SkName,30,'R')
!call openUnit(ckName,40,'R')
call openUnit(numName,50,'R')
do i=1,max_local,1
   pt=dble(Thermblock*blockstep+i*meastep)*dt

   call err_anal(kin_l(1:Nsamples,i),Nsamples,k_t,k_e)
   call err_anal(v_l(1:Nsamples,i),Nsamples,v_t,v_e)
   call err_anal(e_l(1:Nsamples,i),Nsamples,e_t,e_e)
   call err_anal_c(var_l(1:Nsamples,i),Nsamples,var_t,var_e)
   write(10,'(10f15.8)') pt,k_t,k_e,v_t,v_e,e_t,e_e!,abs(var_t),abs(var_e),abs(sig(i)/absig(i))

   call err_anal(nu_l(1:Nsamples,i),Nsamples,nu_t,nu_e)
   call err_anal(nd_l(1:Nsamples,i),Nsamples,nd_t,nd_e)
   write(50,'(5f15.8)') pt,nu_t,nu_e,nd_t,nd_e

   do j=1,Nbravais,1
      call err_anal(sisj_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(20,'(I4,3f15.8)') j,pt,mean,err
   end do
   call fourier_sij(i)
   do j=1,Nsite,1
      call err_anal(sk_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(30,'(I4,3f15.8)') j,pt,mean,err
   end do
!   call fourier_cij(i)
!   do j=1,Nbravais,1
!      call err_anal(ck_l(1:Nsamples,j,i),Nsamples,mean,err)
!      write(40,'(I4,3f15.8)') j,pt,mean,err
!   end do
!   do j=Nsite+1,Nsite+Nbravais,1
!      call err_anal(ck_l(1:Nsamples,j,i),Nsamples,mean,err)
!      write(40,'(I4,3f15.8)') j,pt,mean,err
!   end do
   do ipair=1,Npair
     do j=1,Nbravais,1
       call err_anal_complex(didj_l(1:Nsamples,j,ipair,i),Nsamples,zmean,err)
     enddo
   enddo
   if(I_onebf.eq.1)then
   do i_beta=1,Nbeta,1
      do j=1,2*Nsite,1
         call err_anal_complex(cicj_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
      enddo
   end do
   endif
   if(Nbeta.gt.0)then
     if(I_twob.eq.1)then
       do i_beta=0,Nbeta,1
         do j=1,Nbravais,1
       !  call err_anal_complex(rho_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
           call err_anal_complex(nupnup_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
           call err_anal_complex(ndnnup_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
         enddo
       end do
     endif
   endif
end do
close(10)
close(20)
close(30)
close(40)
close(50)
end subroutine data_step_mc



!-------------------------------------------------------------
!Get the mean value of the MC measurement of different i_local
!-------------------------------------------------------------
subroutine mean_meas()
use method_param
use lattice_param
use mc_loop_param
use meas_param
use model_param
implicit none
integer::max_n
integer::sample,i,sitei,sitej,i_beta,ipair,ib,jb,alpha,beta,idir
max_n=max(max_crn,0)
if(fw_bk.NE.'FW') max_n=0 !for release back propogation.
do sample=1,Nsamples,1
   do i=0,max_n,1
      do idir=1,Dimen
        ZetaN_l(sample,idir,i)=ZetaN_l(sample,idir,i)/dble(max_local)
      enddo
      kin_l(sample,i)=kin_l(sample,i)/dble(max_local)
      v_l(sample,i)=v_l(sample,i)/dble(max_local)
      e_l(sample,i)=e_l(sample,i)/dble(max_local)
      eE_l(sample,i)=eE_l(sample,i)/dble(max_local)
      var_l(sample,i)=var_l(sample,i)/dble(max_local)
      nu_l(sample,i)=nu_l(sample,i)/dble(max_local)
      nd_l(sample,i)=nd_l(sample,i)/dble(max_local)
      if(I_obdm.eq.1)then
        do sitei=1,2*Nsite
          do sitej=1,2*Nsite
            obdm_l(sample,sitej,sitei,i)=obdm_l(sample,sitej,sitei,i)/dble(max_local)
          enddo
        enddo
      endif
      if(ipinn.eq.1)then
      do beta=1,2
        do alpha=1,2
          do sitei=1,Nsite,1
            nofr_l(sample,sitei,alpha,beta,i)=nofr_l(sample,sitei,alpha,beta,i)/max_local
          enddo
        enddo
      enddo
      endif
      if(ipinn.eq.0)then
!      do ib=1,Nbands
!        do jb=1,Nbands
!          do sitei=1,Nbravais
!            cacb_l(sample,sitei,jb,ib,i)=cacb_l(sample,sitei,jb,ib,i)/max_local
!          enddo
!        enddo
!      enddo
      do sitei=1,Nbravais,1
         sisj_l(sample,sitei,i)=sisj_l(sample,sitei,i)/dble(max_local)
      end do
      do ib=1,Nbands
        do jb=1,Nbands
         do sitei=1,Nbravais
           ninj_l(sample,sitei,jb,ib,i)=ninj_l(sample,sitei,jb,ib,i)/max_local
         enddo
       enddo
      enddo
      do ib=1,Nbands
        do jb=1,Nbands
         do sitei=1,Nbravais
           szsz_l(sample,sitei,jb,ib,i)=szsz_l(sample,sitei,jb,ib,i)/max_local
         enddo
       enddo
      enddo

      endif
!      do sitei=1,2*Nsite,1
!         cicj_l(sample,sitei,i)=cicj_l(sample,sitei,i)/max_local
!      end do
      do ipair=1,Npair
        do sitei=1,Nbravais,1
          didj_l(sample,sitei,ipair,i)=didj_l(sample,sitei,ipair,i)/max_local
        end do
      end do
      if(I_onebf.eq.1)then
        do i_beta=1,Nbeta,1
          do sitei=1,2*Nsite,1
            cicj_t_l(sample,sitei,i_beta,i)=cicj_t_l(sample,sitei,i_beta,i)/max_local
          end do
        end do
        do i_beta=1,Nbeta,1
          do sitei=1,2*Nsite,1
            cicjh_t_l(sample,sitei,i_beta,i)=cicjh_t_l(sample,sitei,i_beta,i)/max_local
          end do
        end do
      endif
      if(Nbeta.gt.0)then
        if(I_twob.eq.1)then
          do i_beta=0,Nbeta,1
            do sitei=1,Nbravais,1
         ! rho_t_l(sample,sitei,i_beta,i)=rho_t_l(sample,sitei,i_beta,i)/max_local
              nupnup_t_l(sample,sitei,i_beta,i)=nupnup_t_l(sample,sitei,i_beta,i)/max_local
              ndnnup_t_l(sample,sitei,i_beta,i)=ndnnup_t_l(sample,sitei,i_beta,i)/max_local
            end do
          end do
        endif
      do i_beta=0,Nbeta,1
        GreenP_t_l(sample,i_beta,i)=GreenP_t_l(sample,i_beta,i)/max_local
      end do
      do i_beta=0,Nbeta,1
        GreenH_t_l(sample,i_beta,i)=GreenH_t_l(sample,i_beta,i)/max_local
      end do
      endif
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
use model_param
use mpi_serial_param
!use phiT_param
implicit none
integer::i,j,i_beta,ipair,ib,jb,sitei,sitej,beta,alpha,idir
real(kind=8)::e_t,e_e,v_t,v_e,k_t,k_e,eE_T,eE_e
complex(kind=8)::var_t,var_e,zmean
real(kind=8)::nu_t,nu_e,nd_t,nd_e
real(kind=8)::mean,err
real(kind=8)::pt
integer::max_n
call get_filename()
call openUnit(EnergyName,10,'R')
call openUnit(ZetaNName,11,'R')
call openUnit(ObdmName,15,'R')
call openUnit(ScorrName,20,'R')
call openUnit(NcorrName,25,'R')
call openUnit(SzcorrName,26,'R')
call openUnit(SkName,30,'R')
!call openUnit(ckName,40,'R')
!call openUnit(cabName,45,'R')
call openUnit(nofrName,35,'R')
call openUnit(numName,50,'R')
if(Nbeta.gt.0)then
  if(I_onebf.eq.1)call openUnit(cijtName,60,'R')
  if(I_onebf.eq.1)call openUnit(cijhtName,61,'R')
!  call openUnit(rhotName,61,'R')
  call openUnit(GreenPtName,64,'R')
  call openUnit(GreenHtName,65,'R')
  if(I_twob.eq.1)call openUnit(nupnuptName,62,'R')
  if(I_twob.eq.1)call openUnit(ndnnuptName,63,'R')
endif
if(Npair.gt.0)call openUnit(DcorrName,70,'R')
max_n=max(max_crn,0)
if(fw_bk.NE.'FW') max_n=0 !for release back propogation.
do i=0,max_n,1
   pt=dble(i)*dt
   do idir=1,Dimen
     call err_anal_complex(ZetaN_l(1:Nsamples,idir,i),Nsamples,zmean,err)
     write(11,'(4f15.8)')pt,dble(zmean),aimag(zmean),err
   enddo
   call err_anal(kin_l(1:Nsamples,i),Nsamples,k_t,k_e)
   call err_anal(v_l(1:Nsamples,i),Nsamples,v_t,v_e)
   call err_anal(e_l(1:Nsamples,i),Nsamples,e_t,e_e)
   call err_anal(eE_l(1:Nsamples,i),Nsamples,eE_t,eE_e)
!   call err_anal_c(var_l(1:Nsamples,i),Nsamples,var_t,var_e)
   write(10,'(10f15.8)') pt,k_t,k_e,v_t,v_e,e_t,e_e,eE_t,eE_e !,abs(var_t),abs(var_e),abs(sig(i)/absig(i))

   call err_anal(nu_l(1:Nsamples,i),Nsamples,nu_t,nu_e)
   call err_anal(nd_l(1:Nsamples,i),Nsamples,nd_t,nd_e)
   write(50,'(5f15.8)') pt,nu_t,nu_e,nd_t,nd_e

   if(I_obdm.eq.1)then
   do sitei=1,2*Nsite
     do sitej=1,2*Nsite
       call err_anal_complex(obdm_l(1:Nsamples,sitej,sitei,i),    &
     &                       Nsamples,zmean,err)
!       if(I_wavefun.eq.1)then
!         write(15,'(2I4,3f15.8)')sitej,sitei,dble(zmean),aimag(zmean),err
!       elseif(I_wavefun.eq.2)then 
         if(sitei.le.Nsite.and.sitej.le.Nsite)then
           write(15,'(2I4,3f15.8)')sitej,sitei,dble(zmean),aimag(zmean),err
         elseif(sitei.gt.Nsite.and.sitej.gt.Nsite)then
           write(15,'(2I4,3f15.8)')sitej,sitei,dble(zmean),aimag(zmean),err
         else
           write(15,'(2I4,3f15.8)')sitej,sitei,0.d0,0.d0,0.d0
         endif
!       endif
     enddo
   enddo
   endif


!   if(ipinn.eq.0)then
!   do ib=1,Nbands
!     do jb=1,Nbands
!       do sitei=1,Nbravais
!         call err_anal_complex(cacb_l(1:Nsamples,sitei,jb,ib,i),    &
!     &                         Nsamples,zmean,err)
!         write(45,'(3I4,3f15.8)')sitei,jb,ib,dble(zmean),aimag(zmean),err
!       enddo
!     enddo
!   enddo 
!   endif

   if(ipinn.eq.1)then
   do beta=1,2
     do alpha=1,2
       do sitei=1,Nsite
         call err_anal_complex(nofr_l(1:Nsamples,sitei,alpha,beta,i),    &
     &                         Nsamples,zmean,err)
         write(35,'(3I4,3f15.8)')sitei,alpha,beta,dble(zmean),aimag(zmean),err
       enddo
     enddo
   enddo
   endif

   if(ipinn.eq.0)then
   do j=1,Nbravais,1
      call err_anal(sisj_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(20,'(I4,3f15.8)') j,pt,mean,err
   end do 
   call fourier_sij(i)
   do beta=1,Nbands
     do alpha=1,Nbands
       do sitei=1,Nbravais
         call err_anal(ninj_l(1:Nsamples,sitei,alpha,beta,i),    &
     &                 Nsamples,mean,err)
         write(25,'(I4,3f15.8,2I4)')sitei,pt,mean,err,alpha,beta
       enddo
       write(25,*)
       write(25,*)
     enddo
   enddo
   do beta=1,Nbands
     do alpha=1,Nbands
       do sitei=1,Nbravais
         call err_anal(szsz_l(1:Nsamples,sitei,alpha,beta,i),    &
     &                 Nsamples,mean,err)
         write(26,'(I4,3f15.8,2I4)')sitei,pt,mean,err,alpha,beta
       enddo
       write(26,*)
       write(26,*)
     enddo
   enddo
   do j=1,Nbravais,1
      call err_anal(sk_l(1:Nsamples,j,i),Nsamples,mean,err)
      write(30,'(I4,3f15.8)') j,pt,mean,err
   end do
   endif

!   call fourier_cij(i)
!   do j=1,Nbravais,1
!      call err_anal(ck_l(1:Nsamples,j,i),Nsamples,mean,err)
!      write(40,'(I4,3f15.8)') j,pt,mean,err     
!   end do
!   write(40,*)
!   write(40,*)
!   do j=Nsite+1,Nsite+Nbravais,1
!      call err_anal(ck_l(1:Nsamples,j,i),Nsamples,mean,err)
!      write(40,'(I4,3f15.8)') j,pt,mean,err
!   end do
!   do j=1,2*Nsite,1
!      call err_anal(ck_l(1:Nsamples,j,i),Nsamples,mean,err)
!      write(40,'(I4,3f15.8)') j,pt,mean,err
!   end do
   do ipair=1,Npair
     do j=1,Nbravais,1
       call err_anal_complex(didj_l(1:Nsamples,j,ipair,i),Nsamples,zmean,err)
       write(70,'(I4,4f15.8)')j,pt,dble(zmean),aimag(zmean),err
     enddo
     write(70,*)
     write(70,*)
   enddo
   if(I_onebf.eq.1)then
     do j=1,2*Nsite,1
       do i_beta=1,Nbeta,1
           call err_anal_complex(cicj_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
           write(60,'(4f15.8,I4)')i_beta*dt,dble(zmean),aimag(zmean),err,j
        enddo
        if(Nbeta.gt.0)write(60,*)
        if(Nbeta.gt.0)write(60,*)      
     end do
     do j=1,2*Nsite,1
       do i_beta=1,Nbeta,1
           call err_anal_complex(cicjh_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
           write(61,'(4f15.8,I4)')i_beta*dt,dble(zmean),aimag(zmean),err,j
        enddo
        if(Nbeta.gt.0)write(61,*)
        if(Nbeta.gt.0)write(61,*)
     end do
   endif
   if(Nbeta.gt.0)then
     if(I_twob.eq.1)then
       do j=1,Nbravais,1
         do i_beta=0,Nbeta,1
       
       !  call err_anal_complex(rho_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
       !  write(61,'(4f15.8,I4)')i_beta*dt,dble(zmean),aimag(zmean),err,j
           call err_anal_complex(nupnup_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
           write(62,'(4f15.8,I4)')i_beta*dt,dble(zmean),aimag(zmean),err,j
           call err_anal_complex(ndnnup_t_l(1:Nsamples,j,i_beta,i),Nsamples,zmean,err)
           write(63,'(4f15.8,I4)')i_beta*dt,dble(zmean),aimag(zmean),err,j
         enddo
         if(Nbeta.gt.0)write(62,*)
         if(Nbeta.gt.0)write(62,*)
         if(Nbeta.gt.0)write(63,*)
         if(Nbeta.gt.0)write(63,*)
       end do
     endif
   do i_beta=0,Nbeta,1
     call err_anal_complex(GreenP_t_l(1:Nsamples,i_beta,i),Nsamples,zmean,err)
     write(64,'(4f15.8)')i_beta*dt,dble(zmean),aimag(zmean),err
   enddo
   do i_beta=0,Nbeta,1
     call err_anal_complex(GreenH_t_l(1:Nsamples,i_beta,i),Nsamples,zmean,err)
     write(65,'(4f15.8)')i_beta*dt,dble(zmean),aimag(zmean),err
   enddo
   endif
end do
close(10)
close(20)
close(30)
close(40)
close(50)
end subroutine data_cpmc_rcpmc



!---------------------------
!Get the error bar of dat(N)
!---------------------------
subroutine err_anal_complex(dat,N,m,er)
implicit none
integer,intent(IN)::N
complex(kind=8),intent(IN)::dat(N)
complex(kind=8),intent(OUT)::m
real(kind=8),intent(OUT)::er
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
m=dcmplx(0.d0,0.d0)
do i=1,N,1
   m=m+dat(i)
end do
m=m/dble(N)

er=0.d0
do i=1,N,1
   er=er+dble(dat(i))**2
end do
er=er/dble(N)

er=er-dble(m)**2

if(ABS(er).GT.1.d-10) then
  if(er.lt.0.d0) then
    write(*,*) "Something is wrong in err_anal",er
  end if
else
  er=0.d0
end if

er=sqrt(er/dble(N-1))
end subroutine err_anal_complex






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
