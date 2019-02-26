!This subroutine periodically adjust the ET 
!This subroutine is only used in therm and only for cpmc and rcpmc
!Otherwise we can set the ET from measurement.
!Keep in mind: If we are giving the code a very good try wave function,
!Be sure to give a very good ET.
subroutine adjustET(ad_ET)
use param
use adET_param
use model_param
use mc_loop_param
use mpi_serial_param
use project_param
use phi_param
use method_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
logical,intent(inout)::ad_ET
integer::i,j
real(kind=8)::tmp

if(max_crn.GE.0) then !CPMC or RCPMC

    if(.not.ad_ET) then !the first time
      do i=1,Nwalkers,1
         E_sum=E_sum+weight(i)
      end do

!DEBUG
!      write(*,*)
!      write(*,*)'In adjustET first option '
!      write(*,*)'ad_ET ',ad_ET
!      write(*,*)'weight ',weight
!      write(*,*)'E_sum = ',E_sum
!      write(*,*)

    else
      i_ad=i_ad+1

      if(PopContrlstep.eq.1)then
        do i=1,Nwalkers,1
          E_sum=E_sum+weight(i)
        end do
        if(i_ad.eq.1)then
          E_sum_old=E_sum
        endif
      endif 

#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(E_sum,tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      E_sum=tmp
#endif

!DEBUG
!if(i_ad.gt.1)then
!  write(*,*)
!  write(*,*)'ET BEFORE ADJUSTING ',ET
!  write(*,*)
!endif



      if(i_ad.EQ.1) then
        !The first time do not need to adjust
      else
        !adjust
        ET=ET-dlog(E_sum*m_w_old/E_sum_old)/(dble(PopContrlstep)*dt)
      end if

!DEBUG
!if(i_ad.gt.1)then
!      write(*,*)
!      write(*,*)'ET=ET-dlog(E_sum*m_w_old/E_sum_old)/(dble(PopContrlstep)*dt) '
!      write(*,*)
!      write(*,*)'In adjustET '
!      write(*,*)'ad_ET ',ad_ET
!      write(*,*)'weight ',weight
!      write(*,*)'E_sum = ',E_sum
!      write(*,*)'m_w_old, m_w = ',m_w_old,m_w
!      write(*,*)'E_sum_old = ',E_sum_old
!      write(*,*)'PopContrlstep = ',PopContrlstep
!      write(*,*)'dt = ',dt
!      write(*,*)'dble(PopContrlstep)*dt ',dble(PopContrlstep)*dt
!      write(*,*)'dlog(E_sum*m_w_old/E_sum_old)',dlog(E_sum*m_w_old/E_sum_old)
!      write(*,*)'dlog(E_sum*m_w_old/E_sum_old)/(dble(PopContrlstep)*dt) = ',dlog(E_sum*m_w_old/E_sum_old)/(dble(PopContrlstep)*dt) 
!      write(*,*)
!endif

      if(rank.eq.0) then
        write(*,*) "ADJUST ET CPMC OR RCPMC:",ET
        write(*,*) ""
        write(*,*) ""
        write(*,*) ""
      end if
      E_sum_old=E_sum
      E_sum=0.d0
      ad_ET=.false.
    end if!for ad_ET

else

  !We adjust ET in each step of measure when free-projection

end if !for max_crn
end subroutine adjustET


!Initial the ET
subroutine init_ET()
use lattice_param
use adET_param
use mpi_serial_param
use one_meas_param 
implicit none

 call measure()
 ET=dble(E_one) !Initial the ET
 etrial=ET
 if(rank.eq.0) then
    write(*,*) "INIT ET:",ET
    write(*,*) ""
    write(*,*) ""
    write(*,*) ""
 end if
end subroutine init_ET
