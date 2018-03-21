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

if(quando.eq.1)then
  if(rank.eq.0)then
    open(90,file='energy',status='unknown',access='append')
    rewind(90)
    write(90,'(4f15.5)')etrial,etrial,0.d0,0.d0
    flush(90)
  endif
endif
if(rank.eq.0)then
  write(90,'(4f15.5)')e_l(quando,0)/max_local,eE_l(quando,0)/max_local,kin_l(quando,0)/max_local,v_l(quando,0)/max_local
  flush(90)
endif

if(ipinn.eq.1)then
  iunit=91
  if(quando.eq.1)then
    if(rank.eq.0)then
      open(iunit,file='local_spin_density',status='unknown',access='append')
      rewind(iunit)
    endif  
  endif
  if(rank.eq.0)then
    do beta=1,2
      do alpha=1,2
        do sitei=1,Nsite
          x=dble(nofr_l(quando,sitei,alpha,beta,0))/max_local
          y=aimag(nofr_l(quando,sitei,alpha,beta,0))/max_local 
          write(iunit,'(3I4,3f15.8)')sitei,alpha,beta,x,y,0.d0
          flush(iunit)
        enddo
      enddo
    enddo
    write(iunit,*)
    flush(iunit)
    write(iunit,*)
    flush(iunit)
  endif
endif

if(ipinn.eq.0)then
  iunit=91
  if(quando.eq.1)then
    if(rank.eq.0)then
      open(iunit,file='spin_correlation',status='unknown',access='append')
      rewind(iunit)
    endif    
  endif
  if(rank.eq.0)then
    do sitei=1,Nbravais
      x=sisj_l(quando,sitei,0)/max_local
      write(iunit,'(I4,3f15.8)')sitei,0.d0,x,0.d0
      flush(iunit)
    enddo
    write(iunit,*)
    flush(iunit)
    write(iunit,*)
    flush(iunit)
  endif
endif


if(ipinn.eq.0)then
!  iunit=89
!  if(quando.eq.1)then
!    if(rank.eq.0)then
!!      open(iunit,file='cacb_correlation',status='unknown',access='append')
!      rewind(iunit)
!    endif
!  endif
!  if(rank.eq.0)then
!   do ib=1,Nbands
!     do jb=1,Nbands
!       do sitei=1,Nbravais
!         x=dble(cacb_l(quando,sitei,jb,ib,0))/max_local
!         y=aimag(cacb_l(quando,sitei,jb,ib,0))/max_local
!         write(iunit,'(3I4,2f15.8)')sitei,jb,ib,x,y
!         flush(iunit)
!       enddo
!     enddo 
!   enddo
!   write(iunit,*)
!   flush(iunit)
!   write(iunit,*)
!   flush(iunit)
!  endif
  
  iunit=88
  if(quando.eq.1)then
    if(rank.eq.0)then
      open(iunit,file='density_correlation',status='unknown',access='append')
      rewind(iunit)
    endif
  endif
  if(rank.eq.0)then
   do ib=1,Nbands
     do jb=1,Nbands
       do sitei=1,Nbravais
         x=ninj_l(quando,sitei,jb,ib,0)/max_local
         write(iunit,'(3I4,1f15.8)')sitei,jb,ib,x
         flush(iunit)
       enddo
     enddo
   enddo
   write(iunit,*)
   flush(iunit)
   write(iunit,*)
   flush(iunit)
  endif

  iunit=101
  if(quando.eq.1)then
    if(rank.eq.0)then
      open(iunit,file='szsz_correlation',status='unknown',access='append')
      rewind(iunit)
    endif
  endif
  if(rank.eq.0)then
   do ib=1,Nbands
     do jb=1,Nbands
       do sitei=1,Nbravais
         x=szsz_l(quando,sitei,jb,ib,0)/max_local
         write(iunit,'(3I4,1f15.8)')sitei,jb,ib,x
         flush(iunit)
       enddo
     enddo
   enddo
   write(iunit,*)
   flush(iunit)
   write(iunit,*)
   flush(iunit)
  endif

  iunit=87
  if(quando.eq.1)then
    if(rank.eq.0)then
      open(iunit,file='ZetaN',status='unknown',access='append')
      rewind(iunit)
    endif
  endif
  if(rank.eq.0)then
   do ib=1,Dimen
         x=dble(ZetaN_l(quando,ib,0))/max_local
         y=aimag(ZetaN_l(quando,ib,0))/max_local
         write(iunit,'(1I2,2f15.8)')ib,x,y
         flush(iunit)
   enddo
   write(iunit,*)
   flush(iunit)
   write(iunit,*)
   flush(iunit)
  endif

endif



if(Npair.eq.0)go to 1

!write(*,*)'print_observables ',quando,rank
if(quando.eq.1)then
  if(rank.eq.0)then
    do ipair=1,Npair
      call xifs(sfix,ipair)
      filename='pairing.'//sfix
      iunit=91+ipair
      open(iunit,file=filename,status='unknown',access='append')       
      rewind(iunit)
    enddo
  endif
endif
if(rank.eq.0)then
  do ipair=1,Npair
    iunit=91+ipair
    do sitei=1,Nbravais,1
      write(iunit,*)sitei,dble(didj_l(quando,sitei,ipair,0))/max_local,aimag(didj_l(quando,sitei,ipair,0))/max_local
      flush(iunit)
    enddo
  enddo
endif

1 if(Nbeta.eq.0)return

!iunit=94
!if(quando.eq.1)then
!  if(rank.eq.0)then
!    open(iunit,file='green_p_dynamical',status='unknown',access='append')
!    rewind(iunit)
!  endif
!endif
!if(rank.eq.0)then
!    do i_beta=1,Nbeta,1
!      do sitei=1,2*Nsite,1
!        x=dble(cicj_t_l(quando,sitei,i_beta,0))/max_local
!        y=aimag(cicj_t_l(quando,sitei,i_beta,0))/max_local
!        write(iunit,*)sitei,i_beta,x,y
!      enddo
!    enddo
!endif
!
!iunit=89
!if(quando.eq.1)then
!  if(rank.eq.0)then
!    open(iunit,file='green_h_dynamical',status='unknown',access='append')
!    rewind(iunit)
!  endif
!endif
!if(rank.eq.0)then
!    do i_beta=1,Nbeta,1
!      do sitei=1,2*Nsite,1
!        x=dble(cicjh_t_l(quando,sitei,i_beta,0))/max_local
!        y=aimag(cicjh_t_l(quando,sitei,i_beta,0))/max_local
!        write(iunit,*)sitei,i_beta,x,y
!      enddo
!    enddo
!endif

if(I_twob.eq.1)then

iunit=95
if(quando.eq.1)then
  if(rank.eq.0)then
    open(iunit,file='nupnup_dynamical',status='unknown',access='append')
    rewind(iunit)
  endif
endif
if(rank.eq.0)then
    do i_beta=0,Nbeta,1
      do sitei=1,Nbravais,1
        x=dble(nupnup_t_l(quando,sitei,i_beta,0))/max_local
        y=aimag(nupnup_t_l(quando,sitei,i_beta,0))/max_local
        write(iunit,*)sitei,i_beta,x,y
        flush(iunit)
      enddo
    enddo
endif

iunit=96
if(quando.eq.1)then
  if(rank.eq.0)then
    open(iunit,file='ndnnup_dynamical',status='unknown',access='append')
    rewind(iunit)
  endif
endif
if(rank.eq.0)then
    do i_beta=0,Nbeta,1
      do sitei=1,Nbravais,1
        x=dble(ndnnup_t_l(quando,sitei,i_beta,0))/max_local
        y=aimag(ndnnup_t_l(quando,sitei,i_beta,0))/max_local
        write(iunit,*)sitei,i_beta,x,y
        flush(iunit)
      enddo
    enddo
endif

endif

!iunit=97
!if(quando.eq.1)then
!  if(rank.eq.0)then
!    open(iunit,file='nup_k_nup_-k_dynamical',status='unknown',access='append')
!    rewind(iunit)
!  endif
!endif
!if(rank.eq.0)then
!    do i_beta=0,Nbeta,1
!      do k=1,Nbravais,1
!        do n=1,Dimen,1
!         kl(n)=coor(k,n)-1
!        end do
!
!        !fourier of the momentum k
!        kk=zero
!        do m=1,Nbravais,1
!          ph=0.d0
!          do n=1,Dimen,1
!            ph=ph+2.d0*Pi*dble(kl(n)*(coor(m,n)-1))/dble(Nl(n))
!          enddo
!          kk=kk+exp(Xi*ph)*nupnup_t_l(quando,m,i_beta,0)/max_local
!        enddo
!        x=dble(kk)
!        y=aimag(kk)
!        write(iunit,*)k,i_beta,x,y
!        flush(iunit)
!      enddo
!    enddo
!endif

!iunit=98
!if(quando.eq.1)then
!  if(rank.eq.0)then
!    open(iunit,file='ndn_k_nup_-k_dynamical',status='unknown',access='append')
!    rewind(iunit)
!  endif
!endif
!if(rank.eq.0)then
!    do i_beta=0,Nbeta,1
!      do k=1,Nbravais,1
!        do n=1,Dimen,1
!         kl(n)=coor(k,n)-1
!        end do
!
!        !fourier of the momentum k
!        kk=zero
!        do m=1,Nbravais,1
!          ph=0.d0
!          do n=1,Dimen,1
!            ph=ph+2.d0*Pi*dble(kl(n)*(coor(m,n)-1))/dble(Nl(n))
!          enddo
!          kk=kk+exp(Xi*ph)*ndnnup_t_l(quando,m,i_beta,0)/max_local
!        enddo
!        x=dble(kk)
!        y=aimag(kk)
!        write(iunit,*)k,i_beta,x,y
!        flush(iunit)
!      enddo
!    enddo
!endif


iunit=99

if(quando.eq.1)then
  if(rank.eq.0)then
    open(iunit,file='cnuc+mu_dynamical',status='unknown',access='append')
    rewind(iunit)
  endif
endif
if(rank.eq.0)then
    do i_beta=0,Nbeta,1
        x=dble(GreenP_t_l(quando,i_beta,0))/max_local
        y=aimag(GreenP_t_l(quando,i_beta,0))/max_local
        write(iunit,*)1,i_beta,x,y
        flush(iunit)
    enddo
endif

iunit=100

if(quando.eq.1)then
  if(rank.eq.0)then
    open(iunit,file='c+nucmu_dynamical',status='unknown',access='append')
    rewind(iunit)
  endif
endif
if(rank.eq.0)then
    do i_beta=0,Nbeta,1
        x=dble(GreenH_t_l(quando,i_beta,0))/max_local
        y=aimag(GreenH_t_l(quando,i_beta,0))/max_local
        write(iunit,*)1,i_beta,x,y
        flush(iunit)
    enddo
endif



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
