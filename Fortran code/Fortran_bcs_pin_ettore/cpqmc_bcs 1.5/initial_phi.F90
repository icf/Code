!-----------------------------------------
!Decide what kind of phi initial methods.
!-----------------------------------------
subroutine get_phi()
use param
use lattice_param
use model_param
use phiT_param
use phi_param
use method_param
use mc_loop_param
use mpi_serial_param
implicit none
integer::Nw
complex(kind=8),allocatable::phi_init(:,:,:)
complex(kind=8),allocatable::phi_coe(:)
real(kind=8),allocatable::prob_init(:)
complex(kind=8),allocatable::r_init(:) !The complex phase

integer::table(Nwalkers)
integer::i,j,k
real(kind=8)::norm
real(kind=8),external::dasum
complex(kind=8)::tot_tmp

!allocated the arrays
call deallocate_phi()
call allocate_phi()

!initial the weight dlogw
call initial_w_imp()

if(I_wavefun.eq.2)then
  PP=1
  if(rank.eq.0)then
    write(*,*)
    write(*,*)'BCS option sets PP=1, build phi from phT'
    write(*,*)
  endif
endif

!Get Nw
if(PP.EQ.0) then
  !read Nw 
  open(unit=10,file='phi_N.dat',status='old')
  read(10,*) Nw
  close(10)
else if (PP.EQ.1) then
  Nw=Dtot
else
  write(*,*) "Something is wrong with PP input"
  call mystop
end if

if(rank.eq.0) then
  write(*,*) "Nw of the phi:",Nw
end if


!Do not do population control, read to weight directly.
!MNwalkers = total number of walkers (includes number of cores)
!Nw = (if written from file) number of Slater Determinants in initial wf
#ifdef MPI
 if(Nw.EQ.MNwalkers) then
#else
 if(Nw.EQ.Nwalkers) then
#endif

     if(PP.EQ.0) then
       !Read phi from file
       !Read coe for weight and rx also dlogw
       call read_phi_rank()
       if(rank.eq.0) write(*,*) "read phi from file directly."
     else
       write(*,*) "Nw.EQ.(M)Nwalkers: we should read from the file."
       call mystop
     end if


!Do population control, read to prob_init then control it do different walkers.
!MNwalkers = total number of walkers (includes number of cores)
!Nw = (if written from file) number of Slater Determinants in initial wf
#ifdef MPI
 else if(Nw.LT.MNwalkers) then
#else
 else if(Nw.LT.Nwalkers) then
#endif

   !Get phi_init,prob_init and r_init
   allocate(phi_init(2*Nsite,Ntot,Nw),phi_coe(Nw),prob_init(Nw),r_init(Nw))
   if(PP.EQ.0) then
     !read phi_init
     !read phi_coe
     call read_phi(Nw,Nsite,phi_coe,phi_init)
     if(rank.eq.0) write(*,*) "read phi from file and pop contr."
   else if (PP.EQ.1) then !Get from phiT
     if(rank.eq.0) write(*,*) "read phi from phiT and pop contr."
     call zcopy(2*Nsite*Ntot*Nw,phiT(1,1,1),1,phi_init(1,1,1),1)
     call zcopy(Nw,coe_multi(1),1,phi_coe(1),1)
   else
     write(*,*) "Something is wrong with PP input"
     call mystop
   end if

   do k=1,Nw,1
      if(crn.GT.0.d0) then !fp
        prob_init(k)=abs(phi_coe(k))
        r_init(k)=phi_coe(k)/abs(phi_coe(k)) !The first rx(i)
      else   !CPMC
        call get_tot_tmp(tot_tmp,phi_init(1,1,k))
        prob_init(k)=abs(phi_coe(k)*tot_tmp)
        r_init(k)=phi_coe(k)*tot_tmp/prob_init(k)
        !if(rank.eq.0) write(*,*) r_init(k)
      end if

    end do


 



   norm=1.d0/dasum(Nw,prob_init,1)
   call dscal(Nw,norm,prob_init,1)


   !--------------------------------------------
   !Get phi from phi_init according to prob_init
   !and give out the table
   !--------------------------------------------
   call initial_phi(Nw,prob_init,phi_init,table)


   !------
   !Get rx
   !------
   do i=1,Nwalkers,1
      rx(i)=r_init(table(i))
   end do




   if(allocated(phi_init)) deallocate(phi_init)
   if(allocated(phi_coe)) deallocate(phi_coe)
   if(allocated(prob_init)) deallocate(prob_init)
   if(allocated(r_init)) deallocate(r_init)
 else
   if(rank.eq.0) then 
     write(*,*) "Nw can not be larger than number of walkers",Nw
     call mystop
   end if
 end if


if(crn.LT.0)  call cal_imp_ovlap()


!if(rank.eq.0) then
!  write(*,*)
!  write(*,*)
!  write(*,*)'Initial importance function ',PhiT
!  write(*,*)
!  write(*,*)
!endif

end subroutine get_phi



!-----------------------------
!read the phi_coe and phi_init
!-----------------------------
subroutine read_phi(Nw,Nsite,phi_coe,phi_init)
use param
use io_module
use model_param
implicit none
integer,intent(IN)::Nw
integer,intent(IN)::Nsite
complex(kind=8),intent(OUT)::phi_coe(Nw)
complex(kind=8),intent(OUT)::phi_init(2*Nsite,Ntot,Nw)
integer::i,j,k,j1,j2,npartu,npartd
real(kind=8)::eps
complex(kind=8)::dummy,dummy2
complex(kind=8),dimension(:,:,:),allocatable::phi_read
character(len=300)::filen

eps=1.d-3

!Get phi_coe:
open(unit=10,file='phi_coe.dat',status='old')
  do i=1,Nw,1
     read(10,*) phi_coe(i)
  end do
close(10)
!Read phi_init
allocate(phi_read(2*Nsite,Ntot,Nw))
open(unit=10,file='phi.dat',status='old')
do i=1,Nw,1
  do k=1,Ntot,1
     do j=1,2*Nsite,1
        read(10,*) phi_read(j,k,i)
     end do
  end do
end do
close(10)


if(dtype.eq.'c')then
  phi_init=phi_read
elseif(dtype.eq.'d'.or.dtype.eq.'m')then
  do k=1,Nw,1
    npartu=0
    npartd=0
    do j2=1,Ntot
      dummy=zero
      do j1=1,Nsite,1
        dummy=dummy+conjg(phi_read(j1,j2,k))*phi_read(j1,j2,k)
      enddo
      if(abs(dummy-one).lt.eps)then
        npartu=npartu+1
        do j1=1,Nsite,1
           phi_init(j1,npartu,k)=phi_read(j1,j2,k)
        enddo
        do j1=Nsite+1,2*Nsite,1
           phi_init(j1,npartu,k)=zero
        enddo
      else
        dummy2=zero
        do j1=Nsite+1,2*Nsite,1
          dummy2=dummy2+conjg(phi_read(j1,j2,k))*phi_read(j1,j2,k)
        enddo
        if(abs(dummy2-one).lt.eps)then
          npartd=npartd+1
          do j1=1,Nsite,1
            phi_init(j1,Nspin(1)+npartd,k)=zero
          enddo
          do j1=Nsite+1,2*Nsite,1
            phi_init(j1,Nspin(1)+npartd,k)=phi_read(j1,j2,k)
          enddo
        endif
      endif
    enddo
  enddo
endif

deallocate(phi_read)

!call createFileName(filen,'./')
!call appendBaseName(filen,'l',Nsite)
!call appendBaseName(filen,'nt',Ntot)
!call appendBaseName(filen,'basis.dat')
!call openUnit(filen,10,'B')
!  do i=1,Nw,1
!    read(10) phi_init(1:2*Nsite,1:Ntot,i)
!  end do
!close(10)
end subroutine read_phi



!---------------------------------------------------
!write the coe and phi to the file according to rank
!---------------------------------------------------
subroutine write_phi_rank()
use param
use mpi_serial_param
use io_module
use mc_loop_param
use phi_param
use lattice_param
use model_param
use method_param
use project_param
implicit none
character(len=300)::phi_name,coe_name
complex(kind=8)::phtmp(2*Nsite,Ntot)
integer::i,j,k
call createFileName(phi_name,'phi')
call appendBaseName(phi_name,'_',rank)
call appendBaseName(phi_name,'.dat')
call createFileName(coe_name,'coe')
call appendBaseName(coe_name,'_',rank)
call appendBaseName(coe_name,'.dat')
!call openUnit(phi_name,10,'R')
!call openUnit(coe_name,20,'R')
call openUnit(phi_name,10,'R')
call openUnit(coe_name,20,'R')
do i=1,Nwalkers,1

   !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi(1,1,i),2*Nsite,zero,phtmp,2*Nsite)
   call k_to_ph_dc(exp_halfK,phi(1,1,i),phtmp(1,1))

   do k=1,Ntot,1
      do j=1,2*Nsite,1
         write(10,*) phtmp(j,k)
      end do
   end do
   !write(10) phtmp(1:2*Nsite,1:Ntot)

   if(crn.GT.0.d0) then !fp
     write(20,*) rx(i)*weight(i)
   else
     write(20,*) rx(i)*weight(i)/tot_imp(i)
   end if

end do
close(10)
close(20)
end subroutine write_phi_rank



!--------------------------------------------------
!read the coe and phi to the file according to rank
!--------------------------------------------------
subroutine read_phi_rank()
use mpi_serial_param
use io_module
use mc_loop_param
use phi_param
use lattice_param
use model_param
use method_param
implicit none
character(len=300)::phi_name,coe_name
complex(kind=8)::tmp,tot_tmp
integer::i,j,k

call createFileName(phi_name,'phi')
call appendBaseName(phi_name,'_',rank)
call appendBaseName(phi_name,'.dat')
call createFileName(coe_name,'coe')
call appendBaseName(coe_name,'_',rank)
call appendBaseName(coe_name,'.dat')
!call openUnit(phi_name,10,'O')
!call openUnit(coe_name,20,'O')
call openUnit(phi_name,10,'B')
call openUnit(coe_name,20,'B')
do i=1,Nwalkers,1

   !do k=1,Ntot,1
   !   do j=1,2*Nsite,1
   !      read(10,*) phi(j,k,i)
   !   end do
   !end do
   read(10) phi(1:2*Nsite,1:Ntot,i)

   read(20) tmp
   if(crn.GT.0.d0) then !fp
     weight(i)=abs(tmp)
     rx(i)=tmp/weight(i)
   else  !CPMC
     call get_tot_tmp(tot_tmp,phi(1,1,i))
     weight(i)=abs(tmp*tot_tmp)
     rx(i)=tmp*tot_tmp/weight(i)
     !write(*,*) rx(i)
   end if

   !set dlogw
   if(weight(i).GT.0.d0) then
     dlogw(i)=dlog(weight(i))
   else if(weight(i).EQ.0.d0) then
     dlogw(i)=-1d100
   else
     write(*,*) "Something is wrong with weight init.",weight(i),i
     call mystop
   end if

end do
close(10)
close(20)
end subroutine read_phi_rank





!-----------------------
!We initial weight dlogw
!-----------------------
subroutine initial_w_imp
use param
use phi_param
implicit none
weight=1.d0   
dlogw=0.d0   !All the dlogw here is dlog(weight)
end subroutine initial_w_imp




!---------------------------------------------------------------
!We initial phi and inverse overlap function in this subroutine.
!---------------------------------------------------------------
subroutine initial_phi(Nw,prob_init,phi_init,table)
use param
use model_param
use lattice_param
use phi_param
use mc_loop_param
use rand_num
use method_param
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer,intent(IN)::Nw
real(kind=8),intent(IN)::prob_init(Nw)
complex(kind=8),intent(IN)::phi_init(2*Nsite,Ntot,Nw)
integer,intent(OUT)::table(Nwalkers)
integer,allocatable::Mtable(:)
real(kind=8)::eta!Get the probability of initial wave function
integer::i,j,k,npartu,npartd
integer::j1,j2
real(kind=8)::eps
complex(kind=8)::dummy,dummy2

eps=1.d-3

#ifdef MPI
 if(rank.eq.0) allocate(Mtable(MNwalkers))
#endif

#ifdef MPI
   if(rank.eq.0) call distri_p(Nw,prob_init,MNwalkers,Mtable)
#else
   call distri_p(Nw,prob_init,Nwalkers,table)
#endif


#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_SCATTER(Mtable(1),Nwalkers,MPI_INTEGER,table(1),Nwalkers,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#endif

 do i=1,Nwalkers,1
   !---------------------------------------------------------
   !Set the initial wavefunction proportional to prob_init(k)
   !---------------------------------------------------------
   eta=1d-3
   !eta=100.d0 !test
   if(Nw.eq.1) eta=0.d0

   !Only use decoupled term here, for we do not want mixed term here.
   !Be careful here!
   if(dtype.EQ.'c') then
   
     do j1=1,2*Nsite,1
        do j2=1,Ntot,1
           phi(j1,j2,i)=phi_init(j1,j2,table(i))+eta*(2.d0*rndm()-1.d0)
        end do
     end do
   
   else if(dtype.EQ.'d'.or.dtype.EQ.'m') then

     !Just put in the up and dn term
     do j1=1,Nsite,1
        do j2=1,Nspin(1),1
           phi(j1,j2,i)=phi_init(j1,j2,table(i))+eta*(2.d0*rndm()-1.d0)
        end do
        !test
        do j2=Nspin(1)+1,Ntot,1
        !   write(*,*) rndm()
           phi(j1,j2,i)=zero
        end do
        !end test
     end do
     
     
     do j1=Nsite+1,2*Nsite,1
        !test
        do j2=1,Nspin(1),1
        !   write(*,*) rndm()
           phi(j1,j2,i)=zero
        end do
        !end test
        do j2=Nspin(1)+1,Ntot,1
           phi(j1,j2,i)=phi_init(j1,j2,table(i))+eta*(2.d0*rndm()-1.d0)
        end do
     end do

   end if
 end do


if(rank.eq.0)then
open(2,file='psi_init.info',status='unknown')
do i=1,1,1
  do k=1,Ntot,1
     do j=1,2*Nsite,1
        write(2,*) phi(j,k,i)
     end do
  end do
end do
close(2)
endif

!stop 'DEBUG'

#ifdef MPI
 if(rank.eq.0) deallocate(Mtable)
#endif

end subroutine initial_phi






!----------------------------------------------------------------
!This subroutine generate n walkers by the distribution of p(1:m)
!Similar to the subroutine reconfiguration, while here m/=n+++++
!----------------------------------------------------------------
subroutine distri_p(m,p,n,table)
use rand_num
implicit none
integer,intent(IN)::m,n
real(kind=8),intent(IN)::p(m)
integer,intent(OUT)::table(n)
real(kind=8)::wsum,z
integer::i,j
wsum=0.d0
do i=1,m,1
   wsum=wsum+p(i)
end do
if(abs(wsum-1.d0).GT.1d-10) then
  write(*,*) "p input is not normalized!"
  call mystop
end if

j=1
wsum=p(1)
do i=1,n,1
   z=dble((i-1)+rndm())/dble(n)
   do while(z>wsum)
      j=j+1
      wsum=wsum+p(j)
   enddo
   table(i)=j
enddo

end subroutine distri_p


!-------------------------------
!Get the tot_tmp=<phiT|phi_init>
!-------------------------------
subroutine get_tot_tmp(tot_tmp,phi_init)
use param
use lattice_param
use model_param
use phiT_param
implicit none
complex(kind=8),intent(OUT)::tot_tmp
complex(kind=8),intent(IN)::phi_init(2*Nsite,Ntot)
complex(kind=8)::imp
integer::m,i,k,j
tot_tmp=zero
do m=1,Dtot,1
   !call deter_overlap_imp(2*Nsite,Ntot,phiT(1,1,m),phi_init(1,1),imp)
   if(I_wavefun.eq.1)then
     call imp_fun_dc(phiT(1,1,m),phi_init(1,1),imp)
   elseif(I_wavefun.eq.2)then
     if(Dtot.ne.1)then
       write(*,*)'Something is wrong in BCS case of get_tot_tmp '
       stop
     endif
     if(Nzeta.eq.0)then
       call bcs_imp_fun_dc(FPairing(1,1),phi_init(1,1),imp)
     else
       call unp_bcs_imp_fun_dc(FPairing(1,1),DUnpaired(1,1),phi_init(1,1),imp)
     endif
   endif
   !if(abs(imp).LT.1d-15) then
   !   write(*,*) "Small overlap phiT with phi_init",m
   !  call mystop
   !end if
   tot_tmp=tot_tmp+conjg(coe_multi(m))*imp
end do

end subroutine get_tot_tmp



