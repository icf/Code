!--------------------------------------------------
!Initial rank and Nsize, if mpi then inital the MPI
!--------------------------------------------------
subroutine start_code()
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
#ifdef MPI
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,Nsize,ierr)
#else
 rank=0
 Nsize=1
#endif
end subroutine start_code


!----------------------------------------------------
!This subroutine is used to read parameter from param
!at the beginning of running the code++++++++++++++++
!----------------------------------------------------
subroutine readparam
use lattice_param
use model_param
use mpi_serial_param
use xhf_param
implicit none
open(unit=10,file='param',status='old')

!read the lattice parameter
read(10,*) set
read(10,*) Nsite
read(10,*) Nhop
read(10,*) Dimen
read(10,*) Nl(1)
read(10,*) Nl(2)
read(10,*) Nl(3)
read(10,*) kbound(1)
read(10,*) kbound(2)
read(10,*) kbound(3)

!read the model parameter
read(10,*) t1
read(10,*) onsitU
read(10,*) Ntot
read(10,*) dtype
read(10,*) Nspin(1)
read(10,*) Nspin(2)

if(dtype.EQ.'d') then
  Ntot=Nspin(1)+Nspin(2)
end if

!read for the x-hf
read(10,*) dmin
read(10,*) dmax
read(10,*) dstep

close(10)
end subroutine readparam


!------------------------------------------------
!This subroutine initial the mpi htype and phtype
!------------------------------------------------
subroutine init_mpi_type()
use model_param
use mpi_serial_param
use lattice_param
implicit none
integer::bl(2*Nsite)
integer::disp(2*Nsite)
integer::i
#ifdef MPI
include "mpif.h"
#endif


#ifdef MPI
 if(dtype.EQ.'c') then
   call MPI_TYPE_CONTIGUOUS(4*Nsite*Nsite,MPI_DOUBLE_COMPLEX,htype,ierr)
   call MPI_TYPE_CONTIGUOUS(2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,phtype,ierr)
 else if(dtype.EQ.'d') then

   bl=Nsite
   do i=1,Nsite,1
      disp(i)=(i-1)*2*Nsite
   end do
   do i=1,Nsite,1
      disp(i+Nsite)=(Nsite+i-1)*2*Nsite+Nsite
   end do
   call MPI_TYPE_INDEXED(2*Nsite,bl,disp,MPI_DOUBLE_COMPLEX,htype,ierr)


   bl=Nsite
   do i=1,Nspin(1),1
      disp(i)=(i-1)*2*Nsite
   end do
   do i=1,Nspin(2),1
      disp(i+Nspin(1))=(Nspin(1)+i-1)*2*Nsite+Nsite
   end do
   call MPI_TYPE_INDEXED(Ntot,bl,disp,MPI_DOUBLE_COMPLEX,phtype,ierr)

 end if
 call MPI_TYPE_COMMIT(htype,ierr)
 call MPI_TYPE_COMMIT(phtype,ierr)
#else
return
#endif

end subroutine init_mpi_type


!-----------------------------------------------------------------
!This subroutine allocate the arrays we need to use in set lattice
!-----------------------------------------------------------------
subroutine allocate_lattice_array()
use lattice_param
implicit none
allocate(hopt(Nhop))
allocate(sit(Nhop,2))
allocate(coor(Nsite,Dimen))
allocate(Hzero(2*Nsite,2*Nsite))
allocate(Tmatrix(Nsite,Dimen))
end subroutine allocate_lattice_array


!-------------------------------------------------------------------
!This subroutine deallocate the arrays we need to use in set lattice
!-------------------------------------------------------------------
subroutine deallocate_lattice_array()
use lattice_param
implicit none
if(allocated(hopt)) deallocate(hopt)
if(allocated(sit)) deallocate(sit)
if(allocated(coor)) deallocate(coor)
if(allocated(Hzero)) deallocate(Hzero)
if(allocated(Tmatrix)) deallocate(Tmatrix)
end subroutine deallocate_lattice_array

!------------------------------
!allocate the array used in xhf
!------------------------------
subroutine allocate_xhf_array()
use xhf_param
use lattice_param
use model_param
implicit none
allocate(h0xhf(2*Nsite,2*Nsite))
allocate(ph(2*Nsite,Ntot))
allocate(ph_save(2*Nsite,Ntot))
end subroutine allocate_xhf_array


!--------------------------------
!deallocate the array used in xhf
!--------------------------------
subroutine deallocate_xhf_array()
use xhf_param
implicit none
if(allocated(h0xhf)) deallocate(h0xhf)
if(allocated(ph)) deallocate(ph)
if(allocated(ph_save)) deallocate(ph_save)
end subroutine deallocate_xhf_array

!-----------------------------------------------------------------
!This subroutine deallocate all the arrays at the end of the code.
!-----------------------------------------------------------------
subroutine deallocatearray()
implicit none
call deallocate_lattice_array()
call deallocate_xhf_array()
end subroutine deallocatearray



!-------------------------------------------------
!This subroutine clean call the things before stop
!-------------------------------------------------
subroutine clean
use timing_module
use rand_num
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
call end_genrand()
call deallocatearray()
call EndTiming()
#ifdef MPI
call MPI_TYPE_FREE(htype,ierr)
call MPI_TYPE_FREE(phtype,ierr)
call MPI_Finalize(ierr)
#endif
end subroutine clean



!----------------------------------------------------------------
!This subroutine end the program before deallocate all the arrays
!----------------------------------------------------------------
subroutine mystop
implicit none
call clean()
stop
end subroutine mystop
