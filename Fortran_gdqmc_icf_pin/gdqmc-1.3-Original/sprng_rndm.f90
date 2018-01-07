!--------------------------------------------------
!This random number use the sprng 2.0b to generate.
!--------------------------------------------------
MODULE rand_num
IMPLICIT NONE
#include "sprng_f.h"

 SPRNG_POINTER stream
 real(kind=8),parameter::pi_rn=2.d0*asin(1.d0)
 PRIVATE
 PUBLIC ::rndm,init_genrand,end_genrand,RN,gauss_rndm
CONTAINS

!------------------------------
!Initial the random number seed
!------------------------------
subroutine init_genrand

#ifdef MPI
#define USE_MPI
#include <mpif.h>
#endif


 integer streamnum, nstreams, seed
 integer myid, nprocs, ierror, junk
 integer gtype 

#ifdef MPI
 call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
 call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)
#else
 myid=0
 nprocs=1
#endif
 streamnum = myid
 nstreams = nprocs
 seed = 985456376
 !seed = make_sprng_seed()
 gtype=1
 if (myid .eq. 0) then
!#include "genf_types_menu.h"
     write(*,*) "genf_types:",gtype
     write(*,*) "seed:",seed
 endif

 stream = init_sprng(gtype,streamnum,nstreams,seed,SPRNG_DEFAULT)
 !Print the information of each stream
 !write(*,"('Process', i2, ', print information about stream:')") myid
 !junk = print_sprng(stream)
end subroutine init_genrand


!---------------------
!End the random number
!---------------------
subroutine end_genrand
integer::junk,ierror
 junk = free_sprng(stream)
end subroutine end_genrand 


!Generate a random number from 0~1
real(kind=8) function rndm()
 rndm=sprng(stream)
end function rndm

!Generate a random number from 1~nmax
integer function RN(nmax)
integer,intent(IN):: nmax
 RN=nmax*rndm()+1.d0; if(RN>nmax) RN=nmax
end function RN

!Generate a random number of gauss distribution
real(kind=8) function gauss_rndm()
real(kind=8)::tmp
 !tmp=1.d0-rndm()
 !if(tmp.LE.0.d0) tmp=tmp+1.d-15  !If we get rndm()=1.d0,tmp will be zero,
 !log(tmp) will be NA,use 1.d-15 to prevent this.
 !tmp=0.5d0*(rndm()+rndm()) !if u choose this, it will have bug,do not know why
 tmp=rndm()
 if(tmp.eq.0.d0) tmp=tmp+1.d-15
 gauss_rndm=sqrt(-2.d0*Log(tmp))*Cos(2.d0*pi_rn*rndm())
end function gauss_rndm

END MODULE rand_num
