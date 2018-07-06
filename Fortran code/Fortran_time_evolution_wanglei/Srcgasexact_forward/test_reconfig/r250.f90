!==============================================================================
module r250_params
implicit none
!----------------r250 parameter --------------------------------------------------------------------------
	integer(kind=4):: Kite
        integer(kind=4), dimension(0:255):: Nran
        integer(kind=4):: Iseed 
end module 
!==========================================================================



subroutine initialize_r250
use r250_params
implicit none

integer:: i,j,ici
real(kind=8),external::r250

!----------r250-initialization---------------------------------^M
!--- use a new iseed for a new sequence -----^M
!--- after running, 1< iseed < 2^(31). ------^M
!-- fixed initial ------------------^M
      iseed=1
!-- unfixed initial: ran from time, ----------^M
!--- like in T.Pang's book p48, use time. ----^M
!      iseed=time()^M
!--- iseed must use an odd number ------------^M
      if (iseed .le. 0) then
            print *, 'Wrong iseed! Stop.'
                stop
       endif
      if(mod(iseed,2)==0) iseed=iseed+1
!----------------------------------------------
    do j=0,255
           ici=0
           do i=1, 32
         ici=ISHFT(ici, 1)
         iseed=iseed*16807
         if (iseed .lt. 0) ici=ici+1
           enddo
           Nran(j)=ici
    enddo
        Kite=256        !required as initials.^M
!-------------------------------------------------^M
!     do i=1, 10
!      t=r250()
!     print *,myid, i, t
!     enddo
!=======================================================================^M

end subroutine

!========+=========+=========+=========+=========+=========+=========+=$
!    Program:  r250
!    TYPE   :  function
!    PURPOSE:  To generate random numbers in [0,1], using 
!              the random generator R250.
!    I/O    :
!    VERSION: 10-Apr-2007
!    COMMENT:  Note: IAND(255, x) is mod(x,256) operation.
!              See <<Compuational Physics>> by K.H. Hoffmann
!              and M. Schreiber, P4.
!========+=========+=========+=========+=========+=========+=========+=$
     function r250()
	 use r250_params
	 implicit none
	 integer:: i
	 real(kind=8):: r250
     i=IEOR(Nran(IAND(255,Kite-250)), Nran(IAND(255,Kite-103)))     
     Nran(IAND(255,Kite))=i
	 Kite=Kite+1
!-- 2^(-31)=4.656612873077393e-10 , random number in [0,1] ---
     r250=(1.0_8+4.656612873077393e-10_8*real(i,8))/2.0_8 
	 
	 end function r250
!===============================
