

!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: buildbasis 
! NOTICE :
! TYPE   : subroutine
! PURPOSE: produce Table(nup, ndo, num) which stores the decimal 
!          representation of the num_th basis in (nup, ndo) block ;
!          and also produce the invTable: sequential number of the  
!          decimal rep of a configuration , the sequential number is 
!          not unique for each subspace
! I/O    : nup, ndo/ Table, invTable
! VERSION: 30-Dec-2007
! AUTHOR : Lei Wang, IOP,CAS
! COMMENT: 0 <= nup <= Nsite: the number of up spin electrons (input)
!          0 <= ndo <= Nsite: the number of down spin electrons (input)
!          1 <= num <=Nstate(nup,ndo): an integer number (input)
!          Table : the decimal integer representing the num_th state.
!          A  basis state is represented by 2*Nsite digits in this way:
!          | 1up,1do,2up,2do,...,Nsup,Nsdo >, where 1up,1do...=0,1.
!          It is uniquely recorded by an integer:
!          vbas=1up*(2**0)+1do*(2**1)+2up*(2**2)+2do*(2**3)+...
!          vbas->Table(nup, ndo, num), where num is 1,2,...,Nstate(nup,ndo).\
!========+=========+=========+=========+=========+=========+=========+=$
     subroutine buildbasis(nup, ndo)
     use hubbard_param
     use flags_param
     implicit none

	 integer,intent(IN):: nup, ndo 
     
	 integer,dimension(Nsite):: vup,vdo
     integer,dimension(2*Nsite)::vtot
     
	 integer::flagdo, flagup,i,kk
     integer(kind=4)::count


     if (2*Nsite > 28) then
       print *, 'The number of digits of Table may be insufficient!'
       stop
     endif

     Table=0
    
     count=1
!----- initializing the config of eup and edo.---
!--- initial config is |1,1,..,1, 0,0,...,0> ---
 	 do i=1, nup
  	   vup(i)=1
 	 enddo
	 do i=nup+1, Nsite
	   vup(i)=0
 	 enddo

10	 do i=1, ndo
	   vdo(i)=1
 	 enddo
	 do i=ndo+1, Nsite
	   vdo(i)=0
 	 enddo
      
     do i=1, Nsite
        vtot(2*i-1)=vup(i)
        vtot(2*i)=vdo(i)
     enddo
!--- 1st state counts=1 ------
     do i=1,2*Nsite
        Table(count)=Table(count)+vtot(i)*(2**(i-1))
	 enddo	
		if(invTable_flag==1) invTable(Table(count))=count
!---------------------------------
20   call next(vdo,flagdo)        !a new config for down spin elec.

       if (flagdo .eq.  1) then    !new config obtained successfully.
         count=count+1
         do i=1, Nsite
           vtot(2*i-1)=vup(i)
           vtot(2*i)=vdo(i)
         enddo
         
         do i=1,2*Nsite
           Table(count)=Table(count)+vtot(i)*(2**(i-1))
	     enddo  
		  if (invTable_flag==1) invTable(Table(count))=count
        
         goto 20
       else if (flagdo .eq. 2) then   !error occured.
         write (*,*) 'error in "buildbasis"!'
	     stop
      endif                          !no more new config for down elec.
!-------------
      call next(vup,flagup)          !a new config for up spin elec.
      if (flagup .eq. 1) then
        count=count+1	
        goto 10                      !update up elec config, go back to scan down elec.
      else if (flagup .eq. 2) then
        write (*,*) 'error in "buildbasis"!'
        stop
      endif
!--- scan over, check ----
      if (count .ne. Nstate(nup,ndo)) then
        write (*,*) 'nstates is inconsistent with counting.'
        stop
      endif
      if (count .gt. Nsmax) then
        print *, 'Nsmax is smaller than the maxmum Nstate, wrong;!'
        stop
      endif
!      write (*,*) 'nup,ndo=',nup,' ',ndo,' count=',count
!      write (*,*) Nstate(nup,ndo)
!---------------------------
	 
    end subroutine buildbasis
!=======================================================================

!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: ini_code 
! NOTICE : 
!
! TYPE   : subroutine
! PURPOSE: initialize a code , for an given pair of (nup, ndo)
! I/O    :
! VERSION: 30-Dec-2007
! AUTHOR : Lei Wang, IOP, CAS
! COMMENT: 0 <= nup <= Nsite: the number of up spin electrons (input)
!          0 <= ndo <= Nsite: the number of down spin electrons (input)
!          1 <= num <=Nstate(nup,ndo): an integer number (input)
!          Table : the decimal integer representing the num_th state.
!          A  basis state is represented by 2*Nsite digits in this way:
!          | 1up,1do,2up,2do,...,Nsup,Nsdo >, where 1up,1do...=0,1.
!          It is uniquely recorded by an integer:
!          vbas=1up*(2**0)+1do*(2**1)+2up*(2**2)+2do*(2**3)+...
!          vbas->Table(nup, ndo, num), where num is 1,2,...,Nstate(nup,ndo).
!========+=========+=========+=========+=========+=========+=========+=$
     subroutine ini_code(nup, ndo, vtot)
     use hubbard_param
     implicit none

	 integer,intent(IN):: nup, ndo 
     integer,dimension(Nsite):: vup,vdo
     integer,dimension(2*Nsite)::vtot
     integer::flagdo, flagup,i,kk
     integer(kind=4)::count


     if (2*Nsite .ge. 28) then
       print *, 'The number of digits of Table may be insufficient!'
       stop
     endif

   
!----- initializing the config of eup and edo.---
!--- initial config is |1,1,..,1, 0,0,...,0> ---
 	 do i=1, nup
  	   vup(i)=1
 	 enddo
	 do i=nup+1, Nsite
	   vup(i)=0
 	 enddo

	 do i=1, ndo
	   vdo(i)=1
 	 enddo
	 do i=ndo+1, Nsite
	   vdo(i)=0
 	 enddo
      
     do i=1, Nsite
        vtot(2*i-1)=vup(i)
        vtot(2*i)=vdo(i)
     enddo


    end subroutine ini_code
!=======================================================================


!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: codes_shift
! NOTICE :
!
! TYPE   : subroutine
! PURPOSE: shift code code(2*Nsite)---> voto(2*Nsite)
! I/O    :
! VERSION: 28-Dec-2007
! AUTHOR : Lei Wang, IOP, CAS
! COMMENT: 0 <= nup <= Nsite: the number of up spin electrons (input)
!          0 <= ndo <= Nsite: the number of down spin electrons (input)
!          1 <= num <=Nstate(nup,ndo): an integer number (input)
!          Table : the decimal integer representing the num_th state.
!          A  basis state is represented by 2*Nsite digits in this way:
!          | 1up,1do,2up,2do,...,Nsup,Nsdo >, where 1up,1do...=0,1.
!          It is uniquely recorded by an integer:
!          vbas=1up*(2**0)+1do*(2**1)+2up*(2**2)+2do*(2**3)+...
!          vbas->Table(nup, ndo, num), where num is 1,2,...,Nstate(nup,ndo).
!========+=========+=========+=========+=========+=========+=========+=$
     subroutine code_shift(nup, ndo, code, vtot, flagtot)
     use hubbard_param
     implicit none

     integer,intent(IN),dimension(2*Nsite)::code   
	 integer,intent(IN):: nup, ndo 
     integer,intent(OUT), dimension(2*Nsite)::vtot  ! new code
     integer,intent(OUT)::flagtot

	 integer,dimension(Nsite):: vup,vdo, vupnew, vdonew
     integer::flagdo, flagup,i,kk
     integer(kind=4)::count


     if (2*Nsite .ge. 28) then
       print *, 'The number of digits of Table may be insufficient!'
       stop
     endif

!--------code split-------------------------
     do i=1, Nsite
        vup(i)=code(2*i-1)
        vdo(i)=code(2*i)
     enddo
!---------------------------------
  call next1(vdo,vdonew,flagdo)        !a new config for down spin elec.

       if (flagdo .eq.  1) then    !new config obtained successfully.
         do i=1, Nsite
           vtot(2*i-1)=vup(i)
           vtot(2*i)=vdonew(i)
		  enddo
		   flagtot=1
		   return 
       else if (flagdo .eq. 2) then   !error occured.
         write (*,*) 'error in "buildbasis"!, do sector'
	     stop
      endif                          !no more new config for down elec.
            
			! flagdo==0 then shift spin up
!----------------------
      call next1(vup,vupnew,flagup)          !a new config for up spin elec.
      if (flagup .eq. 1) then
         do i=1, Nsite
           vtot(2*i-1)=vupnew(i)
           vtot(2*i)=vdo(i)
         enddo
		   flagtot=1                              !update up elec config, go back to scan down elec.
      else if (flagup .eq. 2) then
        write (*,*) 'error in "buildbasis"!, upsector'
        stop
      endif

!---------------------------
	flagtot=0 

    end subroutine code_shift
!=======================================================================

!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: next  
! NOTICE : 
!
! TYPE   : subroutine 
! PURPOSE: produce a new configuration of (0100110...) with Nsite digit
!          from previous one.(of Nsite size)
! I/O    :
! VERSION: 04-Nov-02
! AUTHOR : N. H. Tong
! COMMENT: This subroutine is used to build basis for Hilbert space
!          where nup and ndo are conserved.
!          vec01(dim):(input)---- integer array, old config (0 1 1 0...). 
!                     (output)--- produced new config (0 1 1 ...).           
!          flag:      (output)--- 0:   no new configuration is ontained. 
!   		                  --- 1:   new config is obtained.
!      			              --- 2:   error occurs. 
!========+=========+=========+=========+=========+=========+=========+=$
	 subroutine next(vec01, flag)
     use hubbard_param
     implicit none
     integer,intent(OUT)::flag
	 integer,dimension(Nsite),intent(INOUT):: vec01
     integer,dimension(Nsite):: v
     integer::dim,i,j,k

     dim=Nsite
	 flag=2
	 do i=1,dim
	   v(i)=vec01(i)
	 enddo

     j=dim
20   if (v(j) .eq. 0) then
       j=j-1
       if (j .eq. 0) then
         flag=0	
! 	(no new config left) 
         goto 50
       endif
       goto 20
     else if(j .ne. dim) then
	   v(j)=0
	   v(j+1)=1
	   flag=1  
!     (new config obtained)
       goto 50
     endif

     k=1
30   j=j-1
     if (j .eq. 0) then
	   flag=0
	   goto 50
     endif
     if (v(j) .ne. 0) then 
       k=k+1
	   goto 30
     endif
!-- now v(j) is 0 on the most right side.-----
40   if (v(j) .eq. 0) then
	   j=j-1
	   if (j .eq. 0) then
          flag=0
          goto 50
	   endif
       goto 40
     endif
     v(j)=0
     do i=1,k+1
	   v(j+i)=1
     enddo
     do i=j+K+2,dim
	   v(i)=0
     enddo
!   (new config obtained)
     flag=1
50   continue
     do i=1, dim
	   vec01(i)=v(i)
     enddo
     end subroutine next
!=====================================================================


!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: next1  
! NOTICE : differ with subroutine 'next' by keep the old vector
!
! TYPE   : subroutine 
! PURPOSE: 
! I/O    :
! VERSION: 28-Dec-07
! AUTHOR : Lei Wang
! COMMENT: This subroutine is used to build basis for Hilbert space
!          where nup and ndo are conserved.
!          vec01(dim):(input)---- integer array, old config (0 1 1 0...). 
!                     (output)--- produced new config (0 1 1 ...).           
!          flag:      (output)--- 0:   no new configuration is ontained. 
!   		                  --- 1:   new config is obtained.
!      			              --- 2:   error occurs. 
!========+=========+=========+=========+=========+=========+=========+=$
	 subroutine next1(vec01, vecout, flag)
     use hubbard_param
     implicit none
     integer,intent(OUT)::flag
	 integer,dimension(Nsite),intent(In):: vec01
     integer,dimension(Nsite),intent(out):: vecout
     integer,dimension(Nsite):: v
     integer::dim,i,j,k

     dim=Nsite
	 flag=2
	 do i=1,dim
	   v(i)=vec01(i)
	 enddo

     j=dim
20   if (v(j) .eq. 0) then
       j=j-1
       if (j .eq. 0) then
         flag=0	
! 	(no new config left) 
         goto 50
       endif
       goto 20
     else if(j .ne. dim) then
	   v(j)=0
	   v(j+1)=1
	   flag=1  
!     (new config obtained)
       goto 50
     endif

     k=1
30   j=j-1
     if (j .eq. 0) then
	   flag=0
	   goto 50
     endif
     if (v(j) .ne. 0) then 
       k=k+1
	   goto 30
     endif
!-- now v(j) is 0 on the most right side.-----
40   if (v(j) .eq. 0) then
	   j=j-1
	   if (j .eq. 0) then
          flag=0
          goto 50
	   endif
       goto 40
     endif
     v(j)=0
     do i=1,k+1
	   v(j+i)=1
     enddo
     do i=j+K+2,dim
	   v(i)=0
     enddo
!   (new config obtained)
     flag=1
50   continue

     do i=1, dim
	   vecout(i)=v(i)
     enddo
     end subroutine next1
!=====================================================================

