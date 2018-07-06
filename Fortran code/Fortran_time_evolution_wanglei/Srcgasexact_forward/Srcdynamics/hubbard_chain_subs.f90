
!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: binary
! NOTICE :
!
! TYPE   : subroutine
! PURPOSE: 
! I/O    : k: integer in [1, 2^{Nsite}], representing basis state |k> (input)
!          bas: 1-dimensional array, containing {0,1} in each element, 
!               the binary expression of the basis state |k>. (output)
! VERSION: 20-Dec-2006
! AUTHOR : Wang Lei
! COMMENT:  copy from the refrence code of Quanising_lanczos.h by N.H.Tong, modified a little.
!          generate the binary expression 
!          |b_1, b_2, ..., b_{Nsite-1}, b_{Nsite} > for a given decimal
!          integer K, where b_n=0 or 1, spin down or up in sigma_x basis.
!          Note: k=1, ==> bas=(0, 0, ..., 0)
!                k=Dim ==> bas=(1, 1, ..., 1)
!========+=========+=========+=========+=========+=========+=========+=$
     subroutine binary(k,bas)
     use hubbard_param
     implicit none
     integer(kind=4), intent(IN):: k
	 integer, intent(OUT), dimension(2*Nsite):: bas
	 integer:: i, mm
	 mm=k
     if ((mm .gt. 4**Nsite) .or. (mm .lt. 1)) then
	    print *, 'Wrong m!, stop!'
		stop
     endif
	 bas=0
	 mm=mm-1                     !for k=1 ==> (0, 0, .., 0)
     i=1
10	 bas(i)=mod(mm,2)
     mm=(mm-bas(i))/2
	 if (mm .eq. 0) then
	    return
     else
	   i=i+1 
	   goto 10
	 endif   
    end subroutine binary


!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: numstate  
! NOTICE :
!
! TYPE   : subroutine 
! PURPOSE: produce the table of number of basises for (nup,ndo) block,
!          which is stored in array nstate(nup,ndo), direct mul 
! I/O    : 
! VERSION: 25-Jan-2007
! AUTHOR : N. H. Tong
! COMMENT: nstate(nup,ndo)=C(nup,Nsite)*C(ndo,Nsite)
!========+=========+=========+=========+=========+=========+=========+=	
     subroutine numstate()
     use hubbard_param
     implicit none
     integer::nup,ndo,i
	 integer(kind=4)::aup,ado,facup,facdo
     do 10 nup=0,Nsite
       do 20 ndo=0,Nsite
         aup=1         
         ado=1
         facup=1
         facdo=1

         do i=Nsite,Nsite-nup+1,-1
           aup=aup*i
         enddo
         do i=1,nup          
           facup=facup*i
         enddo

         do i=Nsite,Nsite-ndo+1,-1
           ado=ado*i
         enddo
         do i=1,ndo          
           facdo=facdo*i
         enddo
         if ((aup/facup*facup .ne. aup) .or. (ado/facdo*facdo .ne. ado))then
           print *, 'Problem in numstate().', nup, ndo
           stop
         endif     
         Nstate(nup,ndo)=(aup/facup)*(ado/facdo)
20    continue
10   continue
    end subroutine numstate
!=======================================================================



!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: numstate1  
! NOTICE :
!
! TYPE   : subroutine 
! PURPOSE: produce the table of number of basises for (nup,ndo) block,
!          which is stored in array nstate(nup,ndo)           
! I/O    : 
! VERSION: 29-Dec-2007
! AUTHOR : Lei Wang, IOP, CAS
! COMMENT: nstate(nup,ndo)=C(nup,Nsite)*C(ndo,Nsite),NO direct mul here
!========+=========+=========+=========+=========+=========+=========+=	
     subroutine numstate1()          ! 
     use hubbard_param
     implicit none
     integer::nup,ndo,i,j
	 integer(kind=4)::aup(Nsite),ado(Nsite),facup(Nsite),facdo(Nsite)
     do 10 nup=0,Nsite
       do 20 ndo=0,Nsite
         
            aup=1
			ado=1
			facup=1
			facdo=1

         do i=1,nup
           aup(i)=Nsite-i+1
		   facup(i)=i
         enddo

         do i=1,ndo  
		   ado(i)=Nsite-i+1
		   facdo(i)=i
         enddo

         do i=1, nup
           do j=1, nup
             if(aup(i)==facup(j)) then 
                 aup(i)=1
				 facup(j)=1        ! reduction common factors 
              endif
            enddo
		  enddo
 
         do i=1, ndo
           do j=1, ndo
             if(ado(i)==facdo(j)) then 
                 ado(i)=1
				 facdo(j)=1
              endif
            enddo
		  enddo
    
  Nstate(nup,ndo)=(PRODUCT(aup(:))/product(facup(:)))*(product(ado(:))/product(facdo(:)))
20    continue
10   continue
    end subroutine numstate1
!=======================================================================

!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: adagger
! TYPE   : subroutine
! PURPOSE: operation of adagger on the basis state.
! I/O    : 
! VERSION: 04-Nov-02
! COMMENT: operation of adag(i,sigma) on the state (nupold,ndoold,jold)
!          produces the state (nupnew,ndonew,jnew) with sign isign.
!          isign=1,and -1. If isign=0, the new state do not exist,nupnew,
!          ndonew and jnew will be 0,0, and 1.
!          iadag: i=1,2, .., Nsite. (input)
!          sigma: +1 (up), -1 (down). (input)
!          nupold, ndoold, jold: the basis state (nupold, ndoold, jold). (input)
!          nupnew, ndonew, jnew: the produced state (also a basis). (output)
!          isign: the sign of the new stata. (output) 
!========+=========+=========+=========+=========+=========+=========+=$
    subroutine adag(iadag,sigma,nupold,ndoold,jold,nupnew,ndonew,jnew,isign)
    use hubbard_param
    implicit none
    integer(kind=4),intent(IN):: jold
    integer(kind=4),intent(OUT):: jnew

	integer,intent(IN)::iadag,nupold,ndoold,sigma
    integer,intent(OUT)::nupnew,ndonew,isign
    integer,dimension(2*Nsite)::basis    
    integer::posit,i
    integer(kind=4),external::lookup
    if ((iadag .lt. 1) .or. (iadag .gt. Nsite)) then
      print *, 'adag i: i exceeds the bound in adag()!'
      stop
    endif
    if ((nupold .lt. 0) .or. (nupold .gt. Nsite)) then
      print *, 'nupold exceeds the bound in adag()!'
      stop
    endif
    if ((ndoold .lt. 0) .or. (ndoold .gt. Nsite)) then
      print *, 'ndoold exceeds the bound in adag()!'
      stop
    endif
    if ((jold .lt. 0) .or. (jold .gt. Nstate(nupold,ndoold))) then
      print *, 'jold exceeds the bound in adag()!'
      stop
    endif
!----assign the nupnew and ndonew. ----
    if (sigma .eq. 1) then
      if (nupold+1 .le. Nsite) then
        nupnew=nupold+1
        ndonew=ndoold
      else
        nupnew=0
        ndonew=0
        jnew=1
        isign=0
        return
      endif
    elseif (sigma .eq. -1) then
      if (ndoold+1 .le. Nsite) then
        ndonew=ndoold+1
        nupnew=nupold
      else
        nupnew=0
        ndonew=0
        jnew=1
        isign=0
        return
      endif 
    else
      print *, 'Problem in value of sigma in adag().'
      stop
    endif
!----find out the jnew----
    if (sigma .eq. 1) then
      posit=2*iadag-1
    else
      posit=2*iadag
    endif
	allocate(Table(Nstate(nupold, ndoold)))
    call buildbasis(nupold, ndoold)
	call findbas(basis,nupold,ndoold,jold) 
	deallocate(Table)   
    if (basis(posit) .eq. 1) then
       nupnew=0
       ndonew=0
       jnew=1
       isign=0
       return
    endif  
     basis(posit)=1
	 allocate(Table(Nstate(nupnew, ndonew)))
      call buildbasis(nupnew, ndonew)
	 jnew=lookup(basis,nupnew,ndonew)
	 deallocate(Table)
!----find out the isign----
    isign=1
    do i=1,posit-1
     if (basis(i) .eq. 1) then
       isign=isign*(-1)
     endif
    enddo
    end subroutine adag  
!=======================================================================


!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: a
! TYPE   : subroutine
! PURPOSE: operation of a
! I/O    : 
! VERSION: 04-Nov-02
! COMMENT: operation of a(i,sigma) on the state (nupold,ndoold,jold)
!          produces the state (nupnew,ndonew,jnew) with sign isign.
!          
!          ia: i=1,2, .., Nsite. (input)
!          sigma: +1 (up), -1 (down). (input)
!          nupold, ndoold, jold: the basis state (nupold, ndoold, jold). (input)
!          nupnew, ndonew, jnew: the produced state (also a basis). (output)
!          isign: the sign of the new stata. (output) 
!                 isign=1 or -1. isign=0 means the new state do not exist.
!                 Then nupnew=0, ndonew=0 and jnew=1.
!========+=========+=========+=========+=========+=========+=========+=$
    subroutine a(ia,sigma,nupold,ndoold,jold,nupnew,ndonew,jnew,isign)
    use hubbard_param
    implicit none

    integer(kind=4),intent(IN):: jold
    integer(kind=4),intent(OUT):: jnew

    integer,intent(IN)::ia,nupold,ndoold,sigma
    integer,intent(OUT)::nupnew,ndonew,isign
    integer,dimension(2*Nsite)::basis    
    integer::posit,i
    integer(kind=4),external::lookup
    if ((ia .lt. 1) .or. (ia .gt. Nsite)) then
      print *, 'adag i: i exceeds the bound in adag()!'
      stop
    endif
    if ((nupold .lt. 0) .or. (nupold .gt. Nsite)) then
      print *, 'nupold exceeds the bound in adag()!'
      stop
    endif
    if ((ndoold .lt. 0) .or. (ndoold .gt. Nsite)) then
      print *, 'ndoold exceeds the bound in adag()!'
      stop
    endif
    if ((jold .lt. 0) .or. (jold .gt. Nstate(nupold,ndoold))) then
      print *, 'jold exceeds the bound in adag()!'
      stop
    endif
    
!--- assign the nupnew and ndonew.---
    if (sigma .eq. 1) then
      if (nupold-1 .ge. 0) then
        nupnew=nupold-1
        ndonew=ndoold
      else
        nupnew=0
        ndonew=0
        jnew=1
        isign=0
        return
      endif
    elseif (sigma .eq. -1) then
      if (ndoold-1 .ge. 0) then
        ndonew=ndoold-1
        nupnew=nupold
      else
        nupnew=0
        ndonew=0
        jnew=1
        isign=0
        return
      endif 
    else
      print *, 'Problem in value of sigma in a().'
      stop
    endif
!--- find out the jnew ----
    if (sigma .eq. 1) then
      posit=2*ia-1
    else
      posit=2*ia
    endif
	allocate(Table(Nstate(nupold, ndoold)))
    call buildbasis(nupold, ndoold)
	call findbas(basis,nupold,ndoold,jold) 
	deallocate(Table)   
    if (basis(posit) .eq. 0) then
       nupnew=0
       ndonew=0
       jnew=1
       isign=0
       return
    endif  
     basis(posit)=0
	 allocate(Table(Nstate(nupnew, ndonew)))
     call buildbasis(nupnew, ndonew)
	 jnew=lookup(basis,nupnew,ndonew)
	 deallocate(Table)
!---- find out the isign ----
    isign=1
    do i=1,posit-1
     if (basis(i) .eq. 1) then
       isign=isign*(-1)
     endif
    enddo
    end subroutine a
!=====================================================


!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: adaggera
! TYPE   : subroutine
! PURPOSE: operation of (a^dagger_i_spin*a_j_spin )
! I/O    : 
! VERSION: 28-Dec-07
! COMMENT: isign is depended on the number of 1s between i and j-site
!          
!          ia: i=1,2, .., Nsite. (input)
!          sigma: +1 (up), -1 (down). (input)
!          nupold, ndoold, jold: the basis state (nupold, ndoold, jold). (input)
!          nupnew, ndonew, jnew: the produced state (also a basis). (output)
!          isign: the sign of the new stata. (output) 
!                 isign=1 or -1. isign=0 means the new state do not exist.
!========+=========+=========+=========+=========+=========+=========+=$
    subroutine adaggera(ia,jadagger,sigma,nupold,ndoold,vecold,vecnew,isign)
    use hubbard_param
    implicit none

     integer,dimension(2*Nsite),intent(IN):: vecold
     integer,dimension(2*Nsite),intent(OUT):: vecnew

    integer,intent(IN)::ia,jadagger,nupold,ndoold,sigma
    integer,intent(OUT)::isign
   
    integer::posit,posit_j, i, shift,fir, sec
    
    if ((ia .lt. 1) .or. (ia .gt. Nsite)) then
      print *, 'adag i: i exceeds the bound in adag()!'
      stop
    endif
    if ((nupold .lt. 0) .or. (nupold .gt. Nsite)) then
      print *, 'nupold exceeds the bound in adaggera()!', nupold
      stop
    endif
    if ((ndoold .lt. 0) .or. (ndoold .gt. Nsite)) then
      print *, 'ndoold exceeds the bound in adagggera()!', ndoold, Nsite
      stop
    endif
 
  vecnew=vecold

!--- find out the jnew ----
    if (sigma .eq. 1) then
      posit=2*ia-1
	  posit_j=2*jadagger-1
    else
      posit=2*ia
	  posit_j=2*jadagger
    endif
   
    if(posit==posit_j) then
	 isign=vecold(posit)
	 return
	endif

    if (vecold(posit) .eq. 0  .or. vecold(posit_j)==1 ) then 
       isign=0
       return
    endif  

     vecnew(posit)=0
	 vecnew(posit_j)=1

!---- find out the isign ----
 
    isign=1
    
	fir=min(posit, posit_j)
	sec=max(posit, posit_j) 
 
    do i=fir+1, sec-1
	 if(vecold(i)==1) then
	 isign=isign*(-1)
	 endif
	enddo     


    end subroutine adaggera
!=====================================================



!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: lookup  
! NOTICE :
!
! TYPE   : function 
! PURPOSE: for a given basis state represented by binary array 'basis',
!          lookup in the Table the corresponding numbering of this 
!          state in certain (nup,ndo) block of a given basis.
! I/O    : 
! VERSION: 04-Nov-02
! AUTHOR : N. H. Tong
! COMMENT:  basis(2*Nsite): (input)--- the given basis state in the (nup,ndo) block.
!                           It is described by 2*Nsite digits of 0 and 1.  
!                        ordered as: 1up,1do,2up,2do,...,Nsup,Nsdo.  
!           lookup: (output)---natural numbering of the basis in (nup,ndo) block.
!========+=========+=========+=========+=========+=========+=========+=$
    function lookup(basis,nup,ndo)
    use hubbard_param
    implicit none
    integer,dimension(2*Nsite),intent(IN):: basis
    integer,intent(IN)::nup,ndo
    integer(kind=4)::lookup
    integer::j,vbas
	integer(kind=4)::i
    vbas=0
    do i=1,2*Nsite
      vbas=vbas+basis(i)*(2**(i-1))  !  trans basis to demical number
    enddo

    do i=1,Nstate(nup,ndo)                                       ! this step is time consuming, !!!!!!!!!!
     
	  if (vbas .eq. Table(i)) then                      
	    lookup=i
        return
      endif
    enddo
    print *, 'Lookup() does not find the given basis in (nup,ndo) block.'
    end function lookup
!=======================================================================
     
!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM:array2integer
! NOTICE :
!
! TYPE   : function 
! PURPOSE: 
! I/O    : 
! VERSION: 07-Dec-28
! AUTHOR : Lei Wang
! COMMENT:  array(2*Nsite): (input)--- the given basis state in the (nup,ndo) block.
!                           It is described by 2*Nsite digits of 0 and 1.  
!                        ordered as: 1up,1do,2up,2do,...,Nsup,Nsdo.  
!           integer: (output)---
!========+=========+=========+=========+=========+=========+=========+=$
    function array2integer(basis)
    use hubbard_param
    implicit none
    integer,dimension(2*Nsite),intent(IN):: basis
    integer(kind=4)::array2integer
	integer::i

    array2integer=0
    do i=1,2*Nsite
      array2integer=array2integer+basis(i)*(2**(i-1))  !  trans basis to demical number
    enddo

 
    end function array2integer
!=======================================================================
     

!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: findbas
! TYPE   : subroutine
! PURPOSE: return the basis array basis(2*Nsite) for a given 
!          quantum number:  (nup,ndo,num).
! I/O    : 
! VERSION: 04-Nov-02
! COMMENT: the invers process of lookup().
!          basis(2*Nsite): the basis state found. (output)
!          nup, ndo, num: the num_th basis in the (nup, ndo) block.
!========+=========+=========+=========+=========+=========+=========+=$    
    subroutine findbas(basis,nup,ndo,num)
    use hubbard_param
    implicit none
    integer,dimension(2*Nsite),intent(OUT)::basis
    integer,intent(IN)::nup,ndo
	integer(kind=4),intent(IN)::num
    integer::value,i
    logical,intrinsic::btest
    if ((nup .lt. 0) .or. (nup .gt. Nsite)) then
      print *, 'nup exceeds the bound in findbas()!'
      stop
    endif
    if ((ndo .lt. 0) .or. (ndo .gt. Nsite)) then
      print *, 'ndo exceeds the bound in findbas()!'
      stop
    endif
    if ((num .lt. 0) .or. (num .gt. Nstate(nup,ndo))) then
      print *, 'num exceeds the bound in findbas()!'
      stop
    endif
    
    value=Table(num)        ! demical rep of a configuration 
    do i=0,2*Nsite-1
      if (btest(value,i)) then
         basis(i+1)=1                ! from right to left , from small to large 
      else
         basis(i+1)=0
      endif
    enddo
    end subroutine findbas
!======================================================================

!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: random
! NOTICE :
!
! TYPE   : subroutine
! PURPOSE: 
! I/O    :
! VERSION: 08-Oct-2006
! AUTHOR : N. H. Tong
! COMMENT: To check the random generator *16807. This function produces
!          random number in [0,1]. This algorithm is written in the book 
!          Computational Physics, by K. H. Hoffmann and M.Schreiber. 
!          It was compared with ranw.f, which is taken from the code
!          of A.Georges etal with the RMP DMFT paper. The effect is not bad.
!          ** Note: to call r16807, first set iseed=odd number. After it is run, 
!          iseed is changed to 2. 
!========+=========+=========+=========+=========+=========+=========+=$
   function random(iseed)
   integer, intent(inout):: iseed
   real(kind=8):: random
   integer:: ran
   save
   if ((iseed+1)/2*2 .eq. iseed+1) then
     iseed=2
     ran=2*iseed+1
   endif
   ran=ran*16807                       !random number in [-2^31, 2^31]
   random=(1.0+4.6566128731e-10*real(ran))/2.0   !random number in [0,1]
   end function random

     function large_pickup(m,n)          !pick up M from N 
     implicit none
	 integer, intent(IN):: m,n 
     integer::i,j
	 integer(kind=4)::a(N), fac(N), large_pickup
       
	   if(m.lt.0)then
       print *,"m less than zero in pick()!"
       stop
       end if
         
       if(m.gt.n)then
         print *,"m larger than n in pick()!"
         stop
       end if

            a=1
			fac=1

         do i=1,m
           a(i)=N-i+1
		   fac(i)=i
         enddo

         do i=1, N
           do j=1, N
             if(a(i)==fac(j)) then 
                 a(i)=1
				 fac(j)=1        ! reduction common factors 
              endif
            enddo
		  enddo
 
  large_pickup=(PRODUCT(a(:))/product(fac(:)))

  if(large_pickup<=0) then 
   print *,  'wrong happens in large_pickup'
   stop
 endif

    end function large_pickup




