!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: Hubbard_ed.f90 
! NOTICE :
! TYPE   : 
! PURPOSE:  Exact diagonalization of the Hubbard Hamiltonian on a 2*2 
!           lattice, with full Hilbert space divided according to 
!           (N_up, N_down) quantum numbers.
!           H=-sum_{i,j, sigma} t_{ij} c_{i \sigma}^{\dagger} c_{i \sigma}
!             +U\sum_{i} n_{\up} n_{\down} - \mu \sum_{i \sigma} n_{i \sigma}
!          
!          output the following depending on the index in the param part:
!             
!        index=2:   local Energy minimum,  Eg and the dgree of degeneracy     
!              3:      <n> vs. Mu 
!              4:     1/chi vs Temp                            
!              5:   Matsubara GF in atomic limit  Lei WangO   
!COMMENT:      (1) generate Table by distribution method 
!              (2) caculate thermodynamic quantities by various subroutines
!              (3) Based on N. H. Tong' codes.
!VERSION: 25-4-2007
! AUTHOR : Lei Wang
!========+=========+=========+=========+=========+=========+=========+=$

!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: prec
! NOTICE :
! TYPE   : module
! PURPOSE: For defining the real kind value Rkind. 
! I/O    :
! VERSION: 25-Dec-2006
! AUTHOR : N. H. Tong
! COMMENT: This module is for RS and Param. Here since main program did not
!          use Rkind, so Rkind=8 should always be used, can not use 16.
!          otherwise there will be problem in transfering data to RS().
!========+=========+=========+=========+=========+=========+=========+=$
   module prec
   implicit none
   integer,parameter::Rkind=8
   end module prec
!========================================================================


!========+=========+=========+=========+=========+=========+=========+=$
 module param
   use prec
   implicit none
!--- lattice information: only one section should be used.----
!   integer, parameter:: Nrank=2, Ncol=4    !geometry of the lattice site.
!   integer, parameter:: Nsite=2*4          !total site number =Nrank *Ncol.
!   integer, parameter:: Nsmax=4900         !(C_8^{4})**2,dimension of subspace (4,4).
!----
!   integer, parameter:: Nrank=2, Ncol=3    !geometry of the lattice site.
!   integer, parameter:: Nsite=2*3          !total site number =Nrank *Ncol.
!   integer, parameter:: Nsmax=400          !(C_6^{3})**2, dimension of (3,3).
!----
   integer, parameter:: Nrank=1, Ncol=6    !geometry of the lattice site.
   integer, parameter:: Nsite=Nrank*Ncol         !total site number =Nrank *Ncol.
   integer, parameter:: Nsmax=400           !(C_4^{2})**2, the dimension of the 
                                           !largest subspace (2,2).
!--- model parameters -------------
   real(kind=Rkind):: Thop !hopping amplitude
   real(kind=Rkind):: U    !on-site U 
   real(kind=Rkind):: Mu     !chemical potential 
   real(kind=Rkind):: Temp  !temperatuer Temp>0 only.
   real(kind=Rkind):: Beta                 !inverse temperature.
!--- calculation parameters -------
   integer, parameter:: Iwmax=600           !Matsubara frequencies number.
   integer, parameter:: Ntau=100            !number of points of imaginary time tau for G(tau).
   integer, parameter:: Nomega=301           !number of omega points in DOS, using odd number.
   real(kind=Rkind), parameter:: Omegac=3.0_Rkind   !largest omega of DOS.
   real(kind=Rkind), parameter:: Broaden=0.05_Rkind   !eta, for broadening of DOS peaks.
   integer, parameter:: Site=1              !the site on which Gfs are calculated 1<=Site<=Nsite.
!--- variables ---------------------
   integer(kind=4), dimension(0:Nsite, 0:Nsite, Nsmax):: Table  !The table of numbering of states.
   integer(kind=4), dimension(0:Nsite, 0:Nsite):: Nstate        !containing the dimension D(nup, ndo)
   real(kind=Rkind), dimension(Nsite, Nsite):: Tmn              !hopping matrix.
   real(kind=Rkind), dimension(0:Nsite, 0:Nsite, Nsmax):: Eval  !eigen values from RS().
   real(kind=Rkind), dimension(Nsite):: Avnup, Avndo    !electron number averages.
   real(kind=Rkind):: Chiuni                  !uniform susceptibility.
   real(kind=Rkind), dimension(Nomega):: Dosup, Dosdo      !local DOS on the site Msite.
   complex(kind=Rkind), dimension(0: Iwmax):: Gwup, Gwdo   !local Matsubara GF on the site Msite.         
   real(kind=Rkind), parameter:: Pi=3.1415926535_Rkind     !Pi
   complex(kind=Rkind), parameter:: Xi=(0.0_Rkind, 1.0_Rkind)    !imaginary unit!    
   real(kind=Rkind):: Eg                                     !globle ground state energy 
   real(kind=Rkind),dimension((Nsite+1)**2) :: Egs        !local ground state energys
   integer::Ndeg  
!------------index for assignments------------------------------------    
   integer, parameter:: index=2                                 ! index = 2;  do assignment (2):   find Eg and the dgree of degeneracy          
                                                                  !       3      assignment (3): <n> vs. Mu
                                                                  !       4       assignment (4): 1/chi vs Temp                                                 
  end module param                                                !       5        assignment (5)  Matsubara GF in atomic limit 
!==============================================================


!========+=========+=========+=========+=========+=========+=== ======+=$
! PROGRAM: Hubbard_ed
! NOTICE :
! TYPE   : main
! PURPOSE: For exact diagonalization of the Hubbard model on a 2*4 lattice.
! I/O    :
! VERSION: 16-Apl-2007
! AUTHOR : Wang Lei
! COMMENT: 
!========+=========+=========+=========+=========+=========+=========+=$
     program Hubbard_ed
     use param
     implicit none
     integer::nup,ndo,n,m,i
	 real(kind=Rkind):: avntot, omegan, avnat, taun, sum

	
!--- check parameters -------------
	 if ((Site .lt. 1) .or. (Site .gt. Nsite)) then
	   print *, 'Site is wrong!'
	   stop
     endif


!--------(1) generate Table()--------------------------------------
      call numstate()     !produce Nstate(nup, ndo).
      call distribute()   !produce Table(nup, ndo, m)
!--- check Nstate(nup, ndo) -----------------------------
      do nup=0, Nsite
        do ndo=0, Nsite
         if (Nstate(nup,ndo) .gt. Nsmax) then
           print *,'Nsmax is smaller than Nstate(nup,ndo).'
           print *,'There is problem!'
           stop
         endif
        enddo
      enddo
!------------------------------------------------------
!---------goto assignment dependes on the index-----------
 if (index.eq.2) then
     goto 20                                                      ! assigment (2): produce H and diagonlize it, find Eg and the dgree of degeneracy 
     elseif (index.eq.3)then
	 goto 30                                                      ! assigment (3): <n> vs. Mu
	 elseif (index.eq.4)then
	 goto 40                                                      ! assigment (4): 1/chi vs Temp 
	 elseif (index.eq.5)then
	 goto 50                                                      ! (5): G_{up} ant atom limit
	 
endif	                                                          
!--------------------------------------------------------------------------------

20  continue 
!--- (2) Local ground state energy and the dgree of degenerary -----
!-----parameter-----------------------------------
     Thop=1.0_Rkind
	 Mu=0.0_Rkind
	 U=0.0_Rkind
!----------------------------------------------------
	 call hamdiag()                                                  ! Eg is the global ground state energy

!----------caculate the # of degeneracy---------------
	 Ndeg=0  
     do nup=0,Nsite
       do ndo=0,Nsite
         do i=1,Nstate(nup,ndo)
             if (abs(Eval(nup,ndo,i)) .lt. 1e-4) then              ! note Ndeg at least increse 1 time and all the Evals are substract by Eg already      
             Ndeg=Ndeg+1                                           ! if there is another energy level close enough to Eg, degeneracy +1, 
            endif
           enddo   
         enddo
       enddo 

!----Output local ground states---------
     open (3, file='result/locen.dat', form='formatted')
	     write (3,"('t='(f15.6))")  Thop  
		 write (3,"('Mu='(f15.6))")  Mu  
	     write (3,"('U='(f15.6))")  U
		 write (3,"('-----------------------------------------------------')") 
	  
	     write (3,"('Nup|Ndo|Local Energy minimum')") 
	 do nup=0, Nsite
	   do ndo=0, Nsite
	     write (3,"(i3,i3,f15.6)") nup, ndo ,   Egs((nup+1)*(ndo+1)) 
	   enddo
	 enddo
      write (3,"('-----------------------------------------------------')") 
	  
	  write (3,"('Global Ground State Energy:'(f15.6))") Eg
	  write (3,"('Ndeg='(i3))")    Ndeg
	 
	 endfile 3
	 close(3)
         

goto 300                                                                    ! finish this assignment
  

30 continue
!-------(3) <n>~Mu-------------------------------- 
       open (4, file='result/nave.dat', form='formatted')

	 Temp=0.01_Rkind
	 Thop=1.0_Rkind

 do 90 U=0.0_Rkind, 4.0_Rkind, 2.0_Rkind                       !  scan the parameter space of U  
	  	 write (4,"('-----------------------------------------------------')") 
	       write (4,"('Temp='(f15.6))")  Temp
		   write (4,"('U='(f15.6))")    U
           write (4,"('t='(f15.6))")  Thop     
	  write (4,"('-----------------------------------------------------')") 
      write (4, "('       Mu            <n>')")        

 do 100 Mu=0.15_Rkind, 2.5_Rkind, 0.1_Rkind                   ! scan the parameter space of Mu 
     
	 call hamdiag() 

 
!---- calculate Green's functions and averages(does not depend on dagmat() ).----------- 	
     Beta=1.0_Rkind/Temp    !for finite temperature only!   
	
     m=Site
	 call avn(m)     

       avntot=0.0_Rkind
     do m=1, Nsite
	   avntot=avntot+Avnup(m)+Avndo(m)
	 enddo
      avntot=avntot/real(Nsite, Rkind)

         write (4, "(2(f15.6))") Mu, avntot

100 continue
90 continue

   endfile 4
	close(4)
 
goto 300                                             ! finish this assignment


!---------------------------------------------------------
40 continue 
!------------(4) chi_u vs Temp------------------
!-------parameter for this assingment-------------------------------    
    m=Site
	U=3.0_Rkind
    Mu=1.5_Rkind
!--------------------------------------------------------------------   	 
  open (5, file='result/chi.dat', form='formatted') 
  
  do Thop=0.5_Rkind, 1.0_Rkind, 0.1_Rkind                              ! Thop large than 0.8;  AF 
     call hamdiag()  
     
     call dagmat(m)
 
   write (5,"('-----------------------------------------------------')")  
   write (5,"('Thop='(f15.6))") Thop
 
   do 120 Temp=0.01_Rkind, 1.0_Rkind, 0.05_Rkind                        ! temprature parameter space 
    
     Beta=1.0_Rkind/Temp                                      !for finite temperature only!   
	  call chiu(m)  

   write (5, "(2(f15.6))") Temp, 1.0_Rkind/Chiuni

  120 continue
 
  enddo 
   endfile 5
    close(5)
 

goto 300                                              ! finish this assignment


50 continue

!------------------------------------------------------
!--- (5) Output  atom limit Matsubara GF into file 'gwn.dat'-----------

!-----parameter for this assigment-----------------------------------    
	 Thop=0.0_Rkind
	  U=1.0_Rkind
	  Mu=0.5_Rkind
	 Temp=1.0_Rkind
!--------------------------------------------------------------------
	 Beta=1.0_Rkind/Temp  

     call hamdiag()  
     m=Site  
	 call dagmat(m)
     call matgf(m)                              ! caculate the matsubara GF   
	 
	  open (6, file='result/matgf.dat', form='formatted')
	 do n=0, Iwmax
	   omegan=real(2*n+1,Rkind)*Pi/Beta
	  write (6, 200) omegan, real(gwup(n)), aimag(gwup(n))         ! the real and imagine part of the Matsubara GF
	 enddo
      print *, 'the matusubara GF has been caculate, output to matgf.data'	
	 endfile 6
	 close(6)

!--- Output DOS into file 'dos.dat' -----------
     sum=0.0_Rkind    !check sum rule of DOS.
     open (7, file='result/dos.dat', form='formatted')
	 do n=1, Nomega
	   omegan=-Omegac+real(n-1, Rkind)/real(Nomega-1, Rkind)*2.0*Omegac
	   sum=sum+Dosup(n)/real(Nomega-1, Rkind)*2.0*Omegac
	   write (7, 200) omegan, Dosup(n), Dosdo(n)
     enddo
	 
	 endfile 7
	 close(7)
	 print *, 'DOS sum rule=',sum    !should be 1.0 if the integral is exact.
	 print *, 'DoS output to dos.dat'

200   format (1x, 3(e18.8,2x))
300   continue



!------------------------------------------------------
end program Hubbard_ed









!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: hamdiag
! TYPE   : subroutine
! PURPOSE: construc the hamiltonian matrix for each (nup,ndo) block and
!          diagonalize them. Output eigenvalues and eigenvectors into
!          external files.
! I/O    : 
! VERSION: 26-Jan-2007
! COMMENT:
! AUTHOR:  N.H.Tong
!========+=========+=========+=========+=========+=========+=========+=$ 
     subroutine hamdiag()
     use param
     implicit none
     integer:: nup,ndo,i,j, ms, ns, dim, spin, nupnew, ndonew, jnew, is
     integer:: nupmid, ndomid, jmid, ismid, sum,m,n
	 integer,dimension(2*Nsite):: basis
   
	 

	 real(kind=Rkind):: sym  !for testing
     character(len=10):: filename
!--- for Rs() subroutine. ------
     integer:: cevec=1  !always 1: RS() produces eigen-vector.
	 integer:: ierr
     real(kind=Rkind),dimension(Nsmax):: eigvalue, fv1, fv2
     real(kind=Rkind),dimension(Nsmax,Nsmax):: eigvector, hamilt
!------------------------------------------------------------
!--- assign the hopping matrix ----
     Tmn=0.0_Rkind
     if ((Nrank .eq. 2) .and. (Ncol .eq. 2)) then
       Tmn(1,2)=1.0_Rkind;  Tmn(3,4)=1.0_Rkind
	   Tmn(1,3)=1.0_Rkind;  Tmn(2,4)=1.0_Rkind
     elseif ((Nrank .eq. 2) .and. (Ncol .eq. 3)) then
       Tmn(1,2)=1.0_Rkind;  Tmn(2,3)=1.0_Rkind
	   Tmn(4,5)=1.0_Rkind;  Tmn(5,6)=1.0_Rkind
       Tmn(1,4)=1.0_Rkind;  Tmn(2,5)=1.0_Rkind;  Tmn(3,6)=1.0_Rkind
     elseif ((Nrank .eq. 2) .and. (Ncol .eq. 4)) then
       Tmn(1,2)=1.0_Rkind;  Tmn(2,3)=1.0_Rkind;  Tmn(3,4)=1.0_Rkind
	   Tmn(5,6)=1.0_Rkind;  Tmn(6,7)=1.0_Rkind;  Tmn(7,8)=1.0_Rkind
       Tmn(1,5)=1.0_Rkind;  Tmn(2,6)=1.0_Rkind;  Tmn(3,7)=1.0_Rkind
	   Tmn(4,8)=1.0_Rkind
     elseif ((Nrank .eq. 1) .and. (Ncol .eq. 4)) then
       Tmn(1,2)=1.0_Rkind;  Tmn(2,3)=1.0_Rkind;  Tmn(3,4)=1.0_Rkind
	   Tmn(1,4)=1.0_Rkind
	  elseif ((Nrank .eq. 1) .and. (Ncol .eq. 6)) then
       Tmn(1,2)=1.0_Rkind;  Tmn(2,3)=1.0_Rkind;  Tmn(3,4)=1.0_Rkind
	   Tmn(4,5)=1.0_Rkind; Tmn(5,6)=1.0_Rkind; Tmn(1,6)=1.0_Rkind
	  else
	    print *, 'Nrank and Ncol have problem.'
		stop
      endif
      Tmn=Tmn*Thop
	 
	 do m=1, Nsite
	  do n=1, m-1
	     Tmn(m,n)=Tmn(n,m)
      enddo
	  enddo

     do 10 nup=0,Nsite
     do 20 ndo=0,Nsite
       dim=Nstate(nup,ndo)
	   print *, 'U=', U
	   print *, 'Mu=', Mu
       print *, 'subspace (Nup, Ndo)=', Nup, Ndo
	   print *, 'Dimension=', dim
!--- clear ----
       hamilt=0.0_Rkind
       eigvalue=0.0_Rkind
       eigvector=0.0_Rkind
       fv1=0.0_Rkind
       fv2=0.0_Rkind
!-------- construct hamiltonian matrix ---------
        do 25 j=1,dim
!---- hopping term ----
          do 30 spin=-1,1,2
           do 40 ms=1, Nsite
		   do 50 ns=1, Nsite
			 if (abs(Tmn(ms,ns)) .gt. 1.0e-14) then
               call a(ns,spin,nup,ndo,j,nupmid,ndomid,jmid,ismid)
               if (ismid .eq. 0) then
                   goto 50
               endif
               call adag(ms,spin,nupmid,ndomid,jmid,nupnew,ndonew,jnew,is)
               hamilt(jnew,j)=hamilt(jnew,j)+real(ismid*is)*(-Tmn(ms,ns))
             endif
50         continue
40         continue
30       continue
!---- U term ----
	   call findbas(basis, nup, ndo, j)
	   sum=0
       do ms=1, Nsite
	      sum=sum+basis(2*ms-1)*basis(2*ms)
       enddo
	   hamilt(j,j)=hamilt(j,j)+U*real(sum, Rkind)
!---- Mu term -----
       sum=0
	   do ms=1, Nsite
	      sum=sum+basis(2*ms-1)+basis(2*ms)
       enddo
	   hamilt(j,j)=hamilt(j,j)-Mu*real(sum, Rkind)
25    continue
!----------- H matrix completed! ------------------

!----- check Hermitian of H matrix ---------
     do i=1,dim
      do j=i+1,dim
        sym=sym+dabs(hamilt(i,j)-hamilt(j,i))
      enddo
     enddo
     if (sym .ne. 0.0) then
       print *, 'not symmetric hamilt, --- problem!'
       stop
     endif
!     print *, 'symmetricity=',sym

!------ diagonalize Hamilt -----------------
     call  rs (Nsmax,dim,hamilt,eigvalue, &
             cevec, eigvector, fv1, fv2, ierr)
     if (ierr .ne. 0) then
       print *, 'There is problem in diagonalizing Hamilt &
                & using RS() in subroutine hamdiag().'
       stop
     endif
!--------------------------------------------
print *, 'sub_gsenergy', eigvalue(1)
pause

!--- store eigenvalues and eigenvectors -----
     do i=1,dim
        Eval(nup,ndo,i)=eigvalue(i)
     enddo
	 Egs((nup+1)*(ndo+1))=eigvalue(1)                             ! put the local ground state energys in Egs()
	 
 
     write(filename,'(a2,2I1,a6)') 'zz',nup,ndo,'evfile'
     open(99,FORM='UNFORMATTED',file=filename,status='unknown')
     write(99)((eigvector(i,j),i=1,dim), j=1,dim)
     close(99)
!---------------------------------------------
20    continue
10   continue

	 
!------ find out globle ground state energy Eg ----
!------ and subtract it from En. -----------------
                                                          
	  Eg=minval(Egs)                                                  !ground state energy
 
       do nup=0,Nsite
       do ndo=0,Nsite
         do i=1,Nstate(nup,ndo)
           Eval(nup,ndo,i)=Eval(nup,ndo,i)-Eg                      ! subtract the ground state energy from all the energy levels, crucial for avoiding overflow of exp!
         enddo
       enddo
    enddo	   
end subroutine hamdiag
!========================================================================


!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: distribute 
! NOTICE :
!
! TYPE   : subroutine
! PURPOSE: produce Table(nup, ndo, m) which stores the decimal 
!          representation of the num_th basis in (nup, ndo) block. 
!          for (0, 0, .., 0)(2Nsite bits) ==> k=0 , no electron ; for (1, 1, .., 1) ==> k=4**Nsite-1 , electron up and down on every site 
! I/O    :
! VERSION: 16-Apl-2007
! AUTHOR : Wang Lei
! COMMENT: 1 <= num <=Nstate(nup,ndo): an integer number (input)
!          Table : the decimal integer representing the num_th state.
!          A  basis state is represented by 2*Nsite digits in this way:
!          | 1up,1do,2up,2do,...,Nsup,Nsdo >, where 1up,1do...=0,1.
!          It is uniquely recorded by an integer:
!          vbas=1up*(2**0)+1do*(2**1)+2up*(2**2)+2do*(2**3)+...
!          vbas->Table(nup, ndo, num), where num is 1,2,...,Nstate(nup,ndo).
!========+=========+=========+=========+=========+=========+=========+=$
subroutine distribute()
     use param
     implicit none
     integer,dimension(2*Nsite)::vtot
     integer::nup,ndo,i,count                            !nup,ndo count the number of spin up/down electrons for a state k
	 integer(kind=4)::k
     integer, dimension(0:Nsite, 0:Nsite) :: Ncount      ! the sequence number of the (nup, ndo) subspace
    


     Table=0
	 Ncount=1                                    
	 
  
  do k=1, 4**Nsite                         ! decimal number represent of a state 
     nup=0                                 !!!
	 ndo=0
      call binary(k,vtot)                    ! binary number represent of the state
	 
	  do i=1, Nsite
	    nup=nup+vtot(2*i-1)                  ! count the # of up and down spin of this state
        ndo=ndo+vtot(2*i)
      enddo
  
   Table(nup, ndo, Ncount(nup,ndo))=k-1         !for (0, 0, .., 0) ==> k=0 
   Ncount(nup, ndo)=Ncount(nup, ndo)+1           !sequence number in this space increase 1
   
   !  print * ,'nup', nup
   !  print * ,'ndo', ndo
   ! print * ,'Ncount=', Ncount(nup,ndo)
 enddo

 end subroutine distribute


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
     use param
     implicit none
     integer, intent(IN):: k
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
! PROGRAM: avn
! TYPE   : subroutine
!
! I/O    :  Input: (msite)      output: average particle number on msite.
! VERSION: 23-Apl-2007
!========+=========+=========+=========+=========+=========+=========+=$ 
      subroutine avn(msite)
      use param
      implicit none
	  integer, intent(IN):: msite
      integer:: i, j, n, m, nup, ndo, dim, sr, sc
      real(kind=Rkind),dimension(Nsmax,Nsmax):: ddag, evec
	  integer, dimension(2*Nsite):: basis
	 
	  real(kind=Rkind), dimension(Nsite):: sumup, sumdo
      real(kind=Rkind):: e0, partf, omegan, sumuni, mag, taun, aij
      character(len=10):: filename
!-- partition function at Temp ----------------------
      partf=0.0_Rkind
      do nup=0,Nsite
       do ndo=0,Nsite
        do i=1,Nstate(nup,ndo)
          partf=partf+exp(-Beta*Eval(nup,ndo,i))
        enddo
       enddo
      enddo

      Avnup=0.0_Rkind
	  Avndo=0.0_Rkind
	 

      do 100 nup=0, Nsite
	  do 120 ndo=0, Nsite
	    dim=Nstate(nup, ndo)
!--- read from files ----
        write(filename,'(a2,2I1,a6)') 'zz', nup, ndo, 'evfile'
        open(99,FORM='UNFORMATTED',file=filename,status='old')
        read(99) ((evec(i,j),i=1,dim),j=1,dim)
        close(99)

        do 140 n=1, dim
!--- (1) n_up and n_do for each site: Avnup, Avndo ----------
          sumup=0.0_Rkind
		  sumdo=0.0_Rkind
		  do m=1, dim
		     call findbas(basis, nup, ndo, m)
		     do i=1, Nsite
                sumup(i)=sumup(i)+basis(2*i-1)*evec(m,n)**2
                sumdo(i)=sumdo(i)+basis(2*i)*evec(m,n)**2
             enddo
           enddo

		   do i=1, Nsite
             Avnup(i)=Avnup(i)+exp(-Beta*Eval(nup, ndo, n))*sumup(i)
			 Avndo(i)=Avndo(i)+exp(-Beta*Eval(nup, ndo, n))*sumdo(i)
		   enddo          
!--------------------------------------------------------
140     continue
120   continue
100   continue
!------------------------------------------
      do i=1, Nsite
	    Avnup(i)=Avnup(i)/partf
		Avndo(i)=Avndo(i)/partf
	  enddo
end subroutine avn









!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM:  chiu
! TYPE   : subroutine
!
! I/O    :  Input: (msite)      output: uniform magnetization
! VERSION: 23-Apl-2007
!========+=========+=========+=========+=========+=========+=========+=$ 
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine chiu(msite)
      use param
      implicit none
	  integer, intent(IN):: msite
      integer:: i, j, n, m, nup, ndo, dim
      real(kind=Rkind),dimension(Nsmax,Nsmax):: ddag, evec
	  integer, dimension(2*Nsite):: basis
	 
	  real(kind=Rkind), dimension(Nsite):: sumup, sumdo
      real(kind=Rkind):: e0, partf, omegan, sumuni, mag, taun, aij
      character(len=10):: filename
!-- partition function at Temp ----------------------
      partf=0.0_Rkind
      do nup=0,Nsite
       do ndo=0,Nsite
        do i=1,Nstate(nup,ndo)
          partf=partf+exp(-Beta*Eval(nup,ndo,i))
        enddo
       enddo
      enddo


!--------------------------------------------------------------------     
	  Chiuni=0.0_Rkind
      do 160 nup=0, Nsite
	  do 180 ndo=0, Nsite
	    dim=Nstate(nup, ndo)
!--- read from files ----
        write(filename,'(a2,2I1,a6)') 'zz', nup, ndo, 'evfile'
        open(99,FORM='UNFORMATTED',file=filename,status='old')
        read(99) ((evec(i,j),i=1,dim),j=1,dim)
        close(99)

        do 190 n=1, dim
!-- (2) uniform susceptibility chiuni and staggered susceptibility chistag ------
          sumuni=0.0_Rkind
  		 
          do m=1, dim
		    call findbas(basis, nup, ndo, m)
			mag=0.0_Rkind
		    do i=1, Nsite
			   mag=mag+basis(2*i-1)-basis(2*i)
            enddo
             sumuni=sumuni+(mag*evec(m,n))**2/4.0
			
           enddo
           Chiuni=Chiuni+exp(-Beta*Eval(nup, ndo, n))*sumuni
	     
!--------------------------------------------------------
190    continue
180   continue
160   continue
!------------------------------------------

	  Chiuni=Chiuni*Beta/(partf*real(Nsite))

end subroutine chiu



!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: matgf
! TYPE   : subroutine
!
! I/O    :  Input: (msite)      output: matsubara green's func on msite
! VERSION: 23-Apl-2007
!========+=========+=========+=========+=========+=========+=========+=$ 
subroutine matgf(msite)

      use param
      implicit none
	  integer, intent(IN):: msite
      integer:: i, j, n, m, nup, ndo, dim
      real(kind=Rkind),dimension(Nsmax,Nsmax):: ddag, evec
	  integer, dimension(2*Nsite):: basis
	 
	  real(kind=Rkind), dimension(Nsite):: sumup, sumdo
      real(kind=Rkind):: e0, partf, omegan, sumuni, mag, taun, aij
      character(len=10):: filename
!----------------------------------------------------------------------
!-- partition function at Temp ----------------------
      partf=0.0_Rkind
      do nup=0,Nsite
       do ndo=0,Nsite
        do i=1,Nstate(nup,ndo)
          partf=partf+exp(-Beta*Eval(nup,ndo,i))
        enddo
       enddo
      enddo
!--------- calculating Gwup(i*omega) ----- (0 < omega < omega_c) ;   --------

!----------calculating DOS(omega) ------ (-Omegac0 < omega < Omegac) ; ---------------

      Gwup(:)=(0.0_Rkind,0.0_Rkind)
     
      Dosup(:)=0.0_Rkind
       
      do 10 nup=0,Nsite-1
      do 20 ndo=0,Nsite         
!--- Read from file ----
       write(filename,'(a2,2I1,a6)') 'zz',nup,ndo,'dupfil'
       open(99,FORM='UNFORMATTED',file=filename,status='unknown')
       read(99)((ddag(i,j),i=1,Nstate(nup+1,ndo)), j=1,Nstate(nup,ndo))
       close(99)
!--- produce Gwup and Dosup, together -----
      do j=1,Nstate(nup,ndo)
      do i=1,Nstate(nup+1,ndo)
         aij=ddag(i,j)**2*(exp(-Beta*Eval(nup,ndo,j))+exp(-Beta*Eval(nup+1,ndo,i)))
!--- Gwup ----
	     do n=0, Iwmax	  
		   omegan=real(2*n+1,Rkind)*Pi/Beta 
           Gwup(n)=Gwup(n)+aij/(Xi*omegan+Eval(nup,ndo,j)-Eval(nup+1,ndo,i))
         enddo
!--- Dosup ----         
	     do n=1, Nomega
            omegan=-Omegac+real(n-1, Rkind)/real(Nomega-1)*2.0_Rkind*Omegac
            Dosup(n)=Dosup(n)+aij/((omegan+Eval(nup,ndo,j)-Eval(nup+1,ndo,i))**2+Broaden**2)
         enddo
	   enddo
       enddo

20    enddo
10    enddo

      do n=0,Iwmax
        Gwup(n)=Gwup(n)/partf
      enddo
	 
	  do n=1, Nomega
	    Dosup(n)=Dosup(n)*Broaden/(Pi*partf)
      enddo

 
!---------- calculating Gwdo(i*omega)  -----------
!---------- (0 < omega < omega_c) ----------
      Gwdo(:)=(0.0_Rkind,0.0_Rkind)
    
	  Dosdo(:)=0.0_Rkind

      do 30 nup=0,Nsite
      do 40 ndo=0,Nsite-1         
!--- Read from file ----
       write(filename,'(a2,2I1,a6)') 'zz',nup,ndo,'ddofil'
       open(99,FORM='UNFORMATTED',file=filename,status='unknown')
       read(99)((ddag(i,j),i=1,Nstate(nup,ndo+1)), j=1,Nstate(nup,ndo))
       close(99)

!--- produce Gwdo and Dosdo, together ------
      do j=1,Nstate(nup,ndo)
      do i=1,Nstate(nup,ndo+1)
	      aij=ddag(i,j)**2*(exp(-Beta*Eval(nup,ndo,j))+exp(-Beta*Eval(nup,ndo+1,i))) 
!--- Gwdo -----
         do n=0,Iwmax
           omegan=real(2*n+1, Rkind)*Pi/Beta 
           Gwdo(n)=Gwdo(n)+aij/(xi*omegan+Eval(nup,ndo,j)-Eval(nup,ndo+1,i))
         enddo
!--- Dosdo ----
         do n=1, Nomega
	       omegan=-Omegac+real(n-1, Rkind)/real(Nomega-1)*2.0_Rkind*Omegac
           Dosdo(n)=Dosdo(n)+aij/((omegan+Eval(nup,ndo,j)-Eval(nup,ndo+1,i))**2+Broaden**2)
         enddo
      enddo
      enddo


!--------------
40    enddo
30    enddo
      
	  do n=0,Iwmax
        Gwdo(n)=Gwdo(n)/partf
      enddo
    
	  do n=1, Nomega
	     Dosdo(n)=Dosdo(n)*Broaden/(Pi*partf)
      enddo
!---------------------------------------------------------      
	  
	  
	 
end subroutine matgf
!--------------------------------------------------------------------------



!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: numstate  
! NOTICE :
!
! TYPE   : subroutine 
! PURPOSE: produce the table of number of basises for (nup,ndo) block,
!          which is stored in array nstate(nup,ndo)
! I/O    : 
! VERSION: 25-Jan-2007
! AUTHOR : N. H. Tong
! COMMENT: nstate(nup,ndo)=C(nup,Nsite)*C(ndo,Nsite)
!========+=========+=========+=========+=========+=========+=========+=	
     subroutine numstate()
     use param
     implicit none
     integer::nup,ndo,i,aup,ado,facup,facdo
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
           print *, 'Problem in numstate().'
           stop
         endif     
         Nstate(nup,ndo)=(aup/facup)*(ado/facdo)
20    continue
10   continue
    end subroutine numstate
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
    use param
    implicit none
    integer,intent(IN)::iadag,nupold,ndoold,jold,sigma
    integer,intent(OUT)::nupnew,ndonew,jnew,isign
    integer,dimension(2*Nsite)::basis    
    integer::posit,i
    integer,external::lookup
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
    call findbas(basis,nupold,ndoold,jold)    
    if (basis(posit) .eq. 1) then
       nupnew=0
       ndonew=0
       jnew=1
       isign=0
       return
    endif  
     basis(posit)=1
     jnew=lookup(basis,nupnew,ndonew)
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
    use param
    implicit none
    integer,intent(IN)::ia,nupold,ndoold,jold,sigma
    integer,intent(OUT)::nupnew,ndonew,jnew,isign
    integer,dimension(2*Nsite)::basis    
    integer::posit,i
    integer,external::lookup
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
    call findbas(basis,nupold,ndoold,jold)    
    if (basis(posit) .eq. 0) then
       nupnew=0
       ndonew=0
       jnew=1
       isign=0
       return
    endif  
     basis(posit)=0
     jnew=lookup(basis,nupnew,ndonew)
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
    use param
    implicit none
    integer,dimension(2*Nsite),intent(IN):: basis
    integer,intent(IN)::nup,ndo
    integer::lookup
    integer::j,i,vbas
    vbas=0
    do i=1,2*Nsite
      vbas=vbas+basis(i)*(2**(i-1))
    enddo

    do i=1,Nstate(nup,ndo)
      if (vbas .eq. Table(nup,ndo,i)) then
	    lookup=i
        return
      endif
    enddo
    print *, 'Lookup() does not find the given basis in (nup,ndo) block.'
    end function lookup
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
    use param
    implicit none
    integer,dimension(2*Nsite),intent(OUT)::basis
    integer,intent(IN)::nup,ndo,num
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
    
    value=Table(nup,ndo,num)
    do i=0,2*Nsite-1
      if (btest(value,i)) then
         basis(i+1)=1
      else
         basis(i+1)=0
      endif
    enddo
    end subroutine findbas
!======================================================================


!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: dagmat
! TYPE   : subroutine
! PURPOSE: calculate the matrix of a_{msite, up}^{+} and a_{msite, do}^{+}
!          for the msite_th site in the H-diagonalized space. 
!          Then output matrixes into files.
! I/O    : 
! VERSION: 27-Jan-2007
! COMMENT: msite: the site where the a^{\dagger} is to be calculated.
!          1< msite < Nsite, according to the lattice site numbering fig.  
!========+=========+=========+=========+=========+=========+=========+=$ 
      subroutine dagmat(msite)
      use param
      implicit none
	  integer,intent(IN):: msite   
      integer:: i, j, dim, dimp, nup, ndo, nupnew, ndonew, jnew, isign, jold
      real(kind=8),dimension(Nsmax,Nsmax):: evec, evecp, matrix
      character(len=10):: filename
      
!--- for ddagger(up) -----
      do nup=0, Nsite-1
       do ndo=0, Nsite

         matrix=0.0_Rkind
         dim=Nstate(nup,ndo)
         dimp=Nstate(nup+1,ndo)
!--- read from files ----
         write(filename,'(a2,2I1,a6)') 'zz', nup, ndo, 'evfile'
         open(99,FORM='UNFORMATTED',file=filename,status='old')
         read(99)((evec(i,j),i=1,dim),j=1,dim)
         close(99)
         write(filename,'(a2,2I1,a6)') 'zz', nup+1, ndo, 'evfile'
         open(99,FORM='UNFORMATTED',file=filename,status='old')
         read(99)((evecp(i,j),i=1,dimp),j=1,dimp)
         close(99)

         do j=1, dim
         do i=1, dimp       
          do 10 jold=1, dim
            call adag(msite,1,nup,ndo,jold,nupnew,ndonew,jnew,isign)
            if (isign .eq. 0) then
             goto 10
            endif
            if ((nupnew .ne. nup+1) .or. (ndo .ne. ndonew)) then
             print *, 'There is problem.'
             stop
            endif
            matrix(i,j)=matrix(i,j)+real(isign)*evec(jold,j)*evecp(jnew,i)
10        continue          
         enddo
         enddo
!---- Output into file -----
         write(filename,'(a2,2I1,a6)') 'zz',nup,ndo,'dupfil'
         open(99,FORM='UNFORMATTED',file=filename,status='unknown')
         write(99)((matrix(i,j),i=1, dimp), j=1, dim)
         close(99)
        enddo
       enddo

!---- for ddagger (do) ------
      do nup=0, Nsite
       do ndo=0, Nsite-1
         matrix=0.0_Rkind
         dim=Nstate(nup,ndo)
         dimp=Nstate(nup,ndo+1)
!---- read from files -------
         write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'evfile'
         open(99,FORM='UNFORMATTED',file=filename,status='old')
         read(99)((evec(j,i),j=1,dim),i=1,dim)
         close(99)
         write(filename,'(a2,2I1,a6)')'zz',nup,ndo+1,'evfile'
         open(99,FORM='UNFORMATTED',file=filename,status='old')
         read(99)((evecp(j,i),j=1,dimp),i=1,dimp)
         close(99)

         do j=1, dim
         do i=1, dimp       
          do 20 jold=1, dim
            call adag(msite,-1,nup,ndo,jold,nupnew,ndonew,jnew,isign)
            if (isign .eq. 0) then
             goto 20
            endif
            if ((ndonew .ne. ndo+1) .or. (nup .ne. nupnew)) then
             print *, 'There is problem.'
             stop
            endif
            matrix(i,j)=matrix(i,j)+real(isign)*evec(jold,j)*evecp(jnew,i)
20        continue          
         enddo
         enddo
!---- Output into file -----
         write(filename,'(a2,2I1,a6)') 'zz',nup,ndo,'ddofil'
         open(99,FORM='UNFORMATTED',file=filename,status='unknown')
         write(99)((matrix(i,j),i=1, dimp), j=1, dim)
         close(99)
       enddo
      enddo
      end subroutine dagmat
!======================================================================







!C========+=========+=========+=========+=========+=========+=========+=$
!C PROGRAM: a few subroutines from slatec
!C          which were minimally modified to compute 
!C          double precision eigensystems.
!C TYPE   : main
!C PURPOSE: 
!C I/O    :
!C VERSION: 
!C COMMENT: to learn about netlib, send an otherwise empty
!C          message to netlib@research.att.com
!C          containing 'send index' in the subject header)
!C          The WWW address of netlib is
!C          http://netlib.att.com/netlib/search.html
!Cnoprint=+=========+=========+=========+=========+=========+=========+=$
!*DECK RS
      SUBROUTINE RS (NM, N, A, W, MATZ, Z, FV1, FV2, IERR)
!C***BEGIN PROLOGUE  RS
!C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
!C            of a real symmetric matrix.
!C***LIBRARY   SLATEC (EISPACK)
!C***CATEGORY  D4A1
!C***TYPE      SINGLE PRECISION (RS-S, CH-C)
!C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!C***AUTHOR  Smith, B. T., et al.
!C***DESCRIPTION
!C
!C     This subroutine calls the recommended sequence of
!C     subroutines from the eigensystem subroutine package (EISPACK)
!C     to find the eigenvalues and eigenvectors (if desired)
!C     of a DOUBLE PRECISION SYMMETRIC matrix.
!C
!C     On Input
!C
!C        NM must be set to the row dimension of the two-dimensional
!C          array parameters, A and Z, as declared in the calling
!C          program dimension statement.  NM is an INTEGER variable.
!C
!C        N is the order of the matrix A.  N is an INTEGER variable.
!C          N must be less than or equal to NM.
!C
!C        A contains the real symmetric matrix.  A is a two-dimensional
!C          DOUBLE PRECISION array, dimensioned A(NM,N).
!C
!C        MATZ is an INTEGER variable set equal to zero if only
!C          eigenvalues are desired.  Otherwise, it is set to any
!C          non-zero integer for both eigenvalues and eigenvectors.
!C
!C     On Output
!C
!C        A is unaltered.
!C
!C        W contains the eigenvalues in ascending order.  W is a one-
!C          dimensional DOUBLE PRECISION array, dimensioned W(N).
!C
!C        Z contains the eigenvectors if MATZ is not zero.  The
!C          eigenvectors are orthonormal.  Z is a two-dimensional
!C          DOUBLE PRECISION array, dimensioned Z(NM,N).
!C
!C        IERR is an INTEGER flag set to
!C          Zero       for normal return,
!C          10*N       if N is greater than NM,
!C          J          if the J-th eigenvalue has not been
!C                     determined after 30 iterations.
!C                     The eigenvalues, and eigenvectors if requested,
!C                     should be correct for indices 1, 2, ..., IERR-1.
!C
!C  FV1 and FV2 are one-dimensional DOUBLE PRECISION arrays used for temporary
!C          storage, dimensioned FV1(N) and FV2(N).
!C
!C     Questions and comments should be directed to B. S. Garbow,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C     ------------------------------------------------------------------
!C
!C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!C                 system Routines - EISPACK Guide, Springer-Verlag,
!C                 1976.
!C***ROUTINES CALLED  TQL2, TQLRAT,  TRED2
!C***REVISION HISTORY  (YYMMDD)
!C   760101  DATE WRITTEN
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  RS
!C
!===== test to change it into Fortran90 format.=========
!      INTEGER N,NM,IERR,MATZ
!      DOUBLE PRECISION A(NM,*),W(*),Z(NM,*),FV1(*),FV2(*)      
!===========
    use prec  
    implicit none
    integer,intent(IN):: N, NM, MATZ
    integer,intent(OUT):: IERR
    real(kind=rkind),dimension(NM,N),intent(IN):: A
    real(kind=rkind),dimension(NM,N),intent(OUT)::Z
    real(kind=rkind),dimension(N),intent(OUT)::W
    real(kind=rkind),dimension(N),intent(IN)::FV1,FV2
!========================================================
!   print *, (A(3,i),i=1,8)
!   pause
!C
!C***FIRST EXECUTABLE STATEMENT  RS
      IERR = 10 * N
!C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
      CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
!      END
      end subroutine RS




!*DECK TQL2
!
      SUBROUTINE TQL2 (NM, N, D, E, Z, IERR)
!C***BEGIN PROLOGUE  TQL2
!C***PURPOSE  Compute the eigenvalues and eigenvectors of symmetric
!C            tridiagonal matrix.
!C***LIBRARY   SLATEC (EISPACK)
!C***CATEGORY  D4A5, D4C2A
!C***TYPE      SINGLE PRECISION (TQL2-S)
!C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!C***AUTHOR  Smith, B. T., et al.
!C***DESCRIPTION
!C
!C     This subroutine is a translation of the ALGOL procedure TQL2,
!C     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!C     Wilkinson.
!C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!C
!C     This subroutine finds the eigenvalues and eigenvectors
!C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
!C     The eigenvectors of a FULL SYMMETRIC matrix can also
!C     be found if  TRED2  has been used to reduce this
!C     full matrix to tridiagonal form.
!C
!C     On Input
!C
!C        NM must be set to the row dimension of the two-dimensional
!C          array parameter, Z, as declared in the calling program
!C          dimension statement.  NM is an INTEGER variable.
!C
!C        N is the order of the matrix.  N is an INTEGER variable.
!C          N must be less than or equal to NM.
! C
!C        D contains the diagonal elements of the symmetric tridiagonal
!!  matrix.  D is a one-dimensional DOUBLE PRECISION array, dimensioned D(N).
!C
!C        E contains the subdiagonal elements of the symmetric
!C          tridiagonal matrix in its last N-1 positions.  E(1) is
!C          arbitrary.  E is a one-dimensional DOUBLE PRECISION array,
!C           dimensioned
!C          E(N).
!C
!C        Z contains the transformation matrix produced in the
!C          reduction by  TRED2, if performed.  If the eigenvectors
!C          of the tridiagonal matrix are desired, Z must contain
!C          the identity matrix.  Z is a two-dimensional DOUBLE PRECISION
!C          array dimensioned Z(NM,N).
!C
!C      On Output
!C
!C        D contains the eigenvalues in ascending order.  If an
!C          error exit is made, the eigenvalues are correct but
!C          unordered for indices 1, 2, ..., IERR-1.
!C
!C        E has been destroyed.
!C
!C        Z contains orthonormal eigenvectors of the symmetric
!C          tridiagonal (or full) matrix.  If an error exit is made,
!C          Z contains the eigenvectors associated with the stored
!C          eigenvalues.
!C
!C        IERR is an INTEGER flag set to
!C          Zero       for normal return,
!C          J          if the J-th eigenvalue has not been
!C                     determined after 30 iterations.
!C
!C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!C
!C     Questions and comments should be directed to B. S. Garbow,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C     ------------------------------------------------------------------
!C
!C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!C                 system Routines - EISPACK Guide, Springer-Verlag,
!C                 1976.
!C***ROUTINES CALLED  PYTHAG
!C***REVISION HISTORY  (YYMMDD)
!C   760101  DATE WRITTEN
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  TQL2
!C
!====== test to change it to Fortran 90 format.==========
!      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
!      DOUBLE PRECISION D(*),E(*),Z(NM,*)
!    DOUBLE PRECISION B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
!      DOUBLE PRECISION PYTHAG
!========
      use prec
      implicit none
      integer,intent(IN)::NM,N
      integer,intent(INOUT)::IERR
      integer:: I,J,K,L,M,II,L1,L2,MML
      real(kind=rkind),dimension(N),intent(INOUT):: D,E
      real(kind=rkind),dimension(NM,N),intent(INOUT):: Z
      real(kind=rkind),external::PYTHAG
      real(kind=rkind):: B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
!===========================================================
!C
!C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!C
      F = 0.0_rkind
      B = 0.0_rkind
      E(N) = 0.0_rkind
!C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
!C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
!C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0_rkind * E(L))
         R = PYTHAG(P,1.0_rkind)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
!C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
!C
  145    F = F + H
!C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0_rkind
         C2 = C
         EL1 = E(L1)
         S = 0.0_rkind
         MML = M - L
!C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0_rkind)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0_rkind/ R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0_rkind)
            E(I+1) = S * E(I) * R
            S = 1.0_rkind / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!C
  200    CONTINUE
!C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
!C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
!C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
!C
  300 CONTINUE
!C
      GO TO 1001
!C     .......... SET ERROR -- NO CONVERGENCE TO AN
!C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
!
! 1001 RETURN
!      END
 1001  end subroutine TQL2




!*DECK TRED2
      SUBROUTINE TRED2 (NM, N, A, D, E, Z)
!C***BEGIN PROLOGUE  TRED2
!C***PURPOSE  Reduce a real symmetric matrix to a symmetric tridiagonal
!C            matrix using and accumulating orthogonal transformations.
!C***LIBRARY   SLATEC (EISPACK)
!C***CATEGORY  D4C1B1
!C***TYPE      SINGLE PRECISION (TRED2-S)
!C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!C***AUTHOR  Smith, B. T., et al.
!C***DESCRIPTION
!C
!C     This subroutine is a translation of the ALGOL procedure TRED2,
!C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!C
!C     This subroutine reduces a DOUBLE PRECISION SYMMETRIC matrix to a
!C     symmetric tridiagonal matrix using and accumulating
!C     orthogonal similarity transformations.
!C
!C     On Input
!C
!C        NM must be set to the row dimension of the two-dimensional
!C          array parameters, A and Z, as declared in the calling
!C          program dimension statement.  NM is an INTEGER variable.
!C
!C        N is the order of the matrix A.  N is an INTEGER variable.
!C          N must be less than or equal to NM.
!C
!C        A contains the real symmetric input matrix.  Only the lower
!C          triangle of the matrix need be supplied.  A is a two-
!C          dimensional DOUBLE PRECISION array, dimensioned A(NM,N).
!C
!C     On Output
!C
!C        D contains the diagonal elements of the symmetric tridiagonal
!C          matrix.  D is a one-dimensional DOUBLE PRECISION array,
!C          dimensioned D(N).
!C
!C        E contains the subdiagonal elements of the symmetric
!C          tridiagonal matrix in its last N-1 positions.  E(1) is set
!C          to zero.  E is a one-dimensional DOUBLE PRECISION array,
!C          dimensioned E(N).
!C
!C        Z contains the orthogonal transformation matrix produced in
!C          the reduction.  Z is a two-dimensional DOUBLE PRECISION array,
!C          dimensioned Z(NM,N).
!C
!C        A and Z may coincide.  If distinct, A is unaltered.
!C
!C     Questions and comments should be directed to B. S. Garbow,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C     ------------------------------------------------------------------
!C
!C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!C                 system Routines - EISPACK Guide, Springer-Verlag,
!C                 1976.
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   760101  DATE WRITTEN
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  TRED2
!C
!======= test to change it into Fortran 90 format.======
!      INTEGER I,J,K,L,N,II,NM,JP1
!      DOUBLE PRECISION A(NM,*),D(*),E(*),Z(NM,*)
!      DOUBLE PRECISION F,G,H,HH,SCALE
!============
      use prec
      implicit none
      integer,intent(IN)::NM,N
      integer:: I,J,K,L,II,JP1
      real(kind=rkind),dimension(NM,N),intent(IN):: A
      real(kind=rkind),dimension(NM,N),intent(OUT):: Z
      real(kind=rkind),dimension(N),intent(OUT):: D,E
      real(kind=rkind):: F,G,H,HH,SCALE
!========================================================

!C
!C***FIRST EXECUTABLE STATEMENT  TRED2
      DO 100 I = 1, N
!C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
!C
      IF (N .EQ. 1) GO TO 320
!C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0_rkind
         SCALE = 0.0_rkind
         IF (L .LT. 2) GO TO 130
!C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
!C
         IF (SCALE .NE. 0.0_rkind) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
!C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
!C
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0_rkind
!C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0_rkind
!C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
!C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
!C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
!C
         HH = F / (H + H)
!C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
!C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
!C
  290    D(I) = H
  300 CONTINUE
!C
  320 D(1) = 0.0_rkind
      E(1) = 0.0_rkind
!C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0_rkind) GO TO 380
!C
         DO 360 J = 1, L
            G = 0.0_rkind
!C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
!C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
!C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0_rkind
         IF (L .LT. 1) GO TO 500
!C
         DO 400 J = 1, L
            Z(I,J) = 0.0_rkind
            Z(J,I) = 0.0_rkind
  400    CONTINUE
!C
  500 CONTINUE
!C
!      RETURN
!      END
      end subroutine TRED2




!====== test to change it into Fortran90 format.====
!      double precision function pythag(a,b)
!      double precision a,b
!      double precision p,r,s,t,u
!==========
       function pythag(a,b)
       use prec
       implicit none
       real(kind=rkind),intent(IN):: a,b
       real(kind=rkind):: p,r,s,t,u,pythag
!====================================================
!c
!c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!c

      p = max(abs(a),abs(b))
      if (p .eq. 0.0_rkind) go to 20
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0_rkind + r
         if (t .eq. 4.0_rkind) go to 20
         s = r/t
         u = 1.0_rkind + 2.0_rkind*s
         p = u*p
         r = (s/u)**2.0_rkind * r
      go to 10
   20 pythag = p
!      return
!      end
    end function PYTHAG


