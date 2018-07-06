!=============================================================================
subroutine sub_lanc(nup, ndo,dim,sub_en, phitemp)
use prec
use lanc_param
use hubbard_param
use flags_param
use Hopt_module
!use omp_lib

implicit none

     integer,intent(IN):: nup, ndo      ! (nup, ndo) space
     integer,intent(IN):: dim 
     real(kind=Rkind), intent(OUT):: sub_en(Noflow)           ! gs energy in this subspace 
     real(kind=rkind), intent(OUT):: phitemp(dim)


real (kind=RKind), allocatable:: Phi(:, :)   !Phi(# of bisis for LM, # of basis in Hilbert space)
real (kind=RKind):: Have(0: Lanc_size-1), Bb(0:Lanc_size-1) ! elements of the Lanczos Matrix

real (kind=RKind):: Hlanc(Lanc_size,Lanc_size),energy(Lanc_size) !eigenvalues of the LM
real (kind=RKind):: eigenvector(Lanc_size,Lanc_size) !eigenvectors of the LM 
real (kind=RKind)::entemp(Noflow)               ! used to compare with Eg, serve as the criterion of convergence
real (kind=RKind):: Ortho(0:Lanc_size-2) ! for Schimit orthogonalization

real (kind=RKind):: sum, res
integer :: n,m,k,ii,j, iter, l
integer:: Matsize  ! always equal to Lanc_size
real(kind=RKind), external:: random
integer:: iseed


matsize=Lanc_size
!-------------------------
iseed=12345
entemp=1d2
!---------------------------


allocate(Phi(0:matsize-1 ,dim))
Phi=0.0_RKInd


!----------------random the initial state-------------------------------
res=0.0_Rkind
!$omp parallel do if(dim>Condition_omp) reduction(+:res)
do k=1, dim,1
phi(0,k)=2.0_RKind*(random(iseed)-0.5_RKind)
res=res+phi(0,k)*phi(0,k)
enddo

!$omp parallel do if(dim>Condition_omp)
do k=1, dim,1
Phi(0,k)=Phi(0,k)/dsqrt(res)                                !normalize |Phi_0>
enddo


iter=0
10   continue                                 ! lanczos iteration 
    iter=iter+1
   
	Have=0.0_RKind
	Bb=0.0_Rkind  !initializtion


do n=0, matsize-1
 if(Iterprint==1) print *, 'making lanc_vec', iter, n

 call Hopt(dim, nup, ndo, Phi(n,:), Phitemp(:))    !  H|Phi_n>гнгнгн>|Phi> temp

!--------------H_ave-----------------------------------------
res=0.0_RKind
!$omp parallel do if(dim>Condition_omp) reduction(+:res)  
  do k=1, dim,1
   res=res+Phitemp(k)*Phi(n,k)
  enddo
Have(n)=res
  
!---------------------above is excute even when n=matsize-1,----------------------------
!-----------------------when n=matsize-1, below is not ex ---------------------------
if (n<Matsize-1) then  ! not the last one so there is Bb to built
   
    
	   ortho=0.0_Rkind
!----------|Phi_(n+1)>-----------------------------------------
     if (n.eq.0) then
	 
!$omp parallel do reduction(+: ortho)  
	   do k=1, dim, 1
	      Phi(n+1,k)=Phitemp(k)-Have(n)*Phi(n,k)                    !|Phi_1>=H|Phi_0>-a_0|Phi_0>
	      do m=0, n
		  ortho(m)=ortho(m)+(Phi(n+1, k)*Phi(m,k))
		  enddo
	   enddo

	 
	 else
	 
!$omp parallel do reduction(+:ortho)
          do k=1, dim, 1	   
           Phi(n+1,k)=Phitemp(k)-Have(n)*Phi(n,k)-Bb(n)*Phi(n-1,k)  ! |Phi_n+1>=H|Phi_n>-a_n|Phi_n>-b_n|Phi_n-1>
            do m=0, n
		    ortho(m)=ortho(m)+(Phi(n+1, k)*Phi(m,k))
			enddo
		  enddo   

           
      endif	  
 
!--- Schmit orthogonalization---------------------------------      
	   
            res=0.0_RKind
		 !$omp parallel do reduction(+: res)
	      do k=1, dim,1
		    do m=0,n
		     Phi(n+1,k)=Phi(n+1,k)-Phi(m,k)*ortho(m)
            enddo
		     res=res+Phi(n+1,k)*Phi(n+1,k)
          enddo

    !    print *,iter, n+1, dot_product(PHI(n+1,:), PHi(n-1, :))
	!	  pause



!--------------------------------------------------------------
!-- calculate Bb(n+1) and then normalize |phi_n+1> --------
      !res=dot_product(Phi(n+1, :), Phi(n+1, :))
	  !Bb(n+1)=dsqrt(res)  	  
      !Phi(n+1, :)= Phi(n+1, :)/Bb(n+1)
if(res<1d-10) then 
Matsize=n+1
exit
endif
Bb(n+1)=dsqrt(res)

!$omp parallel do if(dim>Condition_omp)
do k=1, dim,1
Phi(n+1,k)=Phi(n+1,k)/Bb(n+1)                                !normalize |Phi_n+1>
enddo


endif ! end produce Lanczos vectors

enddo  ! end produce Have and Bbs 

!------check the otho------------------------------------------
!j=matsize-1
!do m=0, j
!print *, dot_product(PHi(j,:), PHi(m,:))
!enddo
!pause
!------------------------------------------------------
!print *, 'have' ,have
!print *, 'bb', bb

!--- produce tri-diagonal hamilt matrix ------
!    Hlanc=0.0_Rkind
!     do n=0, Matsize-2
!	    Hlanc(n+1,n+1)=Have(n)
!		Hlanc(n+1,n+2)=Bb(n+1)
!        Hlanc(n+2,n+1)=Hlanc(n+1,n+2)
!     enddo
!	 n=Matsize-1
!     Hlanc(n+1,n+1)=Have(n)
!----------------------------------------------
!call diaMat_realsymm(matsize,Hlanc, energy, eigenvector) !rs routine , eigenvalue from small to large 

call diaMat_tri(matsize,Have, Bb, energy, eigenvector)   ! upon call, Have also gives energy and Bbtemp has been destroyed!



do ii=1, Noflow
     sub_en(ii)=energy(ii) 
enddo

!--------update the initial state |phi_0>---------------- 

!$omp parallel do private(sum, l)   
    do k=1, dim, 1
	   sum=0.0_Rkind
	   do l=1, Matsize
	     sum=sum+eigenvector(l,1)*Phi(l-1,k)                   ! coefficient on working basis
	   enddo
	  Phitemp(k)=sum                                            ! update |Phi_0>
	enddo
!-------------------------------------------------------------------


!-------------------converge or not?---------------------------------
    res=0.0_RKind
	 do ii=1, Noflow
      res=res+dabs(entemp(ii)-sub_en(ii))
	  enddo    

if (iterprint==1) print *, 'iter in lanc', iter, res
 
	 
	    if (res > Pre) then
                     phi(0,k)=phitemp(k)
		     entemp=sub_en                   ! put the current Energy into entemp(j) for the comparation after the next iteration 
		   goto 10                                               !not converge, go on iteration
        endif
!-----------------------------------------------------------------------------


deallocate(Phi)

return

end subroutine sub_lanc


!=============================================================================
subroutine sub_lancdynamics(nup, ndo, dim, phi0,matsize,Hlanc,Phi)
use prec
use lanc_param
use hubbard_param
use flags_param
use Hopt_module
!use omp_lib

implicit none

     integer,intent(IN):: nup, ndo      ! (nup, ndo) space
     integer,intent(IN):: dim           
    
complex (kind=RKind),intent(IN):: Phi0(dim)   
complex (kind=RKind),intent(OUT):: Phi(0:Lanc_size-1,dim)  !Phi(# of bisis for LM, # of basis in Hilbert space)
real (kind=RKind),intent(OUT):: Hlanc(Lanc_size,Lanc_size) 


real (kind=RKind):: Have(0: Lanc_size-1), Bb(0:Lanc_size-1) ! elements of the Lanczos Matrix
complex (kind=RKind):: Ortho(0:Lanc_size-2) ! for Schimit orthogonalization
complex (kind=RKind), allocatable:: Phitemp(:)    !Phi(# of bisis for LM, # of basis in Hilbert space)

real (kind=RKind):: sum, res
integer :: n,m,k,ii,j, iter, l
integer:: Matsize  ! always equal to Lanc_size

matsize=Lanc_size


allocate(Phitemp(dim))
Phi=0.0_RKInd


!$omp parallel do if(dim>Condition_omp)
do k=1, dim,1
Phi(0,k)=Phi0(k)                              !normalize |Phi_0>
enddo



iter=0
10   continue                                 ! lanczos iteration 
    iter=iter+1
   
	Have=0.0_RKind
	Bb=0.0_Rkind  !initializtion


do n=0, matsize-1
 if(Iterprint==1) print *, 'making lanc_vec', iter, n

 call Hopt(dim, nup, ndo, Phi(n,:), Phitemp(:))  

!--------------H_ave-----------------------------------------
res=0.0_RKind
!$omp parallel do if(dim>Condition_omp) reduction(+:res)  
  do k=1, dim,1
   res=res+Phitemp(k)* conjg(Phi(n,k))
  enddo
Have(n)=res
  
!---------------------above is excute even when n=matsize-1,----------------------------
!-----------------------when n=matsize-1, below is not ex ---------------------------
if (n<Matsize-1) then  ! not the last one so there is Bb to built
   
    
	   ortho=0.0_Rkind
!----------|Phi_(n+1)>-----------------------------------------
     if (n.eq.0) then
	 
!$omp parallel do reduction(+: ortho)  
	   do k=1, dim, 1
	      Phi(n+1,k)=Phitemp(k)-Have(n)*Phi(n,k)                    !|Phi_1>=H|Phi_0>-a_0|Phi_0>
	      do m=0, n
		    ortho(m)=ortho(m)+Phi(n+1, k)*conjg(Phi(m,k))
		  enddo
	   enddo

	 
	 else
	 
!$omp parallel do reduction(+:ortho)
          do k=1, dim, 1	   
           Phi(n+1,k)=Phitemp(k)-Have(n)*Phi(n,k)-Bb(n)*Phi(n-1,k)  ! |Phi_n+1>=H|Phi_n>-a_n|Phi_n>-b_n|Phi_n-1>
            do m=0, n
		    ortho(m)=ortho(m)+Phi(n+1, k)*conjg(Phi(m,k))
			enddo
		  enddo   

           
      endif	  
 
!--- Schmit orthogonalization---------------------------------      
	   
            res=0.0_RKind
		 !$omp parallel do reduction(+: res)
	      do k=1, dim,1
		    do m=0,n
		     Phi(n+1,k)=Phi(n+1,k)-Phi(m,k)*ortho(m)
            enddo
		     res=res+abs(Phi(n+1,k))**2
          enddo

    !    print *,iter, n+1, dot_product(PHI(n+1,:), PHi(n-1, :))
	!	  pause



!--------------------------------------------------------------
!-- calculate Bb(n+1) and then normalize |phi_n+1> --------
      !res=dot_product(Phi(n+1, :), Phi(n+1, :))
	  !Bb(n+1)=dsqrt(res)  	  
      !Phi(n+1, :)= Phi(n+1, :)/Bb(n+1)

if(res<1d-10) then 
Matsize=n+1
exit
endif

Bb(n+1)=dsqrt(res)

!$omp parallel do if(dim>Condition_omp)
do k=1, dim,1
Phi(n+1,k)=Phi(n+1,k)/Bb(n+1)                             !normalize |Phi_n+1>
enddo


endif ! end produce Lanczos vectors

enddo  ! end produce Have and Bbs 

!------check the otho------------------------------------------
!j=matsize-1
!do m=0, j
!print *, dot_product(PHi(j,:), PHi(m,:))
!enddo
!pause
!------------------------------------------------------
!--- produce tri-diagonal hamilt matrix ------
    Hlanc=0.0_Rkind
     do n=0, Matsize-2
        Hlanc(n+1,n+1)=Have(n)
	Hlanc(n+1,n+2)=Bb(n+1)
        Hlanc(n+2,n+1)=Hlanc(n+1,n+2)
     enddo
	 n=Matsize-1
     Hlanc(n+1,n+1)=Have(n)
!----------------------------------------------
!call diaMat_realsymm(matsize,Hlanc, energy, eigenvector) !rs routine , eigenvalue from small to large 

deallocate(Phitemp)

return

end subroutine sub_lancdynamics


