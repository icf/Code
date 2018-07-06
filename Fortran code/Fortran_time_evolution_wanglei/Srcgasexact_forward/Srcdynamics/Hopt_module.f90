!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: Hopt
! TYPE   : subroutine
! PURPOSE: compute H |phi> in subspace given by nup, ndo 
!       
! I/O    : 
! VERSION: 26-Dec-2007
! COMMENT:
! AUTHOR:  Lei Wang, IOP, CAS
!========+=========+=========+=========+=========+=========+=========+=$ 
module Hopt_module
     use hubbard_param
     use flags_param
     implicit none

interface Hopt
  module procedure Hopt_r, Hopt_c
end interface Hopt

contains

     subroutine Hopt_r(dim, nup, ndo, Phi_old, Phi_new)

     integer, intent(IN):: nup, ndo      ! (nup, ndo) space
	 integer, intent(IN):: dim             ! dimension of the sub space 
     real(kind=RKind), intent(IN):: Phi_old(dim)
	 real(kind=RKind), intent(OUT):: Phi_new(dim)      
     	  
     integer:: i,ms, ns, spin, is
	 integer:: j, jnew, jtemp
     integer:: sum,m,n
	 integer,dimension(2*Nsite):: codej, codejnew
	 integer(kind=4), external::array2integer, lookup
	 integer::error

   
error=0
allocate(Table(Nstate(nup, ndo)), stat=error)
if(error/=0) then
print *, 'wrong in allocate Table'
stop
endif



call buildbasis(nup, ndo) ! build Table and invTable 



!----------clear-------------------------------
Phi_new=0.0_RKind
!-------- construct hamiltonian matrix ---------

!$omp parallel do private(codej,codejnew,is,jnew,sum, spin, ms, ns) reduction(+ : Phi_new)
        do  j=1,dim,1             ! sequential number of configurations in the subspace 
         call findbas(codej, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej

!---- hopping term ----
          do spin=-1,1,2
           do  ms=1, Nsite
		   do    ns=1, Nsite
			 if (dabs(Tmn(ms,ns)) .gt. 1d-10) then  ! connect by hopping 
 	      
 	      call adaggera(ns,ms,spin,nup,ndo,codej,codejnew,is)	
 	       if (is==0) cycle
 	     
 	     if(invTable_flag==1) then 
		   jnew=invTable(array2integer(codejnew))           ! array rep-->demical rep-----> real inv table--->sequentail rep in subspace
		                                                    ! it is memerory consuming, could be replace by lookup fucntion, but time consuming 
		  else
		      jnew=lookup(codejnew, nup,ndo)    
		 endif
		 ! print *,'===========',array2integer(codejnew),jnew
		  
        
          Phi_new(jnew)=Phi_new(jnew)+dble(is)*(-Tmn(ms,ns))* Phi_old(j)
             endif

enddo !ns
 enddo ! ms
   enddo ! spin 
!---- U term ----  
	   sum=0
       do ms=1, Nsite
	      sum=sum+codej(2*ms-1)*codej(2*ms)      ! doubly occupied contribute energy U 
       enddo
	   Phi_new(j)=Phi_new(j)+U*dble(sum)*Phi_old(j)
!---- Mu term -----
       sum=0
       do ms=1, Nsite
       sum=sum+codej(2*ms-1)+codej(2*ms)
       enddo
      Phi_new(j)=Phi_new(j)-Mu*dble(sum)*Phi_old(j)
!---- External potential term -----
      do ns=1,Nsite
      Phi_new(j)=Phi_new(j)+ Vext*(dble(ns)-dble(Nsite/2))* dble(codej(2*ns-1)-codej(2*ns))*  Phi_old(j)
      enddo 
!---- Staggered magnetic field -----
      do ns=1,Nsite
      Phi_new(j)=Phi_new(j)+ (-1d0)**ns*stagger_h* dble(codej(2*ns-1)-codej(2*ns))*  Phi_old(j)
      enddo 


enddo !j




!----------- Phi_new complete!------------------
deallocate(Table)
!print *, 'Hopt'

end subroutine Hopt_r
!========================================================================

     subroutine Hopt_c(dim, nup, ndo, Phi_old, Phi_new)


     integer, intent(IN):: nup, ndo      ! (nup, ndo) space
	 integer, intent(IN):: dim             ! dimension of the sub space 
     complex(kind=RKind), intent(IN):: Phi_old(dim)
     complex(kind=RKind), intent(OUT):: Phi_new(dim)      
     	  
     integer:: i,ms, ns, spin, is
	 integer:: j, jnew, jtemp
     integer:: sum,m,n
	 integer,dimension(2*Nsite):: codej, codejnew
	 integer(kind=4), external::array2integer, lookup
	 integer::error

   
error=0
allocate(Table(Nstate(nup, ndo)), stat=error)
if(error/=0) then
print *, 'wrong in allocate Table'
stop
endif



call buildbasis(nup, ndo) ! build Table and invTable 



!----------clear-------------------------------
Phi_new=0.0_RKind
!-------- construct hamiltonian matrix ---------

!$omp parallel do private(codej,codejnew,is,jnew,sum, spin, ms, ns) reduction(+ : Phi_new)
        do  j=1,dim,1             ! sequential number of configurations in the subspace 
         call findbas(codej, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej

!---- hopping term ----
          do spin=-1,1,2
           do  ms=1, Nsite
		   do    ns=1, Nsite
			 if (dabs(Tmn(ms,ns)) .gt. 1d-10) then  ! connect by hopping 
 	      
 	      call adaggera(ns,ms,spin,nup,ndo,codej,codejnew,is)	
 	       if (is==0) cycle
 	     
 	     if(invTable_flag==1) then 
		   jnew=invTable(array2integer(codejnew))           ! array rep-->demical rep-----> real inv table--->sequentail rep in subspace
		                                                    ! it is memerory consuming, could be replace by lookup fucntion, but time consuming 
		  else
		      jnew=lookup(codejnew, nup,ndo)    
		 endif
		 ! print *,'===========',array2integer(codejnew),jnew
		  
        
          Phi_new(jnew)=Phi_new(jnew)+dble(is)*(-Tmn(ms,ns))* Phi_old(j)
             endif

enddo !ns
 enddo ! ms
   enddo ! spin 
!---- U term ----  
	   sum=0
       do ms=1, Nsite
	      sum=sum+codej(2*ms-1)*codej(2*ms)      ! doubly occupied contribute energy U 
       enddo
	   Phi_new(j)=Phi_new(j)+U*dble(sum)*Phi_old(j)
!---- Mu term -----
       sum=0
       do ms=1, Nsite
       sum=sum+codej(2*ms-1)+codej(2*ms)
       enddo
      Phi_new(j)=Phi_new(j)-Mu*dble(sum)*Phi_old(j)
!---- External potential term -----
      do ns=1,Nsite
      Phi_new(j)=Phi_new(j)- Vext*(dble(ns)-dble(Nsite/2))* dble(codej(2*ns-1)-codej(2*ns))*  Phi_old(j)
      enddo 
!---- Staggered magnetic field -----
      do ns=1,Nsite
      Phi_new(j)=Phi_new(j)+ (-1d0)**ns*stagger_h* dble(codej(2*ns-1)-codej(2*ns))*  Phi_old(j)
      enddo 

enddo !j


!----------- Phi_new complete!------------------
deallocate(Table)
!print *, 'Hopt'

end subroutine Hopt_c

end module Hopt_module
