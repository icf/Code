module calphy_module
     use hubbard_param
     use flags_param
     implicit none


interface cal_den
module procedure cal_den_r, cal_den_c
end interface

contains

     subroutine cal_den_r(dim, nup, ndo, Phi, denup,dendo,dbocc)
     integer, intent(IN):: dim           ! dimension of the sub space 
     integer, intent(IN):: nup, ndo      ! (nup, ndo) space
     real(kind=RKind), intent(IN):: Phi(dim)
     real(kind=RKind), intent(OUT):: denup(Nsite), dendo(Nsite), dbocc
     integer:: j, ns, error
     integer,dimension(2*Nsite):: codej

   
error=0
allocate(Table(Nstate(nup, ndo)), stat=error)
if(error/=0) then
print *, 'wrong in allocate Table'
stop
endif


call buildbasis(nup, ndo) ! build Table and invTable 

denup=0d0
dendo=0d0
dbocc=0d0
!$omp parallel do private(codej,codejnew,is,jnew,sum, spin, ms, ns) reduction(+ : Phi_new)
        do  j=1,dim,1             ! sequential number of configurations in the subspace 
         call findbas(codej, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej

          do ns=1, Nsite
             denup(ns)=denup(ns)+codej(2*ns-1)* abs(phi(j))**2
             dendo(ns)=dendo(ns)+codej(2*ns)*   abs(phi(j))**2
             dbocc=dbocc+codej(2*ns-1)*codej(2*ns) *abs(phi(j))**2
          enddo 

        enddo 
deallocate(Table)
!print *, 'Hopt'
end subroutine cal_den_r

!========================================================================
     subroutine cal_den_c(dim, nup, ndo, Phi, denup,dendo,dbocc)

     integer, intent(IN):: dim           ! dimension of the sub space 
     integer, intent(IN):: nup, ndo      ! (nup, ndo) space
     complex(kind=RKind), intent(IN):: Phi(dim)
     real(kind=RKind), intent(OUT):: denup(Nsite), dendo(Nsite), dbocc
     integer:: j, ns, error
     integer,dimension(2*Nsite):: codej

   
error=0
allocate(Table(Nstate(nup, ndo)), stat=error)
if(error/=0) then
print *, 'wrong in allocate Table'
stop
endif


call buildbasis(nup, ndo) ! build Table and invTable 


denup=0d0
dendo=0d0
dbocc=0d0
!$omp parallel do private(codej,codejnew,is,jnew,sum, spin, ms, ns) reduction(+ : Phi_new)
        do  j=1,dim,1             ! sequential number of configurations in the subspace 
         call findbas(codej, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej

          do ns=1, Nsite
             denup(ns)=denup(ns)+codej(2*ns-1)* abs(phi(j))**2
             dendo(ns)=dendo(ns)+codej(2*ns)*   abs(phi(j))**2
             dbocc=dbocc+codej(2*ns-1)*codej(2*ns) *abs(phi(j))**2
          enddo 

        enddo 
deallocate(Table)
!print *, 'Hopt'
end subroutine cal_den_c

!========================================================================
     subroutine cal_den_mix(dim, nup, ndo,phil,phir,Occ,Dbocc)
     implicit none

     integer, intent(IN):: dim           ! dimension of the sub space 
     integer, intent(IN):: nup, ndo      ! (nup, ndo) space
     complex(kind=RKind), intent(IN):: phil(dim), phir(dim)
     integer:: j, ns, error
     integer,dimension(2*Nsite):: codej

     real(kind=RKind), intent(OUT):: Occ(Nsite,2),Dbocc

     complex(kind=8):: Occ_dummy(Nsite,2),Dbocc_dummy
     complex(kind=8):: overlap


error=0
allocate(Table(Nstate(nup, ndo)), stat=error)
if(error/=0) then
print *, 'wrong in allocate Table, cal_den_mix', error
stop
endif


call buildbasis(nup, ndo) ! build Table and invTable 


Occ_dummy=0d0
Dbocc_dummy=0d0
overlap=0d0
!$omp parallel do private(codej,codejnew,is,jnew,sum, spin, ms, ns) reduction(+ : Phi_new)
        do  j=1,dim,1             ! sequential number of configurations in the subspace 
         call findbas(codej, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej

          do ns=1, Nsite
        Occ_dummy(ns,1)=Occ_dummy(ns,1)+codej(2*ns-1)* conjg(phil(j)) * phir(j)
        Occ_dummy(ns,2)=Occ_dummy(ns,2)+codej(2*ns)*   conjg(phil(j)) * phir(j)
        dbocc_dummy=dbocc_dummy+codej(2*ns-1)*codej(2*ns) * conjg(phil(j)) * phir(j)
          enddo 

        overlap=overlap+  conjg(phil(j)) * phir(j)
        enddo 

       occ=occ_dummy/overlap
       dbocc=dbocc_dummy/overlap

deallocate(Table)
!print *, 'Hopt'
end subroutine cal_den_mix


end module calphy_module
