module calphy_module
     use hubbard_param
     use flags_param
     implicit none


interface cal_den
module procedure cal_den_r, cal_den_c
end interface

contains

     subroutine cal_den_r(dim, nup, ndo, Phi, denup,dendo,dbocc,kin)
     integer, intent(IN):: dim           ! dimension of the sub space 
     integer, intent(IN):: nup, ndo      ! (nup, ndo) space
     real(kind=RKind), intent(IN):: Phi(dim)
     real(kind=RKind), intent(OUT):: denup(Nsite), dendo(Nsite), dbocc, kin
     integer:: j,jnew, error
     integer,dimension(2*Nsite):: codej, codejnew

     integer:: is, ms, ns, spin
     integer(kind=4), external::array2integer

   
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
kin=0d0
!$omp parallel do private(codej,codejnew,is,jnew,sum, spin, ms, ns) reduction(+ : Phi_new)
        do  j=1,dim,1             ! sequential number of configurations in the subspace 
         call findbas(codej, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej

          do ns=1, Nsite
             denup(ns)=denup(ns)+codej(2*ns-1)* abs(phi(j))**2
             dendo(ns)=dendo(ns)+codej(2*ns)*   abs(phi(j))**2
             dbocc=dbocc+codej(2*ns-1)*codej(2*ns) *abs(phi(j))**2
          enddo 


          do spin=-1,1,2
            do  ms=1, Nsite
             do  ns=1, Nsite
             if (dabs(Tmn(ms,ns)) < 1d-10) cycle ! connect by hopping 
       
           call adaggera(ns,ms,spin,nup,ndo,codej,codejnew,is)	
           if (is==0) cycle
 
        jnew=invTable(array2integer(codejnew))           ! array rep-->demical rep-----> real inv table--->sequentail rep in subspace
        kin=kin+(phi(jnew))* phi(j) * dble(is) * (-Tmn(ms,ns))
            enddo 
          enddo 
         enddo 




        enddo 
deallocate(Table)
!print *, 'Hopt'
end subroutine cal_den_r

!========================================================================
     subroutine cal_den_c(dim, nup, ndo, Phi, denup,dendo,dbocc,kin)

     integer, intent(IN):: dim           ! dimension of the sub space 
     integer, intent(IN):: nup, ndo      ! (nup, ndo) space
     complex(kind=RKind), intent(IN):: Phi(dim)
     real(kind=RKind), intent(OUT):: denup(Nsite), dendo(Nsite), dbocc, kin
     integer:: j,jnew, error
     integer,dimension(2*Nsite):: codej, codejnew

     integer:: is, ms, ns, spin
     integer(kind=4), external::array2integer
   
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
kin=0d0
!$omp parallel do private(codej,codejnew,is,jnew,sum, spin, ms, ns) reduction(+ : Phi_new)
        do  j=1,dim,1             ! sequential number of configurations in the subspace 
         call findbas(codej, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej

          do ns=1, Nsite
             denup(ns)=denup(ns)+codej(2*ns-1)* abs(phi(j))**2
             dendo(ns)=dendo(ns)+codej(2*ns)*   abs(phi(j))**2
             dbocc=dbocc+codej(2*ns-1)*codej(2*ns) *abs(phi(j))**2
          enddo 

          do spin=-1,1,2
            do  ms=1, Nsite
             do  ns=1, Nsite
             if (dabs(Tmn(ms,ns)) < 1d-10) cycle ! connect by hopping 
       
           call adaggera(ns,ms,spin,nup,ndo,codej,codejnew,is)	
           if (is==0) cycle
 
        jnew=invTable(array2integer(codejnew))           ! array rep-->demical rep-----> real inv table--->sequentail rep in subspace
        kin=kin+conjg(phi(jnew))* phi(j) * dble(is) * (- Tmn(ms,ns) )
            enddo 
          enddo 
         enddo 


        enddo 
deallocate(Table)
!print *, 'Hopt'
end subroutine cal_den_c


!========================================================================
     subroutine cal_den_mix(dim, nup, ndo,phil,phir,Occ,Dbocc,overlap,kin)
     implicit none

     integer, intent(IN):: dim           ! dimension of the sub space 
     integer, intent(IN):: nup, ndo      ! (nup, ndo) space
     complex(kind=8), intent(IN):: phil(dim), phir(dim)
     integer:: j,jnew, error
     integer,dimension(2*Nsite):: codej, codejnew

     complex(kind=8), intent(OUT):: Occ(Nsite,2),Dbocc

     complex(kind=8),intent(OUT):: overlap,kin
     integer:: is, spin, ms, ns
     integer(kind=4), external::array2integer


error=0
allocate(Table(Nstate(nup, ndo)), stat=error)
if(error/=0) then
print *, 'wrong in allocate Table, cal_den_mix', error
stop
endif


call buildbasis(nup, ndo) ! build Table and invTable 


Occ=0d0
Dbocc=0d0
overlap=0d0
kin=0d0
!$omp parallel do private(codej,codejnew,is,jnew,sum, spin, ms, ns) reduction(+ : Phi_new)
        do  j=1,dim,1             ! sequential number of configurations in the subspace 
         call findbas(codej, nup, ndo, j)      ! j--> read table---> demical rep--->array rep codej

          do ns=1, Nsite
        Occ(ns,1)=Occ(ns,1)+codej(2*ns-1)* conjg(phil(j)) * phir(j)
        Occ(ns,2)=Occ(ns,2)+codej(2*ns)*   conjg(phil(j)) * phir(j)
        dbocc=dbocc+codej(2*ns-1)*codej(2*ns) * conjg(phil(j)) * phir(j)
          enddo 

        overlap=overlap+  conjg(phil(j)) * phir(j)


        do spin=-1,1,2
           do  ms=1, Nsite
               do    ns=1, Nsite
             if (dabs(Tmn(ms,ns)) < 1d-10) cycle  ! connect by hopping 
       
       call adaggera(ns,ms,spin,nup,ndo,codej,codejnew,is)
        if (is==0) cycle

        jnew=invTable(array2integer(codejnew))           ! array rep-->demical rep-----> real inv table--->sequentail rep in subspace
        kin=kin+conjg(phil(jnew))* phir(j) * dble(is) * (- Tmn(ms,ns) )

           enddo 
          enddo 
         enddo 
        enddo 

       occ=occ/overlap
       dbocc=dbocc/overlap
       kin=kin/overlap

deallocate(Table)
!print *, 'Hopt'
end subroutine cal_den_mix


end module calphy_module
