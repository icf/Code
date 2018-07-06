! translater a slater wave function to the fock basis
subroutine slater2fock(nup,ndo,dim,phi_slater_up,phi_slater_do,phi_fock)
use hubbard_param
implicit none

integer,intent(IN):: nup,ndo,dim
complex(kind=8),intent(IN):: phi_slater_up(Nsite,nup)
complex(kind=8),intent(IN):: phi_slater_do(Nsite,ndo)

complex(kind=8),intent(OUT):: phi_fock(dim)

complex(kind=8):: cdummy
integer:: num_perm_up, num_perm_do
integer:: perm0_up(nup),perm0_do(ndo)
integer:: perm_up(nup),perm_do(ndo)
integer:: npar,site,i,j,l,m
integer:: perm_sign , perm_sign_up,perm_sign_do
integer:: codej(2*Nsite), posup(nup), posdo(ndo), posupdo(nup+ndo)
!integer:: posuptemp(nup), posdotemp(ndo)
integer:: error

!real(kind=8):: norm

!print *,nup,ndo,dim

!print *, 'phiup', phi_slater_up
!print *, 'phido',phi_slater_do
!pause


error=0
allocate(Table(Nstate(nup, ndo)), stat=error)
if(error/=0) then
print *, 'wrong in allocate Table'
stop
endif


call buildbasis(nup, ndo) ! build Table and invTable 



!----initialize---------------
do i=1,nup
perm0_up(i)=i
enddo 

do i=1,ndo
perm0_do(i)=i
enddo 

num_perm_up=product(perm0_up(:))
num_perm_do=product(perm0_do(:))

perm_up=perm0_up
perm_do=perm0_do


phi_fock=0d0
do j=1,dim ! loop over configuration

  call findbas(codej, nup, ndo, j)  ! j--> read table---> demical rep--->array rep codej

  ! Get postion of particles 
     npar=1
  do site=1, Nsite
   if (codej(2*site-1)==1) then 
     posup(npar) = site
     npar=npar+1
   endif
  enddo 

     npar=1
  do site=1, Nsite
   if (codej(2*site)==1) then 
     posdo(npar) = site
     npar=npar+1
   endif
  enddo 

 ! for each permutation, multiply of phi_slater gives the wave function
 perm_sign_up=1
 do l=1,num_perm_up    
 cdummy=1d0

!  print *, 'j, codej', j,',',codej
!  print *, 'posup',posup
!  print *, 'posdo',posdo
!  print *, 'perm0_up', perm0_up
!  print *, 'perm_up', perm_up
 
! print *, 'codej', codej
 do npar=1,nup 
 cdummy=cdummy* phi_slater_up(posup(perm_up(npar)),npar) 
! print *,'multiply',(posup(perm_up(npar))),npar, phi_slater_up(posup(perm_up(npar)),npar)
 enddo 

 cdummy=perm_sign_up*cdummy

! print *, 'posup,sign', posup, ',', perm_sign_up

 call perm_cy(l,nup, perm0_up, perm_up) 

! posuptemp=posup
! do npar=1,nup
! posup(npar)=posuptemp(perm_up(npar))
! enddo 

 call find_parity(nup,perm_up,perm_sign_up)

 perm_sign_do=1
 do m=1,num_perm_do

! print *, 'l,m' ,l, m
! print *, 'posdo,sign', posdo,',', perm_sign_do

 do npar=1,ndo
 cdummy=cdummy* phi_slater_do(posdo(perm_do(npar)),npar) 
 enddo 

 cdummy=perm_sign_do*cdummy

 call perm_cy(m,ndo, perm0_do, perm_do) 

! posdotemp=posdo
! do npar=1,ndo
! posdo(npar)=posdotemp(perm_do(npar))
! enddo 

 call find_parity(ndo,perm_do,perm_sign_do)


 ! combine spin up and spin down, find the sign due to the combining 
  do npar=1,nup
   posupdo(npar)=2*posup(npar)-1
  enddo 
  do npar=1,ndo
   posupdo(nup+npar)=2*posdo(npar)
  enddo 

  call find_parity(nup+ndo,posupdo,perm_sign)
!  print *, 'posupdo,perm_sign', posupdo, perm_sign

  phi_fock(j)=phi_fock(j)+cdummy*perm_sign

! print *, 'posupdo,sign', posupdo, perm_sign
 !print *, 'cdummy', cdummy
! pause

   enddo  ! m 
 enddo  ! l

 enddo 

!print *, 'phi_slater_up',phi_slater_up
!print *, 'phi_slater_do',phi_slater_do

!print *, 'phi_fock'
!do j=1,dim

! call findbas(codej, nup, ndo, j)  ! j--> read table---> demical rep--->array rep codej
!print *, phi_fock(j)
!print *, codej(:) 
!enddo 
!pause

!norm=0d0
!do j=1,dim
!norm=norm+abs(phi_fock(j))**2
!enddo 
!print *, 'norm', norm

deallocate(Table)

end subroutine slater2fock
