!subroutine propagationHF_onestep(expKHF) ! propagation HF state one step
!use system_parameters
!use hubbard_param
!implicit none

!complex(kind=8), intent(IN):: expKHF(Nsite,Nsite,2)
!integer:: spin

!   do spin=1,2
!   phiT(1:Nsite,1:numupdo(spin),spin)=matmul(expKHF(1:Nsite,1:Nsite,spin),phiT(1:Nsite,1:numupdo(spin),spin)) 
!   enddo 

!end subroutine propagationHF_onestep


subroutine propagation_onestep1(step,direction,xbar,x2bar,expfactor,fit_error) ! propagation one step
use system_parameters
use hubbard_param
implicit none

integer, intent(IN):: step, direction
complex(kind=8), intent(OUT):: xbar(Nsite,Nwalkers), expfactor(Nwalkers)
complex(kind=8):: x2bar(Nsite,Nwalkers)

complex(kind=8):: expV
integer:: spin,i,site

real(kind=8):: ran
complex(kind=8):: Auxfield ! is complex due to the shift 
real(kind=8),external:: gasdev, rannyu
!integer:: ising

complex(kind=8)::overlap1(2), overlap2(2),overlap12(2)
complex(kind=8)::phi_1(Nsite,max(nup,ndo),2), phi_2(Nsite,max(nup,ndo),2)
!real(kind=8):: offset

real(kind=8):: offset1, offset2

real(kind=8):: fit_error



expfactor=0d0
xbar=0d0
x2bar=0d0


fit_error=0d0
    do  i=1,Nwalkers

          call get_xbar(i,direction, xbar(:,i), x2bar(:,i))
            if ((abs(xbar(1,i))-Maxshift)>1d-8 ) then
                print *, 'wrong xshift', abs(xbar(1,i))
                stop
             endif 

         !-----------exp(-i*Kt/2)-----------------------
         do spin=1,2
         phi(1:Nsite,1:numupdo(spin),spin,i)=matmul(expK(1:Nsite,1:Nsite),phi(1:Nsite,1:numupdo(spin),spin,i)) 
         enddo 


         ! copy the state
         !-----------exp(-iVt)-----------------------
         do site=1,Nsite

         ran=gasdev()
         Auxfield=ran - xbar(site,i)

          if (Local_energy) then 
         expfactor(i)=expfactor(i) + (x2bar(site,i))
          else
         expfactor(i)=expfactor(i) + (ran*xbar(site,i)-0.5d0*(xbar(site,i))**2)
          endif

         expV=exp(gamma*Auxfield) 
         phi(site,1:nup,1,i)=expV*phi(site,1:nup,1,i)
          
         expV=exp(-gamma*Auxfield) 
         phi(site,1:ndo,2,i)=expV*phi(site,1:ndo,2,i)
         enddo
 

         !-----------exp(-i*Kt/2)-----------------------

         if (TrotterOrder==2) then 

         do spin=1,2
         phi(1:Nsite,1:numupdo(spin),spin,i)=matmul(expK(1:Nsite,1:Nsite),phi(1:Nsite,1:numupdo(spin),spin,i)) 
         enddo 


!         print *, 'K/2 2 passed' 
         endif

!       print *, 'xbar', xbar(:,i)
!       print *, 'expfactor', expfactor(i)
!       pause

   enddo !walker

!fit_error=fit_error/dble(Nwalkers)

!write (56,*) step, fit_error

end subroutine propagation_onestep1



!function phaseless_proj(z)
!implicit none
!complex(kind=8),intent(IN)::z
!real(kind=8)::phaseless_proj

!complex(kind=8):: lnz
!complex(kind=8),external:: clog

!lnz=clog(z)
!phaseless_proj=max(cos(dimag(lnz)),0d0)

!if (dble(z)<0d0) then 
!print *, 'dble(z)', dble(z)
!pause
!endif 

!end function  phaseless_proj



