subroutine xhf()
use lattice_param
use xhf_param
use model_param
implicit none
integer::i,j,k

!allocate the array
 call allocate_xhf_array()

 open(unit=10,file='data.dat',status='replace')
 dmax=dmax+lit
 do delt=dmin,dmax,dstep
    call set_h0xhf()

    call input_phiT(2*Nsite,Ntot,h0xhf,ph)

    call esitmate()
    
    write(10,*) delt,dble(e_var)
    if(dble(e_var).LT.E_save) then
      E_save=dble(e_var)
      ph_save=ph
    end if
 end do
 write(10,*) "min:",E_save
 close(10)

 open(unit=10,file='phiT.dat',status='replace')
 do k=1,Ntot,1
    do j=1,2*Nsite,1
       write(10,*) ph_save(j,k)
    end do
 end do
 close(10)
end subroutine xhf


subroutine set_h0xhf()
use lattice_param
use xhf_param
use model_param
implicit none
integer::i,j,k
!set h0xhf
 h0xhf(1:2*Nsite,1:2*Nsite)=hzero(1:2*Nsite,1:2*Nsite)
!set the off diagonal term
do i=1,Nsite,1
   k=0
   do j=1,Dimen,1
     k=k+coor(i,j)
   end do
   h0xhf(Nsite+i,i)=onsitU*delt*((-1)**k)
   h0xhf(i,Nsite+i)=onsitU*delt*((-1)**k)
end do
!call check_Hermite_c(h0xhf,2*Nsite) 
end subroutine set_h0xhf


!Diagonialize h0, give the lowest eigenstate to ph
subroutine input_phiT(nl,nt,h0,ph)
implicit none
integer,intent(IN)::nl
integer,intent(IN)::nt
complex(kind=8),intent(IN)::h0(nl,nl)
complex(kind=8),intent(OUT)::ph(nl,nt)
complex(kind=8)::hu(nl,nl)
real(kind=8)::ev(nl)


call check_Hermite_c(h0,nl)
call zcopy(nl*nl,h0,1,hu,1)
call eigen(Hu,nl,ev)
call zcopy(nl*nt,hu,1,ph,1)

!write(*,*) ev ;call mystop
end subroutine input_phiT


!get the variational energy
subroutine esitmate()
use param
use lattice_param
use xhf_param
use model_param
implicit none
complex(kind=8)::phl(2*Nsite,Ntot)
complex(kind=8)::Amat(2*Nsite,2*Nsite)
integer::sitei,sitej

phl(1:2*Nsite,1:Ntot)=ph(1:2*Nsite,1:Ntot)
call cal_Amat(2*Nsite,Ntot,phl,ph,Amat)

e_var=zero
do sitei=1,2*Nsite,1
   do sitej=1,2*Nsite,1
      e_var=e_var+Hzero(sitei,sitej)*Amat(sitei,sitej)
   enddo
enddo

do sitei=1,Nsite,1
   e_var=e_var+onsitU*Amat(sitei,sitei)*Amat(sitei+Nsite,sitei+Nsite) &
              & -onsitU*Amat(sitei,sitei+Nsite)*Amat(sitei+Nsite,sitei)
end do

end subroutine esitmate
