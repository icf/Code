!This page contains the basic subroutine at the beginning of code and the end of
!the code

!1.readparam:Read parameter from 'param' at the beginning of running the code
!2.allocatearray:allocate the arrays at the beginning of the code
!3.deallocatearray:deallocate all the arrays which is allocated.
!4.mystop:deallocate all the arrays and print the running time before stop


!----------------------------------------------------
!This subroutine is used to read parameter from param
!at the beginning of running the code
!---------------------------------------------------- 
subroutine readparam
use param
implicit none
!---------------------------------------
!set the param used in the hubbard model
!---------------------------------------
 open(unit=10,file='param',status='old')
 !set:integer input 
 !    1: we should input the hoping matrix by hand
 !    not 1: the computer set the hoping matrix
 read(10,*) set


 !Nlattice:the whole number of lattice,
 !         when set=1 the Nlattice is useful
 read(10,*) Nlattice


 !Nhop:the total number of the connect lattice,such 
 !     as T12,T14,T21,T23,T32,T34,T41,T43 in 1D,four
 !     lattice model
 !    when set=1 the Nhop is useful
 read(10,*) Nhop


 !Dimen:the Dimension of the lattice
 !    when set/=1 the Dimen is useful
 read(10,*) Dimen


 !Nl(3):the Dimension of the x,y,z axis lattice
 !    when Dim<3,only Nl(1:Dimen) is useful
 !    when set/=1 the Nl(3) is useful
 read(10,*) Nl(1)
 read(10,*) Nl(2)
 read(10,*) Nl(3)


 !kbound(3):the twist boundary condition which is 
 !          twist by exp(I*kbound*Pi) when Nl->1
 !          exp(-I*kbound*Pi) when 1->Nl
 !          when Dim<3,only kbound(1:Dimen) is useful
 !          when set=3, the TBC is used
 read(10,*) kbound(1)
 read(10,*) kbound(2)
 read(10,*) kbound(3)

 read(10,*) PS

 read(10,*) PK

 read(10,*) two_rSS
 read(10,*) k_keep(1)
 read(10,*) k_keep(2)
 read(10,*) k_keep(3)
 !Nup,Ndown:the up and down particle in the lattice
 !          it is always useful whatever the set is.
 read(10,*) Nup
 read(10,*) Ndn

 !kin,onsitU:the param appears in Hubbard model
 !          it is always useful whatever the set is.
 read(10,*) kin
 read(10,*) tprime
 read(10,*) onsitU

 !LanM:The matrix in the lanczos matrix m*m,lit difine Nlan
 read(10,*) LanM
 read(10,*) LanMmin
 read(10,*) LanMmax
 read(10,*) Lanadd
 read(10,*) Nlanmatrix_set
 read(10,*) conv
 read(10,*) lit
 read(10,*) lim

 !Nexcit: the excited state need to be calculate
 !wstate: 1-->write down eigenstate, 0--->do not write
 read(10,*) Nexcit
 read(10,*) wstate

 close(10)

end subroutine readparam




!---------------------------------------------------------------------
!This subroutine allocate all the arrays at the beginning of the code.
!---------------------------------------------------------------------
subroutine allocatearray()
use param
implicit none
allocate(hopt(Nhop),sit(Nhop,2))
allocate(coor(Nlattice,Dimen))
allocate(k_int(Dimen))
allocate(Tmatrix(Nlattice,Dimen))
allocate(Phi(Nhilbert,3))
allocate(dig(LanM))
allocate(offdig(LanM+1))
allocate(MatrixU(LanM,LanM))
allocate(eigenvalue(Nexcit))
allocate(eigenstate(Nhilbert,Nexcit))
allocate(Hm_up(Nhup,Nlattice,Nlattice))
allocate(Hm_dn(Nhdn,Nlattice,Nlattice))
allocate(Tm_up(Nhup,3))
allocate(Tm_dn(Nhdn,3))
allocate(Hcoe_up(Nhup,Nlattice,Nlattice))
allocate(Hcoe_dn(Nhdn,Nlattice,Nlattice))
allocate(Tcoe_up(Nhup,3))
allocate(Tcoe_dn(Nhdn,3))
end subroutine allocatearray

!----------------------------------------------------------------
!deallocatearray subroutine deallocate all the arrays in the code
!----------------------------------------------------------------
subroutine deallocatearray()
use param
implicit none
 if(allocated(hopt)) deallocate(hopt)
 if(allocated(sit)) deallocate(sit)
 if(allocated(coor)) deallocate(coor)
 if(allocated(k_int)) deallocate(k_int)
 if(allocated(Tmatrix)) deallocate(Tmatrix)
 if(allocated(dig)) deallocate(dig)
 if(allocated(offdig)) deallocate(offdig)
 if(allocated(MatrixU)) deallocate(MatrixU) !!
 if(allocated(eigenvalue)) deallocate(eigenvalue) !!
 if(allocated(eigenstate)) deallocate(eigenstate) !!
 if(allocated(Phi)) deallocate(Phi)     !!
 if(allocated(Binomial)) deallocate(Binomial)
 if(allocated(Hm_up)) deallocate(Hm_up)
 if(allocated(Hm_dn)) deallocate(Hm_dn)
 if(allocated(Tm_up)) deallocate(Tm_up)
 if(allocated(Tm_dn)) deallocate(Tm_dn)
 if(allocated(Hcoe_up)) deallocate(Hcoe_up)
 if(allocated(Hcoe_dn)) deallocate(Hcoe_dn)
 if(allocated(Tcoe_up)) deallocate(Tcoe_up)
 if(allocated(Tcoe_dn)) deallocate(Tcoe_dn)
end subroutine deallocatearray

!----------------------------------------------------------------
!This subroutine end the program before deallocate all the arrays
!----------------------------------------------------------------
subroutine mystop
use param
use timing_module
use rand_num
implicit none
call end_genrand()
call deallocatearray()
call EndTiming()
stop
end subroutine mystop

