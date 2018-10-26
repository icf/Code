!---------------------------------------------
!This parameter contains the useful parameters
!---------------------------------------------
module param
implicit none
 complex(kind=8),parameter:: Xi=dcmplx(0d0,1d0)
 complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
 complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)
 !real(kind=8)::Pi=dacos(-1d0)
 real(kind=8)::Pi=3.1415926535898d0
end module param



!----------------------------------------------------------
!This parameter contains all the one with lattice and Hzero
!----------------------------------------------------------
module lattice_param
implicit none
 integer::set !set the hopt by yourself or by computer
 integer::Nsite !the number of the whole sites
 integer::Nhop ! the number of the hoping terms need to consider
 integer::Dimen !the dimension
 integer::Nl(3) !the number of different axis
 real(kind=8)::kbound(3) !the twist boundary condition number
 complex(kind=8),allocatable::hopt(:) !hopt and sit record the information
 integer,allocatable::sit(:,:)        !of hoping term in different sites
 integer,allocatable::coor(:,:) ! we label the site, it record the coordinate
 integer,allocatable::Tmatrix(:,:) !Use to store the nearest hopping in different direction.
 complex(kind=8),allocatable::Hzero(:,:)    !2*Nsite,2*Nsite,the Hezo Hamiltonian of the lattice
end module lattice_param



!-------------------------------------------
!This parameter contains all the model_param
!-------------------------------------------
module model_param
implicit none
 complex(kind=8):: t1      !Hubbard hopping t1 in nearest direction.
 real(kind=8):: onsitU     !Hubbard U interaction on the same site.
 integer:: Ntot            !the tot number of electrons
 character(len=1):: dtype  !for the determinate type: d decouple, c couple.
 integer:: Nspin(2)        !Nup and Ndn
end module model_param


module xhf_param
 complex(kind=8)::e_var  !variational energy of the xhf wf
 complex(kind=8),allocatable:: h0xhf(:,:) !The hartree fork hamiltonian
 complex(kind=8),allocatable:: ph(:,:)  ! The wf of the hf hamiltonian
 complex(kind=8),allocatable:: ph_save(:,:)  ! The wf of the hf hamiltonian
 real(kind=8)::E_save=1d20
 real(kind=8)::delt,dmin,dmax,dstep
 real(kind=8)::lit=1d-8
end module xhf_param


!--------------------------
!mpi or serial  parammeters
!--------------------------
module mpi_serial_param
implicit none
 integer::ierr
 integer::rank
 integer::Nsize
 integer::htype
 integer::phtype
end module mpi_serial_param
