module param
implicit none
!--------------
!The parameters
!--------------
complex(kind=8),parameter:: Xi=dcmplx(0d0,1d0)
complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)
real(kind=8)::Pi

!---------------------------------
!Lattice Hubbard model information
!---------------------------------
integer::Dimen  !the dimension
integer::Nl(3)  !the number of different axis
real(kind=8)::kbound(3) ! the twist boundary condition number
integer::Nlattice    !the number of the whole sites
integer::set    !set the hopt by yourself of by computer
integer::Nup   ! the number of up electrons
integer::Ndn ! the number of down electrons
complex(kind=8)::kin ! the kinetic energy of hubbard model
complex(kind=8)::tprime ! the tprime of hubbard model
real(kind=8)::onsitU ! the onsitU interaction of hubbard model
integer::Nhop ! the number of the hoping terms need to consider


complex(kind=8),allocatable::hopt(:) !hopt and sit record the information
integer,allocatable::sit(:,:)        !of hoping term in different sites
integer,allocatable::coor(:,:) ! we label the site, it record the coordinate
integer(kind=8)::Nhilbert ! the size of the hilbert space
integer(kind=8)::Nhup,Nhdn !the size of the hilbert space of up and dn

!Add when update the new table
integer(kind=8),allocatable::Binomial(:,:)



!-------------------------
!Lanczos method parameters
!-------------------------
Character(len=1)::conv ! E or W converge Energy of wave function
real(kind=8)::varE !The variational energy
logical::waveconv !The converge of wave function
real(kind=8)::lit !the little number define the Nlan
real(kind=8)::lim !lit define orthogonalize,project symmetry
integer::Nexcit,Nei ! Nexcit is the total excited state to calculate, Nei is the step
integer::wstate ! whether to write eigensatet into Hard-disk

!parammeters need to change lanM
integer::LanM ! the lanczos matrix doing each time
integer::LanMmem !remember the LanM to calculate the excited state.
integer::LanMmin !the aparam used to to change LanM
integer::LanMmax
integer::Lanadd
integer::Nlanmatrix_set


complex(kind=8),allocatable::Phi(:,:) ! record the Phi of the lanczos matrix
real(kind=8),allocatable::dig(:),offdig(:),MatrixU(:,:)
integer::Nlan ! define the true dimension of MatrixU
integer::Nlanmatrix !Each time we get one lanczos matrix,Nlanmatrix plus one.
complex(kind=8),allocatable::eigenvalue(:),eigenstate(:,:)

complex(kind=8)::S_square !Calcualte the S^2 of the lattice model
Character(len=1)::PS
integer::two_SS !total s times two
integer::two_rSS !real s times two

complex(kind=8),allocatable::k_int(:) !calculate the k point 2*Pi/Nx*k(1),2*Pi/Ny*k(2),2*Pi/Ny*k(3)
integer::k_keep(3) !which will be kept in the code
Character(len=1)::PK
integer,allocatable::Tmatrix(:,:) !Use to store the nearest hopping in different direction.


integer(kind=8),allocatable::Hm_up(:,:,:),Hm_dn(:,:,:),Tm_up(:,:),Tm_dn(:,:)  !Store the matrix of up and dn H T 
integer,allocatable::Hcoe_up(:,:,:),Hcoe_dn(:,:,:),Tcoe_up(:,:),Tcoe_dn(:,:)  !Store the coe of the matrix 


character(len=200)::basename,eigvname,eigsname ! Writing parameters
character(len=200)::basenameU
character(len=200)::corrSzname,corrSdname,corrDdname
character(len=200)::corrSkname
character(len=200)::coorcicjname,coornkname

integer,parameter::chunk=100 !For the openmp code

!for pinning field
integer:: pinningtype
complex(kind=8):: Hpinn

end module param
