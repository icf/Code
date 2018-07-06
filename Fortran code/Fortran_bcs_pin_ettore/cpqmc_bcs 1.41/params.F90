!---------------------------------------------
!This parameter contains the useful parameters
!---------------------------------------------
module param
implicit none
 complex(kind=8),parameter:: Xi=dcmplx(0d0,1d0)
 complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
 complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)
 complex(kind=8),parameter:: Nepero=dcmplx(2.71828182845905d0,0.d0)
 !real(kind=8)::Pi=dacos(-1d0)
 real(kind=8)::Pi=3.1415926535898d0
end module param



!----------------------------------------------------------
!This parameter contains all the one with lattice and Hzero
!----------------------------------------------------------
module lattice_param
implicit none
 integer::set !set the hopt by yourself or by computer
 integer::Nbravais !the number of sites of the Bravais lattice
 integer::I_ReadH0  !read the non interacting H from file
 integer::Nbands !the number of bands
 integer::Nsite !the number of the whole sites
 integer::Nhop ! the number of the hoping terms need to consider
 integer::Dimen !the dimension
 integer::Nl(3) !the number of different axis
 real(kind=8)::kbound(3) !the twist boundary condition number
 character(len=1):: boundary ! = o/c   open (y direction) / closed
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
 complex(kind=8):: tprime         !t prime parameter (diagonal hop term)
 real(kind=8):: onsitU     !Hubbard U interaction on the same site.
 integer:: Ntot            !the tot number of electrons
 character(len=1):: dtype  !for the determinate type: d decouple, c couple.
 integer:: Nspin(2)        !Nup and Ndn
 integer:: Nzeta            !Nup-Ndn
 complex(kind=8):: tpp,tpd  !hopping for 3bands HM
 real(kind=8):: epsp,epsd   !energies for 3bands HM
 real(kind=8):: onsitUp,onsitUd !same site interaction for 3bands HM
 real(kind=8):: Vpd   !pd interaction for 3bands HM
 integer::ipinn  !if zero no pinning
 real(kind=8):: Hpinn  !pinning field
 integer:: pinningtype !type of pinning 
 character(len=1):: pinndir !pinning field direction 
 integer::I_Vext  !if zero no additional external field, if one read Vext 
end module model_param



!-------------------------------------------
!This parameter contains all the model_param
!-------------------------------------------
module method_param
implicit none
 real(kind=8)::crn  !Decide the kind of update and mesurement.
 real(kind=8)::ccoe !control the release in pop :0.d0 total free wi, 1.d0 wi<phiT|phi>
 integer::max_crn   !Decide the kind of methods.
 character(len=1)::dcp    !The decouple method
 integer::kcrn      !Decide the different kind of free projection
 integer::bgset     !0. mean field 1. dynamic background walker 
 character(len=2)::fw_bk   !Which kind of measure method: FW,BK.
 integer::Nstps_fwd       ! After each fwd step we do the back propagate
                          !measurement
 integer::Nbk(2)    !Nbk(1) is the free projection Nbk(2) is the cpqmc
 logical::back_pro !If logical back_pro is .true. we do it.
 integer::i_back   !Use to define the step to meas in back_propagate
 integer::I_obdm   !Compute (1) or not compute (0) averager obdm
 integer::Nbeta    !Imaginary time window for dynamical correlations
 integer::Npair    !Number of pairing correlations
 integer::I_twob   !Compute (1) or not compute (0)  two body dynamical
 integer::I_onebf  !Compute (1) or not compute (0) full Green function old method

 integer::GM_input_flag !GM input: GM_input_flag=1; 
 !decomposition method for BCS
 integer::decMethod !1:DET; 2:Analytic

end module method_param



!----------------------------------------------
!This parameter contains all the QMC loop param
!----------------------------------------------
module mc_loop_param
implicit none
 integer::Nsamples        ! number of sampling do we need in MC process.
 integer::Nwalkers        ! number of  walkers 
 integer::MNwalkers       ! number of MPI walkers
 integer::blockstep       ! Block number which inidcate the basic size in measure
                          !and thermals.
 integer::Thermblock      ! number of blocks do we need in thermal process.
 integer::Neqblock        ! number of blocks do we need in equilibrium step (Where we
                          ! update and measure )
 integer::max_local       !Neqblock*blockstep Use to define the max step when we get the local quantity.
                          !Similar to max_therm,it is the max local meas step
 integer::PopContrlstep   ! number of steps we need to do population control.
                          !We also use each population control step to adjust ET
 integer::StepforGram     ! number of steps when we do modified GS
 integer::StepforPure     ! number of steps when we stabiliza pure estimators
 integer::meastep         ! number of steps when we need to do measure.(The toal
                          !step is meastep plus Nstps_fwd) 
end module mc_loop_param



!-------------------------------------------------
!This parameter contains all the QMC project param
!-------------------------------------------------
module project_param
implicit none
 real(kind=8):: dt        ! each slice of imagine time
 complex(kind=8),allocatable::exp_halfK(:,:)   ! 2*Nsite,2*Nsite.
 complex(kind=8),allocatable::exp_mhalfK(:,:)   ! 2*Nsite,2*Nsite.
 complex(kind=8),allocatable::exp_K(:,:)   ! 2*Nsite,2*Nsite.
 complex(kind=8)::gama      !the param used in H-S transforamtion.
 complex(kind=8)::gamad     !gamma for d orbitals
 complex(kind=8)::gamax     !gamma for px orbitals 
 complex(kind=8)::gamay     !gamma for py orbitals
 complex(kind=8)::gamaf     !the gama for the free-projection
 complex(kind=8)::expln_up(2),expln_dn(2) !for the cpmc update
 complex(kind=8)::explnd_up(2),explnd_dn(2) !for the cpmc update
 complex(kind=8)::explnx_up(2),explnx_dn(2) !for the cpmc update
 complex(kind=8)::explny_up(2),explny_dn(2) !for the cpmc update
 real(kind=8),allocatable::ng(:) !The back ground for free projection
end module project_param



!----------------------------------------------
!This parameter contains all the QMC phiT param
!----------------------------------------------
module phiT_param
implicit none
 integer::I_wavefun   ! 1 for Slater Determinants, 2 for BCS wave function
 integer::PT
 integer::Dtot   !Number of multi determinate try wave function
 complex(kind=8),allocatable:: phiT(:,:,:)  !Nsite,Ntot,Dtot.try wavefunction and initia wavefunction_up
 complex(kind=8),allocatable:: coe_multi(:)    !Dtot.the coefficient of different walkers
 complex(kind=8),allocatable:: FPairing(:,:)  !Nsite,Nsite
 complex(kind=8),allocatable:: DUnpaired(:,:)  !Nsite, Nspin(1)-Nspin(2)
end module phiT_param



!-----------------------------------------------------------------------
!In CPMC: |psi>= sum_i rx(i) weight(i) phi(:,:,i)/tot_imp(i)
!In FPMC and RCPMC: |psi>= sum_i rx(i) weight(i) phi(:,:,i)
!-----------------------------------------------------------------------
module phi_param
implicit none
 integer::PP               !Set phi 1 from phiT or 0 read from file.
 complex(kind=8),allocatable::phi(:,:,:)       !2*Nsite,Ntot,Nwalkers determinate.
 complex(kind=8),allocatable::phi0(:,:,:)      !2*Nsite,Ntot,Nwalkers wave function in BK.
 !complex(kind=8),allocatable::invop(:,:,:)       !Ntot,Ntot,Nwalkers the inverse of overlap
 complex(kind=8),allocatable::impfunc(:,:)     !Dtot,Nwalkers  important function.
 complex(kind=8),allocatable::tot_imp(:)       !Nwalkers The total important function.
 complex(kind=8)::abs_imp_avg,imp_avg,imp_err  !The average of total important and the error
 real(kind=8),allocatable::weight(:)           !Nwalkers The weight of each walkers.
 real(kind=8),allocatable::Mweight(:)          !Mwalkers the total weight of diffrent thread.
 real(kind=8),allocatable::dlogw(:)            !Nwalkers The dlog(weight) of each walkers.
 complex(kind=8),allocatable::rx(:)            !Nwalkers phase of wave function.
 real(kind=8),allocatable::back_store(:,:,:) !Nsite,Nstps_fwd,Nwalkers back projection field.
end module phi_param


!--------------------------------------------------------------
!This module contains the file need to be backup in rcpmc or bk
!--------------------------------------------------------------
module backup_param
implicit none
integer::i_pop_tmp,i_GS_tmp
real(kind=8)::crn_tmp
complex(kind=8),allocatable::phi_tmp(:,:,:)
real(kind=8),allocatable::weight_tmp(:)
real(kind=8),allocatable::dlogw_tmp(:)
complex(kind=8),allocatable::rx_tmp(:)
complex(kind=8),allocatable::impfunc_tmp(:,:)
complex(kind=8),allocatable::tot_imp_tmp(:)
end module backup_param


!----------------------
!Param used in adjustET
!----------------------
module adET_param
implicit none
 integer::i_ad       !number of adjust ET
 real(kind=8)::ET,etrial    !The try energy which is first used in the code.
 integer::max_ad    !The maximum adjust ET in free projection.
 real(kind=8)::E_sum,E_sum_old             !The sum in the adjust ET
 real(kind=8)::m_w,m_w_old                 !The mean weight during one pop control,used
                                          !in adjust ET
end module adET_param



!-----------------------------------
!Parameter contain the measure param
!-----------------------------------
module meas_param
implicit none
 complex(kind=8),allocatable::ZetaN_l(:,:,:)
 real(kind=8),allocatable::kin_l(:,:),v_l(:,:),e_l(:,:),eE_l(:,:)
 complex(kind=8),allocatable::var_l(:,:)
 real(kind=8),allocatable::nu_l(:,:),nd_l(:,:)
 real(kind=8),allocatable::sisj_l(:,:,:)
 real(kind=8),allocatable::ninj_l(:,:,:,:,:)
 real(kind=8),allocatable::szsz_l(:,:,:,:,:)
 real(kind=8),allocatable::sk_l(:,:,:)
 complex(kind=8),allocatable::obdm_l(:,:,:,:)
 complex(kind=8),allocatable::didj_l(:,:,:,:)
 complex(kind=8),allocatable::cacb_l(:,:,:,:,:)
 complex(kind=8),allocatable::nofr_l(:,:,:,:,:)
 complex(kind=8),allocatable::cicj_l(:,:,:)
 complex(kind=8),allocatable::cicj_t_l(:,:,:,:)
 complex(kind=8),allocatable::cicjh_t_l(:,:,:,:)
 complex(kind=8),allocatable::GreenP_t_l(:,:,:)
 complex(kind=8),allocatable::GreenH_t_l(:,:,:)
! complex(kind=8),allocatable::rho_t_l(:,:,:,:)
 complex(kind=8),allocatable::nupnup_t_l(:,:,:,:)
 complex(kind=8),allocatable::ndnnup_t_l(:,:,:,:)
 real(kind=8),allocatable::ck_l(:,:,:)
 complex(kind=8),allocatable::sig(:)!Use to decide the sign of the free project problem
 complex(kind=8),allocatable::absig(:)
!--------------------------------
 complex(kind=8),allocatable::cicj_l_global(:,:,:,:)
 character(len=300)::basename,EnergyName,numName,scorrName,ncorrName,skName,ckName,cijtName
 character(len=300)::DcorrName,cabName,nofrName,obdmName,zetaNName,szcorrName
 character(len=300)::cijhtName,nupnuptName,ndnnuptName,GreenPtName,GreenHtName
!--------------------------------
 character(len=300)::cdensityName,sdensityName

end module meas_param


!--------------------------------------------
!The one meas parameter in measure subroutine
!--------------------------------------------
module one_meas_param
implicit none
complex(kind=8)::K_one
complex(kind=8)::KE_one
complex(kind=8)::V_one
complex(kind=8)::VE_one
complex(kind=8)::E_one
complex(kind=8)::EE_one
complex(kind=8)::var_one
complex(kind=8)::nu_one
complex(kind=8)::nd_one
complex(kind=8),allocatable::ZetaN_one(:)
complex(kind=8),allocatable::Cab_one(:,:,:)
complex(kind=8),allocatable::Obdm_one(:,:)
complex(kind=8),allocatable::Nr_one(:,:,:)
complex(kind=8),allocatable::S_one(:)
complex(kind=8),allocatable::N_one(:,:,:)
complex(kind=8),allocatable::Sz_one(:,:,:)
complex(kind=8),allocatable::D_one(:,:)
complex(kind=8),allocatable::c_one(:)
complex(kind=8)::Sig_one
complex(kind=8)::AbSig_one
complex(kind=8),allocatable::GP_one(:,:)
complex(kind=8),allocatable::GH_one(:,:)
!complex(kind=8),allocatable::RHO_one(:,:)
complex(kind=8),allocatable::nupnup_one(:,:)
complex(kind=8),allocatable::ndnnup_one(:,:)
complex(kind=8),allocatable::mu0(:),mu(:)
complex(kind=8),allocatable::GreenP_one(:)
complex(kind=8),allocatable::GreenH_one(:)
!--------------------------------------
complex(kind=8),allocatable::c_one_global(:,:)
end module one_meas_param


!--------------------------
!S.C. loop para    
!--------------------------
module sc_loop_param
implicit none

complex(kind=8),allocatable::cicj_sc_global(:,:),cicj_sc_global1(:,:)
complex(kind=8),allocatable::phi_sc(:,:,:)       !2*Nsite,Ntot,Dtot determinate.
complex(kind=8),allocatable::BCS_sc(:,:)         !Nsite,Nsite
integer::sc_loop_flag !0-->stop get_cc_phiT in get_phiT, 1-->turn on get_cc_phiT in get_phiT
integer::sc_ite_flag  !0-->don't iterate
integer::sc_loop_num ! number of loop steps
integer::sc_loop_Na  ! size of A= matrix(<phiT|phi>)
integer::sc_step_counter ! sc step counter
!---------------------
integer::sc_initial_choose !1 or 2:HF or free electrons; 0:BCS simulation from input cicj_sc_global 

end module sc_loop_param
    
 
!--------------------------
!S.C. loop para    
!--------------------------
module HF_param
implicit none

complex(kind=8)::HF_a
integer::HF_num ! number of HF steps


end module HF_param
    
    
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
 integer::phtype2
end module mpi_serial_param
