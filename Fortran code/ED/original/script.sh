#!/bin/bash
export set=3              #set:1 by hand 2 PBC 3 TBC
export Nlattice=9         #Nlattice->set=1
export Nhop=36            #Nhop->set=1
export Dimen=1            #Dimen->set/=1
export Nl1=14              #Nl(1)->set/=1 & Dimen>=1
export Nl2=4              #Nl(2)->set/=1 & Dimen>=2
export Nl3=2              #Nl(3)->set/=1 & Dimen=3
export kbound1=0.d0     #kbound(1)->set=3 & Dimen>=1
export kbound2=0.d0    #kbound(2)->set=3 & Dimen>=2
export kbound3=0.d0     #kbound(3)->set=3 & Dimen>=3
export PS='N'             #PS for project S,'Y' or 'N'
export PK='N'             #PK for project K,'Y' or 'N'
export two_rSS=0          #two_rSS useually is ABS(Nup-Ndn)
export k_keep1=0          #k_keep(1)->Dimen>=1
export k_keep2=0          #k_keep(2)->Dimen>=2
export k_keep3=0          #k_keep(3)->Dimen=3
export Nup=7              #Nup
export Ndown=7           #Ndown
export kin=-1.d0          #kin  
export onsitU=4.d0        #onsitU
export LanM=10            #LanM
export LanMmin=10         #LanMmin
export LanMmax=50         #LanMmax
export Lanadd=10          #Lanadd
export Nlanmatrix_set=10  #Nlanmatrix_set
export conv='E'           #conv E or W converge of energy or wave function
export lit=5.d-10         #lit define Nlan 
export lim=1.d-2          #lim define orthogonalize,project symmetry
export Nexcit=1           #Nexcit-->excite state
export wstate=0           #wstate-->1 write,0-->not write


cat >param <<!
${set}  !set:1 by hand 2 PBC 3 TBC
${Nlattice}  !Nlattice->set=1
${Nhop}  !Nhop->set=1
${Dimen}  !Dimen->set/=1
${Nl1}  !Nl(1)->set/=1 & Dimen>=1
${Nl2}  !Nl(2)->set/=1 & Dimen>=2
${Nl3}  !Nl(3)->set/=1 & Dimen=3
${kbound1}  !kbound(1)->set=3 & Dimen>=1
${kbound2}  !kbound(2)->set=3 & Dimen>=2
${kbound3}  !kbound(3)->set=3 & Dimen>=3
${PS}         !PS for project S,'Y' or 'N'
${PK}         !PK for project K,'Y' or 'N'
${two_rSS}    !two_rSS useually is ABS(Nup-Ndn)
${k_keep1}    !k_keep(1)->Dimen>=1
${k_keep2}    !k_keep(2)->Dimen>=2
${k_keep3}    !k_keep(3)->Dimen=3
${Nup}  !Nup
${Ndown}  !Ndown
(${kin},0.d0)  !kin
${onsitU}  !onsitU
${LanM}  !LanM
${LanMmin} !LanMmin
${LanMmax} !LanMmax
${Lanadd}  !Lanadd
${Nlanmatrix_set} !Nlanmatrix_set
${conv} !conv E or W converge of energy or wave function
${lit}  !lit define Nlan
${lim}  !lim define orthogonalize,project symmetry
${Nexcit}  !Nexcit-->excite state
${wstate}  !wstate-->1 write,0-->not write
!

export gz=lanhubbard_openmp4.3_${set}_${Nlattice}_${Nhop}_${Dimen}_${Nl1}_${Nl2}_${Nl3}_${kbound1}_${kbound2}_${kbound3}_${PS}_${PK}_${two_rSS}_${k_keep1}_${k_keep2}_${k_keep3}_${Nup}_${Ndown}_${kin}_${onsitU}_${LanM}_${LanMmin}_${LanMmax}_${Lanadd}_${Nlanmatrix_set}_${conv}_${lit}_${lim}_${Nexcit}_${wstate}.tar.gz
export wkdir=lanhubbard_openmp4.3_${set}_${Nlattice}_${Nhop}_${Dimen}_${Nl1}_${Nl2}_${Nl3}_${kbound1}_${kbound2}_${kbound3}_${PS}_${PK}_${two_rSS}_${k_keep1}_${k_keep2}_${k_keep3}_${Nup}_${Ndown}_${kin}_${onsitU}_${LanM}_${LanMmin}_${LanMmax}_${Lanadd}_${Nlanmatrix_set}_${conv}_${lit}_${lim}_${Nexcit}_${wstate}

#export gz=lanhubbard_openmp4.2_U_0.0_15.0_0.5_${set}_${Nlattice}_${Nhop}_${Dimen}_${Nl1}_${Nl2}_${Nl3}_${kbound1}_${kbound2}_${kbound3}_${PS}_${PK}_${two_rSS}_${k_keep1}_${k_keep2}_${k_keep3}_${Nup}_${Ndown}_${kin}_${onsitU}_${LanM}_${LanMmin}_${LanMmax}_${Lanadd}_${Nlanmatrix_set}_${conv}_${lit}_${lim}_${Nexcit}_${wstate}.tar.gz
#export wkdir=lanhubbard_openmp4.2_U_0.0_15.0_0.5_${set}_${Nlattice}_${Nhop}_${Dimen}_${Nl1}_${Nl2}_${Nl3}_${kbound1}_${kbound2}_${kbound3}_${PS}_${PK}_${two_rSS}_${k_keep1}_${k_keep2}_${k_keep3}_${Nup}_${Ndown}_${kin}_${onsitU}_${LanM}_${LanMmin}_${LanMmax}_${Lanadd}_${Nlanmatrix_set}_${conv}_${lit}_${lim}_${Nexcit}_${wstate}


tar -cvf ${gz} lanh param hop
mv ${gz} ../openmp4.3_result/
cd ../openmp4.3_result/
rm -rf ${wkdir}
mkdir ${wkdir}
mv ${gz} ${wkdir}
cd ${wkdir}
tar -xvf ${gz}
rm -rf ${gz}

nohup ./lanh &
echo "qsub---done"
