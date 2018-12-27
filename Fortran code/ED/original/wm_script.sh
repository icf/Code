#!/bin/bash
export set=3              #set:1 by hand 2 PBC 3 TBC
export Nlattice=9         #Nlattice->set=1
export Nhop=36            #Nhop->set=1
export Dimen=2            #Dimen->set/=1
export Nl1=3              #Nl(1)->set/=1 & Dimen>=1
export Nl2=3              #Nl(2)->set/=1 & Dimen>=2
export Nl3=3              #Nl(3)->set/=1 & Dimen=3
export kbound1=0.62d0     #kbound(1)->set=3 & Dimen>=1
export kbound2=0.20d0    #kbound(2)->set=3 & Dimen>=2
export kbound3=0.38d0     #kbound(3)->set=3 & Dimen>=3
export Nup=2              #Nup
export Ndown=2            #Ndown
export kin=-1.d0          #kin  
export onsitU=2.d0        #onsitU
export LanM=20            #LanM
export LanMmin=20         #LanMmin
export LanMmax=100        #LanMmax
export Lanadd=10          #Lanadd
export Nlanmatrix_set=10  #Nlanmatrix_set
export lit=5.d-13         #lit define Nlan 
export seed=1             #seed
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
${Nup}  !Nup
${Ndown}  !Ndown
(${kin},0.d0)  !kin
${onsitU}  !onsitU
${LanM}  !LanM
${LanMmin} !LanMmin
${LanMmax} !LanMmax
${Lanadd}  !Lanadd
${Nlanmatrix_set} !Nlanmatrix_set
${lit}  !lit define Nlan
${seed}  !seed
${Nexcit}  !Nexcit-->excite state
${wstate}  !wstate-->1 write,0-->not write
!

export gz=lanhubbard_nompi3.1_${set}_${Nlattice}_${Nhop}_${Dimen}_${Nl1}_${Nl2}_${Nl3}_${kbound1}_${kbound2}_${kbound3}_${Nup}_${Ndown}_${kin}_${onsitU}_${LanM}_${LanMmin}_${LanMmax}_${Lanadd}_${Nlanmatrix_set}_${lit}_${seed}_${Nexcit}_${wstate}.tar.gz
export wkdir=lanhubbard_nompi3.1_${set}_${Nlattice}_${Nhop}_${Dimen}_${Nl1}_${Nl2}_${Nl3}_${kbound1}_${kbound2}_${kbound3}_${Nup}_${Ndown}_${kin}_${onsitU}_${LanM}_${LanMmin}_${LanMmax}_${Lanadd}_${Nlanmatrix_set}_${lit}_${seed}_${Nexcit}_${wstate}

tar -cvf ${gz} lanh param hop
mv ${gz} ../lan_result/
cd ../lan_result/
rm -rf ${wkdir}
mkdir ${wkdir}
mv ${gz} ${wkdir}
cd ${wkdir}
tar -xvf ${gz}
rm -rf ${gz}


cat >script <<!
#!/bin/csh

#SBATCH -J hubblan
#SBATCH -n 1
#SBATCH -t 1-00:00:00

./lanh > OUTFILE
!
sbatch script

echo "qsub---done"
