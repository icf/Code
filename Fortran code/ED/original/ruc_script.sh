#!/bin/bash
export set=3              #set:1 by hand 2 PBC 3 TBC
export Nlattice=9         #Nlattice->set=1
export Nhop=36            #Nhop->set=1
export Dimen=2            #Dimen->set/=1
export Nl1=6              #Nl(1)->set/=1 & Dimen>=1
export Nl2=5              #Nl(2)->set/=1 & Dimen>=2
export Nl3=4              #Nl(3)->set/=1 & Dimen=3
export kbound1=0.85d0     #kbound(1)->set=3 & Dimen>=1
export kbound2=0.632d0    #kbound(2)->set=3 & Dimen>=2
export kbound3=0.38d0     #kbound(3)->set=3 & Dimen>=3
export Nup=1              #Nup
export Ndown=1            #Ndown
export kin=\(-1.d0,0.d0\) #kin  
export onsitU=2.d0        #onsitU
export LanM=20            #LanM
export lit=5.d-13         #lit define Nlan 
export seed=1             #seed
export Nexcit=3           #Nexcit-->excite state
export wstate=1           #wstate-->1 write,0-->not write


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
${kin}  !kin
${onsitU}  !onsitU
${LanM}  !LanM
${lit}  !lit define Nlan
${seed}  !seed
${Nexcit}  !Nexcit-->excite state
${wstate}  !wstate-->1 write,0-->not write
!

export dir=/gpfs/fs1/hshi/nompi_result
#export gz=lanhubbardlx_${Nl1}_Nup_${Nup}_Ndown_${Ndown}_U_${onsitU}_kx_${kbound1}.tar.gz
#export workname=lanhubbardlx_${Nl1}_Nup_${Nup}_Ndown_${Ndown}_U_${onsitU}_kx_${kbound1}
export gz=lanhubbardlx_${Nl1}_ly_${Nl2}_Nup_${Nup}_Ndown_${Ndown}_U_${onsitU}_kx_${kbound1}_ky${kbound2}.tar.gz
export workname=lanhubbardlx_${Nl1}_ly_${Nl2}_Nup_${Nup}_Ndown_${Ndown}_U_${onsitU}_kx_${kbound1}_ky${kbound2}
#export gz=lanhubbardlx_${Nl1}_ly_${Nl2}_lz_${Nl3}_Nup_${Nup}_Ndown_${Ndown}_U_${onsitU}_kx_${kbound1}_ky${kbound2}_kz_${kbound3}.tar.gz
#export workname=lanhubbardlx_${Nl1}_ly_${Nl2}_lz_${Nl3}_Nup_${Nup}_Ndown_${Ndown}_U_${onsitU}_kx_${kbound1}_ky${kbound2}_kz_${kbound3}
tar -cvf ${gz} lanh param hop
rm -rf ${dir}/${workname}
mkdir -p ${dir}/${workname}
mv -f ${gz} ${dir}/${workname}
cd ${dir}/${workname}
tar -xvf ${gz}
rm -rf ${gz}

cat >script <<!
#####################################################################
#!/bin/sh
#$ -S /bin/sh
module load intel/intel-11
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
#$ -cwd
#$ -j y
#$ -q seri2.q # queue name
#$ -N lanhubbard # job name
#$ -l h_cpu=72:00:00 # hard cpu time
./lanh
echo 'done'
#####################################################################
!
qsub script
rm -rf script
