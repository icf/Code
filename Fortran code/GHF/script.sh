#!/bin/bash
OLDIFS=$IFS
IFS=$'\n'
export set=2           #set:1 by hand 2 TBC
export Nsite=9         #Nsite->set=1
export Nhop=36         #Nhop->set=1
export Dimen=2         #Dimen->set/=1
export Nl1=20          #Nl(1)->set/=1 & Dimen>=1
export Nl2=20           #Nl(2)->set/=1 & Dimen>=2
export Nl3=2           #Nl(3)->set/=1 & Dimen=3
export kbound1=0.0d0    #kbound(1)->set=3 & Dimen>=1
export kbound2=0.0d0    #kbound(2)->set=3 & Dimen>=2
export kbound3=0.d0    #kbound(3)->set=3 & Dimen>=3
export t1=-1.d0        #t1 Hubbard hopping t1 in nearest direction.
export onsitU=2.d0    #onsitU Hubbard U interaction on the same site.
export Ntot=400          #Ntot the up spin and down spin particle numbers
export dtype='c'       #dtype for the determinate type: d decouple, c couple.
export Nspin1=200        #Nspin(1) means Nup
export Nspin2=200       #Nspin(2) means Ndn
export dmin=0.00d0     #order parameter min
export dmax=0.90d0     #order parameter max 
export dstep=0.01d0    #order parameter add
export Nmpi=1

dir_now=`pwd`
mkdir /storm/data10/mingpu/sxhf_"$Nl1"_"$Nl2"_U"$onsitU"
cp  for_k_point/k_point_1 ./
sort -k2n k_point_1 | uniq > k_point_2
gawk '{$1=$1"d0";$2=$2"d0";print $0}' k_point_2 > k_point
file="k_point"
for kxky in `cat $file`
do
 # echo $kxky| gawk '{kx=$1;ky=$2;print kx;print ky}'
  kbound1=`echo $kxky| gawk '{print $1}'`
  kbound2=`echo $kxky| gawk '{print $2}'`
# this is from the original script
cat >param <<!
${set}        !set:1 by hand 2 TBC
${Nsite}      !Nsite->set=1
${Nhop}       !Nhop->set=1
${Dimen=2}    #Dimen->set/=1
${Nl1}        !Nl(1)->set/=1 & Dimen>=1
${Nl2}        !Nl(2)->set/=1 & Dimen>=2
${Nl3}        !Nl(3)->set/=1 & Dimen=3
${kbound1}    !kbound(1)->set=3 & Dimen>=1
${kbound2}    !kbound(2)->set=3 & Dimen>=2
${kbound3}    !kbound(3)->set=3 & Dimen>=3
(${t1},0.d0)  !t1 Hubbard hopping t1 in nearest direction.
${onsitU}     !onsitU Hubbard U interaction on the same site.
${Ntot}       !Ntot the up spin and down spin particle numbers
${dtype}      !dtype for the determinate type: d decouple, c couple.
${Nspin1}     !Nspin(1) means Nup
${Nspin2}     !Nspin(2) means Ndn
${dmin}       !order parameter min
${dmax}       !order parameter max 
${dstep}      !order parameter add
!

export gz=sxhf-m${Nmpi}-${set}-${Nsite}-${Nhop}-${Dimen}-${Nl1}-${Nl2}-${Nl3}-${kbound1}-${kbound2}-${kbound3}-${t1}-${onsitU}-${Ntot}-${dtype}-${Nspin1}-${Nspin2}-${dmim}-${dmax}-${dstep}.tar.gz
export wkdir=sxhf-m-${Nl1}-${Nl2}-${Nl3}-${kbound1}-${kbound2}-${kbound3}-${t1}-${onsitU}-${Ntot}-${dtype}-${Nspin1}-${Nspin2}

tar -cf ${gz} sxhf param
#tar -cf ${gz} afqmc param
#mv ${gz} ../result_"$Nl1"_"$Nl2"/
mv ${gz} /storm/data10/mingpu/sxhf_"$Nl1"_"$Nl2"_U"$onsitU"
#cd ../result_"$Nl1"_"$Nl2"/
cd /storm/data10/mingpu/sxhf_"$Nl1"_"$Nl2"_U"$onsitU"
rm -rf ${wkdir}
mkdir ${wkdir}
mv ${gz} ${wkdir}
cd ${wkdir}
tar -xf ${gz}
rm -rf ${gz}


cat >script <<!
#!/bin/tcsh
#PBS -N sxhf
#PBS -l walltime=1:00:00
#PBS -l nodes=${Nmpi}:rain:ppn=1
cd /storm/data10/mingpu/sxhf_"$Nl1"_"$Nl2"_U"$onsitU"/${wkdir}
./sxhf  > OUTFILE

!

qsub script
#./sxhf
cd $dir_now
done
IFS=$OLDIFS
echo "qsub---done"

#nohup mpirun -np ${Nmpi} cpqmc &
#SBATCH --exclusive
#SBATCH -C new
