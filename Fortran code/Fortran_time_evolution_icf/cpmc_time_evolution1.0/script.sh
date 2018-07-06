#!/bin/bash
export set=2           #set:1 by hand 2 TBC
export Nsite=9         #Nsite->set=1
export Nhop=36         #Nhop->set=1
export Dimen=1         #Dimen->set/=1
export Nl1=8           #Nl(1)->set/=1 & Dimen>=1
export Nl2=6           #Nl(2)->set/=1 & Dimen>=2
export Nl3=2           #Nl(3)->set/=1 & Dimen=3
export kbound1=0.8d0    #kbound(1)->set=3 & Dimen>=1
export kbound2=0.d0    #kbound(2)->set=3 & Dimen>=2
export kbound3=0.d0    #kbound(3)->set=3 & Dimen>=3
export t1=-1.d0        #t1 Hubbard hopping t1 in nearest direction.
export onsitU=4.d0     #onsitU Hubbard U interaction on the same site.
export Ntot=8         #Ntot the up spin and down spin particle numbers
export dtype='d'       #dtype for the determinate type: d decouple, c couple.
export Nspin1=4        #Nspin(1) means Nup
export Nspin2=4        #Nspin(2) means Ndn
export crn=-1.d0        #crn the number of release node
export ccoe=0.d0      #ccoe control the release in pop :0.d0 total free wi, 1.d0 wi<phiT|phi>
export max_crn=-1     #max_crn control the release node
export dcp='S'         #dcp The decouple method:'S' or 'C'
export kcrn=2          #kcrn 1.d s;2.d c;3.c s;4 c c:for different kinds of release method
export bgset=0         #bgset 0. mean field 1. dynamic background walker
export fw_bk='FW'      #fw_bk FW means no back propogation,BK means back propogation.
export Nbk1=0         #Nbk(1) is the free projection path
export Nbk2=0        #Nbk(2) is the constraint path
export Nsamples=1     #Nsamples number of sampling do we need in MC process.
export Nwalkers=1000   #Nwalkers number of  walkers
export blockstep=1    #blockstep Block number which inidcate the basic size
export Thermblock=0   #Thermblock number of blocks do we need in thermal process.
export Neqblock=500      #Neqblock  number of blocks do we need in equilibrium step(update and measure)
export PopContrlstep=10  #PopContrlstep  we need to do population control and adjust ET
export StepforGram=10  #StepforGram number of steps when we do modified GS
export meastep=5       #meastep number of steps when we need to do measure
export dt=0.01d0       #dt each slice of imagine time
export PT=1            #PT:0 read from file, 1 set by FS: For the phiT.
export PP=1            #PP:0 read from file, 1 set by phiT: For the phi.
export max_ad=100      #max_ad for free-projection
export Nmpi=1

cat >param <<!
${set}        !set:1 by hand 2 TBC
${Nsite}      !Nsite->set=1
${Nhop}       !Nhop->set=1
${Dimen}      !Dimen->set/=1
${Nl1}        !Nl(1)->set/=1 & Dimen>=1
${Nl2}        !Nl(2)->set/=1 & Dimen>=2
${Nl3}        !Nl(3)->set/=1 & Dimen=3
${kbound1}    !kbound(1)->set=3 & Dimen>=1
${kbound2}    !kbound(2)->set=3 & Dimen>=2
${kbound3}    !kbound(3)->set=3 & Dimen>=3
(${t1},0.d0)     !t1 Hubbard hopping t1 in nearest direction.
${onsitU}     !onsitU Hubbard U interaction on the same site.
${Ntot}       !Ntot the up spin and down spin particle numbers
${dtype}      !dtype for the determinate type: d decouple, c couple.
${Nspin1}     !Nspin(1) means Nup
${Nspin2}     !Nspin(2) means Ndn
${crn}        !crn the number of release node
${ccoe}       !ccoe control the release in pop :1.d0 total free wi, 0.d0 wi<phiT|phi>
${max_crn}    !max_crn control the release node
${dcp}        !dcp The decouple method:'S' or 'C'
${kcrn}       !kcrn 1.d s;2.d c;3.c s;4 c c:for different kinds of release method
${bgset}      !bgset 0. mean field 1. dynamic background walker
${fw_bk}      !fw_bk FW means no back propogation,BK means back propogation.
${Nbk1}       !Nbk(1) is the free projection path
${Nbk2}       !Nbk(2) is the constraint path
${Nsamples}   !Nsamples number of sampling do we need in MC process.
${Nwalkers}   !Nwalkers number of  walkers
${blockstep}  !blockstep Block number which inidcate the basic size
${Thermblock} !Thermblock number of blocks do we need in thermal process.
${Neqblock}   !Neqblock  number of blocks do we need in equilibrium step(update and measure)
${PopContrlstep} !PopContrlstep  we need to do population control and adjust ET
${StepforGram}   !StepforGram number of steps when we do modified GS
${meastep}    !meastep number of steps when we need to do measure
${dt}         !dt each slice of imagine time
${PT}         !PT:0 read from file, 1 set by FS: For the phiT.
${PP}         !PP:0 read from file, 1 set by phiT: For the phi.
${max_ad}     !max_ad for free-projection
!

export gz=m${Nmpi}-gdqmc1.3-${set}-${Nsite}-${Nhop}-${Dimen}-${Nl1}-${Nl2}-${Nl3}-${kbound1}-${kbound2}-${kbound3}-${t1}-${onsitU}-${Ntot}-${dtype}-${Nspin1}-${Nspin2}-${crn}-${ccoe}-${max_crn}-${dcp}-${kcrn}-${bgset}-${fw_bk}-${Nbk1}-${Nbk2}-${Nsamples}-${Nwalkers}-${blockstep}-${Thermblock}-${Neqblock}-${PopContrlstep}-${StepforGram}-${meastep}-${dt}-${PT}-${PP}-${max_ad}.tar.gz
export wkdir=m${Nmpi}-gdqmc1.3-${set}-${Nsite}-${Nhop}-${Dimen}-${Nl1}-${Nl2}-${Nl3}-${kbound1}-${kbound2}-${kbound3}-${t1}-${onsitU}-${Ntot}-${dtype}-${Nspin1}-${Nspin2}-${crn}-${ccoe}-${max_crn}-${dcp}-${kcrn}-${bgset}-${fw_bk}-${Nbk1}-${Nbk2}-${Nsamples}-${Nwalkers}-${blockstep}-${Thermblock}-${Neqblock}-${PopContrlstep}-${StepforGram}-${meastep}-${dt}-${PT}-${PP}-${max_ad}


tar -cvf ${gz} cpqmc param *.dat
mv ${gz} ../gdqmc-result
cd ../gdqmc-result
rm -rf ${wkdir}
mkdir ${wkdir}
mv ${gz} ${wkdir}
cd ${wkdir}
tar -xvf ${gz}
rm -rf ${gz}

nohup mpirun -np ${Nmpi} cpqmc &

