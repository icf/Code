#!/bin/tcsh
#PBS -N cpqmc
#PBS -l walltime=00:30:00
#PBS -l nodes=1:x5672:ppn=8
cd $PBS_O_WORKDIR
mvp2run ./cpqmc > out
