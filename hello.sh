#!/bin/sh --login
#PBS -N simu
#PBS -t 1-102
#PBS -l nodes=1:ppn=1,walltime=19:59:59
#PBS -j oe

cd $PBS_O_WORKDIR
./run_mainPro.sh /share/apps/MATLAB/R2015a $PBS_ARRAYID
