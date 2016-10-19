#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=00:59:59
#PBS -q normalQ
#PBS -N gd
#PBS -j oe
#PBS -t 1-50

# Change the following to your working directory
module load r320
cd $PBS_O_WORKDIR
Rscript codes.R $PBS_ARRAYID
