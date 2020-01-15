#!/bin/bash -l
#PBS -l nodes=5:ppn=11
#PBS -l walltime=23:59:59
#PBS -l pmem=15gb
#PBS -N July07_continue_1950
#PBS -m b -M mt@gmail.com
 
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
 
 
cd $PBS_O_WORKDIR
 
 
module load R/3.5.0-iomkl-2018a-X11-20180131
 
./scripts/20190707_reincarnation_continue.R $index