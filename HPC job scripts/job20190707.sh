#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=55:00:00
#PBS -l pmem=12gb
#PBS -N July07_1340
#PBS -m b -M m@gmail.com
 
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
 
 
cd $PBS_O_WORKDIR
 
 
module load R/3.5.0-iomkl-2018a-X11-20180131
 
./scripts/20190707_reincarnation.R $index