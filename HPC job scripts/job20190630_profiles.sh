#!/bin/bash -l
#PBS -l nodes=12:ppn=9
#PBS -l walltime=0:30:00
#PBS -l pmem=20gb
#PBS -N June30_profiles_onlycoefs_2145
#PBS -m b -M m@gmail.com
 
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
 
 
cd $PBS_O_WORKDIR
 
 
module load R/3.5.0-iomkl-2018a-X11-20180131
 
./scripts/20190630_hpc_profile.R $index