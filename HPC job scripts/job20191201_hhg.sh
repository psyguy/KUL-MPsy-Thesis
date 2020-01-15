#!/bin/bash -l
#PBS -l nodes=6:ppn=36
#PBS -l walltime=00:59:59
#PBS -l pmem=5gb
#PBS -N December01_hhg
#PBS -m b -M m@gmail.com
 
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1

 
cd $PBS_O_WORKDIR
 
 
module load R/3.6.0-foss-2018a-bare
 
./scripts/20191201_hpc_hhg.R $counter