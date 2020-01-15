#!/bin/bash -l
#PBS -l nodes=2:ppn=24
#PBS -l walltime=1:00:00
#PBS -l pmem=7gb
#PBS -N June29_netviz_2045
#PBS -m b -M m@gmail.com
 
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
 
 
cd $PBS_O_WORKDIR
 
 
module load R/3.5.0-iomkl-2018a-X11-20180131
 
./scripts/20190629_hpc_netviz.R $rounds