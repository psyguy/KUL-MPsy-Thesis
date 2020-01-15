#!/bin/bash -l
#PBS -l nodes=1:ppn=27
#PBS -l walltime=00:20:59
#PBS -l pmem=7gb
#PBS -N July20_2100
#PBS -m b -M m@gmail.com
 
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
 
 
cd $PBS_O_WORKDIR
 
 
module load R/3.5.0-iomkl-2018a-X11-20180131
 
./scripts/20190712_plot_n_glue.R $ind