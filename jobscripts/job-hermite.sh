#!/bin/bash
#PBS -N run_128
#PBS -l mppwidth=64
#PBS -l mppnppn=32
#PBS -l walltime=01:00:00             
#PBS -M volker.springel@h-its.org
#PBS -m abe 

  
# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

# Launch the parallel job to the allocated compute nodes
aprun -n 64 -N 32 ./Arepo  param.txt  


