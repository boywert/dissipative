#!/bin/bash 
#SBATCH -n 128                    #Number of cores 
#SBATCH -p vogelsberger           #Partition to submit to 
#SBATCH --mem-per-cpu=1500        #Memory per cpu in MB (see also --mem)
#SBATCH --ntasks-per-node=64      #asks for 64 cores per node as in the script hostgen.
#SBATCH --time=48:00:00 
#SBATCH -J arepotest
#SBATCH -o ../log_intel/arepotest.%j.out
#SBATCH -e ../log_intel/arepotest.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fmarinac@mit.edu

if [ -f hosts ]; then
  rm hosts
fi

hostlist=$(scontrol show hostname $SLURM_JOB_NODELIST)

for f in $hostlist
  do
    echo $f':'$SLURM_NTASKS_PER_NODE >> hosts
  done


mpirun -np $SLURM_NTASKS -f hosts ./Arepo param.txt 1

