#!/bin/tcsh
#$ -j n
#$ -cwd
#$ -pe impi4 32
#$ -m be
#$ -M mpetkova@mpa-garching.mpg.de
#$ -N reion_cosmo
#$ -l h_rt=24:00:00

module load intel
module load impi
module load fftw/2.1.5
module load gsl/1.14

mpiexec -np $NSLOTS  ./Arepo  param.txt

if [ -f ../cont ] 
then
 rm -f "../cont"
 qsub job-opa-1.sh
fi
    





