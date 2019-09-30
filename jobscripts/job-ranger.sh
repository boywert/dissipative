#!/bin/tcsh
#$ -j n
#$ -V
#$ -o $JOB_NAME.o$JOB_ID
#$ -cwd
#$ -q  normal
#$ -pe 16way 128
#$ -m  be
#$ -M  volker@mpa-garching.mpg.de
#$ -N  str_128
#$ -l  h_rt=02:00:00


module load gmp/4.2.4
module load gsl/1.13
module load fftw2/2.1.5
module load hdf5/1.6.5

ibrun  ./Arepo   param.txt





