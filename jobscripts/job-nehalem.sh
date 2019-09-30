#!/bin/sh
#$ -j n
#$ -cwd
#$ -pe mvapich 16
#$ -q nehalem-01
#$ -m be
#$ -M volker.springel@h-its.org
#$ -N cluster
#$ -l h_rt=24:00:00
#

cat $PE_HOSTFILE

. /etc/profile.d/modules.sh
module load sge
module load mvapich2/gcc/64/1.6
module load hydra

mpistart -np $NSLOTS  ./Arepo param.txt 


if [ -f "../cont" ] 
then
#    rm -f "../cont"
#    qsub job-magny.sh
fi





