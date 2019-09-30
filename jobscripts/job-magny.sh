#!/bin/sh
#$ -j n
#$ -cwd
#$ -pe mvapich 32
#$ -q standard.q,maxmem.q
#$ -m be
#$ -M volker.springel@h-its.org
#$ -N fof
#$ -binding linear:32
#$ -l h_rt=0:30:00
#

cat $PE_HOSTFILE

. /etc/profile.d/modules.sh
module load sge
module load mvapich2/gcc/64/1.6-qlc
module load tap
module load tap_hdf5
module load tap_gmp
module load tap_hwloc
module load hydra

mpistart -np $NSLOTS  ./Arepo param.txt 


if [ -f "../cont" ] 
then
#    rm -f "../cont"
#    qsub job-magny.sh
fi





