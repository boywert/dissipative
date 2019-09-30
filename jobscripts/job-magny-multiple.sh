#!/bin/tcsh
#$ -S /bin/tcsh
#$ -j n
#$ -cwd
#$ -pe mvapich 16
##$ -q postprocess.q
#$ -m be
#$ -M volker.springel@h-its.org
#$ -N powerspec
#$ -binding linear:16
#$ -l h_rt=24:00:00
#

cat $PE_HOSTFILE

source /etc/profile.d/modules.csh
module load sge
module load mvapich2/gcc/64/1.5.1-qlc
module load hydra




@ I = 0

while ($I <= 56)

    mpistart -np $NSLOTS  ./Arepo param.txt 7 $I

    @ I += 1
end




 
