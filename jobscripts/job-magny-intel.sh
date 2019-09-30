#$ -S /bin/bash
#$ -j n
#$ -cwd
#$ -pe impi48 96
#$ -q standard.q
#$ -m be
#$ -M volker.springel@h-its.org
#$ -N job 

source /etc/profile.d/modules.sh
module load sge
module load intel/impi

mpirun  -n $NSLOTS  ./MyExecutable
