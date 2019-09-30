#$ -S /bin/bash
#$ -j n
#$ -cwd
#$ -pe mvapich 32
#$ -q standard.q,maxmem.q
#$ -m be
#$ -M volker.springel@h-its.org
#$ -binding linear:32
#$ -l h_rt=24:00:00
#$ -N test 

source /etc/profile.d/modules.sh

module load sge
module load mvapich2/gcc/64/1.6-qlc
module load tap
module load tap_hdf5
module load tap_gmp
module load tap_hwloc
module load hydra


export OMP_NUM_THREADS=4
export NTASKS=$((NSLOTS/OMP_NUM_THREADS))
export TASKS_PER_HOST=$((NTASKS/NHOSTS))

mpistart -ppn $TASKS_PER_HOST -n $NTASKS  ./Arepo param.txt

