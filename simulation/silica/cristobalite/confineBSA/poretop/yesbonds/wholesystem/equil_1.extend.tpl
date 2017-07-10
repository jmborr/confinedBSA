#!/bin/bash
#PBS -q batch
#PBS -A MAT049
#PBS -N q0.90_ovl1.5
#PBS -j oe
#PBS -l nodes=4
#PBS -l walltime=02:00:00
#PBS -l gres=atlas1%atlas2

export CRAY_CUDA_MPS=1
source $MODULESHOME/init/bash
module load gromacs/5.1.0

# iteration control
extratime=1000   # extend simulations by this much (in picoseconds)
last=10  # no more than these simulations
prev=_PREV_
curr=$(($prev + 1))

cd $PBS_O_WORKDIR
# Extend the run
aprun -n 1 gmx_mpi convert-tpr -s equil_1.${prev}.tpr -extend $extratime -o equil_1.${curr}.tpr
# 8 MPI ranks per node, each one tied to 2 openMP threads
aprun -n 32 -N 8 -d 2 gmx_mpi mdrun -ntomp 2 -gpu_id 00000000 -s equil_1.${curr}.tpr -cpi equil_1.cpt -deffnm equil_1 -append

# submit next job
if [ $curr -eq $last ];
then
  echo "last simulation done"
  exit
fi
./equil_1.submitnext.sh $curr  # submit next job

