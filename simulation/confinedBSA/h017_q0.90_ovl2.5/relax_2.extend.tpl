#!/bin/bash
#PBS -q batch
#PBS -A MAT049
#PBS -N q0.90_ovl2.5
#PBS -j oe
#PBS -l nodes=4
#PBS -l walltime=02:00:00
#PBS -l gres=atlas1%atlas2

export CRAY_CUDA_MPS=1
source $MODULESHOME/init/bash
module load gromacs/5.1.0

# iteration control
extratime=1000   # extend simulations by this much (in picoseconds)
last=15  # no more than these simulations
prev=_PREV_
curr=$(($prev + 1))

cd $PBS_O_WORKDIR
# Extend the run
aprun -n 1 gmx_mpi convert-tpr -s relax_2.${prev}.tpr -extend $extratime -o relax_2.${curr}.tpr
# 8 MPI ranks per node, each one tied to 2 openMP threads
aprun -n 32 -N 8 -d 2 gmx_mpi mdrun -ntomp 2 -gpu_id 00000000 -s relax_2.${curr}.tpr -cpi relax_2.cpt -deffnm relax_2 -append

# submit next job
if [ $curr -eq $last ];
then
  echo "last simulation done"
  exit
fi
./relax_2.submitnext.sh $curr  # submit next job

