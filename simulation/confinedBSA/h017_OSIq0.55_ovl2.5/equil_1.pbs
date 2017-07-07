#!/bin/bash
#PBS -q batch
#PBS -A MAT049
#PBS -N q0.55_ovl2.5
#PBS -j oe
#PBS -l nodes=4
#PBS -l walltime=02:00:00
#PBS -l gres=atlas1%atlas2

export CRAY_CUDA_MPS=1
source $MODULESHOME/init/bash
module load gromacs/5.1.0

cd $PBS_O_WORKDIR

# continue the previous NPT run but changing options for the control run (now run in NVT)
oldtpr="relax_2.10.tpr"
oldcpt="relax_2.cpt"
aprun -n 1 gmx_mpi grompp -c $oldtpr -t $oldcpt -p confinedBSA.top -f equil_1.mdp -o equil_1.tpr

# 8 MPI ranks per node, each one tied to 2 openMP threads
aprun -n 32 -N 8 -d 2 gmx_mpi mdrun -ntomp 2 -gpu_id 00000000 -s equil_1.tpr -deffnm equil_1

wait
