#!/bin/bash
#PBS -A MAT049
#PBS -N run1
#PBS -j oe
#PBS -l nodes=48
#PBS -l walltime=02:00:00
#PBS -l gres=atlas1%atlas2

export CRAY_CUDA_MPS=1
source $MODULESHOME/init/bash
module load gromacs/5.1.0
cd $PROJWORK/mat049/kslater/nobonds

#aprun -n 1 gmx_mpi grompp -f confinedBSA_run1.mdp -c relax2_em.gro -p confinedBSA_em.top -o run1.tpr
aprun -n 384 -N 8 gmx_mpi mdrun -gpu_id 00000000 -deffnm run2 -cpi run1.cpt

