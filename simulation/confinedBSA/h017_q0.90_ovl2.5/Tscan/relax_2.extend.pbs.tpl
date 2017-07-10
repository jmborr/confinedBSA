#!/bin/bash
#PBS -q batch
#PBS -A MAT049
#PBS -N relax_2
#PBS -j oe
#PBS -l nodes=80
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

for rawt in $(seq 40 20 420); do
  t=`printf T%03d $rawt`
  cd $PBS_O_WORKDIR/$t/
  aprun -n 1 gmx_mpi convert-tpr -s relax_2.${prev}.tpr -extend $extratime -o relax_2.${curr}.tpr &
done
wait

for rawt in $(seq 40 20 420); do
  t=`printf T%03d $rawt`
  cd $PBS_O_WORKDIR/$t/
  # 8 MPI ranks per node, each one tied to 2 openMP threads
  aprun -n 32 -N 8 -d 2 gmx_mpi mdrun -ntomp 2 -gpu_id 00000000 -s relax_2.${curr}.tpr -cpi relax_2.cpt -deffnm relax_2 -append &
done
wait

# submit next job
if [ $curr -eq $last ];
then
  echo "last simulation done"
  exit
fi
cd $PBS_O_WORKDIR/
./relax_2.submitnext.sh $curr  # submit next job

