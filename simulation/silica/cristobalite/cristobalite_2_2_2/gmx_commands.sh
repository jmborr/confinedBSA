#!/bin/sh

# ********************
# minimization run
# ********************
gmx grompp -f minim.mdp -c cristobalite_2_2_2.gro -p cristobalite_2_2_2.top -o minim.tpr
gmx mdrun -ntmpi 4 -ntomp 1 -pin on -deffnm minim

# ********************
# relaxation NVT run
# ********************
gmx grompp -f relax_1.mdp -c minim.gro -p cristobalite_2_2_2.top -o relax_1.tpr
gmx mdrun -ntmpi 4 -ntomp 1 -pin on -deffnm relax_1

# ********************
# relaxation NPT run
# ********************
gmx grompp -f relax_2.mdp -c minim.gro -p cristobalite_2_2_2.top -o relax_2.tpr
gmx mdrun -ntmpi 4 -ntomp 1 -pin on -deffnm relax_2
echo "1 2 3 4 11 13 17 0" | gmx energy -f relax_2.edr -o relax_2.xvg
# **********************
# continuation NPT run
# **********************
# Use checkpoint relax_2.cpt for continuation
gmx grompp -f relax_3.mdp -t relax_2.cpt -c relax_2.tpr -o relax_3.tpr
gmx mdrun -ntmpi 4 -ntomp 1 -pin on -deffnm relax_3
echo "1 2 3 4 11 13 17 0" | gmx energy -f relax_3.edr -o relax_3.xvg
echo "0 0" | gmx msd -f relax_3.trr -s relax_2.gro -o relax_3_msd.xvg -rmcomm