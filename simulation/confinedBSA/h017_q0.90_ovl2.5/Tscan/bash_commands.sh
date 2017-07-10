#!/bin/bash

rootd=`pwd`

# temperatures
for rawt in $(seq 40 20 420); do
  t=`printf T%03d $rawt`; echo $t
  #mkdir -p $t
  #/bin/cp -r ../confinedBSA.top ../SIL.itp ../oplsaa.ff $t/
  ###############################
  # PREPARE JOB relax_2
  ###############################
  #sed -e "s/_TEMP_/$rawt/g" relax_2.mdp.tpl > $t/relax_2.mdp
  #cd $rootd/$t/ ; ln -s ../../relax_2.gro relax_1.gro  #use relaxed NPT at 300K as starting conformation
  ###############################
  # PREPARE JOB extension to relax_2
  ###############################
  #cd $rootd/$t/;  ln -s relax_2.tpr relax_2.0.tpr 
  ###############################
  # PREPARE JOB equil_1
  ###############################
  #sed -e "s/_TEMP_/$rawt/g" equil_1.mdp.tpl > $t/equil_1.mdp
  ###############################
  # PREPARE JOB extension to equil_1
  ###############################
  cd $rootd/$t/;  ln -s equil_1.tpr equil_1.0.tpr 
done