#!/bin/bash
if [ -z "$1" ];
then
  echo "need to pass a number as argument"
  exit 1
fi
sed -e "s/_PREV_/$1/g" equil_1.extend.pbs.tpl > equil_1.extend.pbs
qsub equil_1.extend.pbs

