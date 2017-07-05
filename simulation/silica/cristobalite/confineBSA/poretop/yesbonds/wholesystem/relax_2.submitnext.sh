#!/bin/bash
if [ -z "$1" ];
then
  echo "need to pass a number as argument"
  exit 1
fi
sed -e "s/_PREV_/$1/g" relax_2.extend.tpl > relax_2.extend.pbs
qsub relax_2.extend.pbs

