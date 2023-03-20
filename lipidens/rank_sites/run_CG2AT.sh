#!/usr/bin/env bash

# Bash script to run CG2AT using user suppleied CG frame and atomsitic protein structure.

### Input order for CG2AT script ###
#1- input_CG_frame
#2- save_dir
#3- CG_forcefield

cg2at=`which cg2at`
if [ -f "${cg2at}" ]; then
    ${cg2at} -c $1 -loc ${2}/CG2AT -w tip3p -o de_novo -ff charmm36 -fg ${3} 
else
    echo ${cg2at}
fi
