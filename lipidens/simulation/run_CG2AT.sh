#!/usr/bin/env bash

# Bash script to run CG2AT using user suppleied CG frame and atomsitic protein structure.

### Input order for CG2AT script ###
#1- protein_AT_full
#2- input_CG_frame
#3- save_dir

cg2at=`which cg2at`
if [ -f "${cg2at}" ]; then
    ${cg2at} -h
    ${cg2at} -a $1 -c $2 -loc ${3}/CG2AT
else
    echo ${cg2at}
fi
