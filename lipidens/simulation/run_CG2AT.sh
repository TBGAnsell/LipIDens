#!/usr/bin/env bash

# Bash script to run CG2AT using user suppleied CG frame and atomsitic protein structure.
# Uses cg2at-lite (cg2at_lite) if installed, otherwise the original CG2AT (cg2at). Both prompt for the force field
# and fragment library.

### Input order for CG2AT script ###
#1- protein_AT_full
#2- input_CG_frame
#3- save_dir

if command -v cg2at_lite >/dev/null 2>&1; then
    cg2at=cg2at_lite
elif command -v cg2at >/dev/null 2>&1; then
    cg2at=cg2at
else
    echo "Neither cg2at_lite (cg2at-lite) nor cg2at (CG2AT) was found on the PATH" >&2
    exit 1
fi
${cg2at} -h
${cg2at} -a $1 -c $2 -loc ${3}/CG2AT
