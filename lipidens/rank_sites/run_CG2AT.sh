#!/usr/bin/env bash

# Bash script to run CG2AT using user suppleied CG frame and atomsitic protein structure.
# Uses cg2at-lite (cg2at_lite) if installed, otherwise the original CG2AT (cg2at).

### Input order for CG2AT script ###
#1- input_CG_frame
#2- save_dir
#3- CG_forcefield (fragment libraries for the original CG2AT, e.g. "martini_2-2_charmm36 martini_2-2_charmm36-VS")

if command -v cg2at_lite >/dev/null 2>&1; then
    # cg2at-lite provides the charmm36-cg2at force field and the fragment libraries without the -VS variants
    fragments=`echo $3 | tr ' ' '\n' | grep -v -- '-VS$' | tr '\n' ' '`
    cg2at_lite -c $1 -loc ${2}/CG2AT -w tip3p -o de_novo -ff charmm36-cg2at -fg ${fragments} -silent
elif command -v cg2at >/dev/null 2>&1; then
    cg2at -c $1 -loc ${2}/CG2AT -w tip3p -o de_novo -ff charmm36-jul2020-updated -fg ${3} -silent
else
    echo "Neither cg2at_lite (cg2at-lite) nor cg2at (CG2AT) was found on the PATH" >&2
    exit 1
fi
