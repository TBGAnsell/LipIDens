#!/usr/bin/env bash

# Bash script to orient the protein along the principle axis and rotate by 90 degrees. Retruns orientated atomistic structure.

### Input order for script ###
#1 - protein_AT_full
#2 - protein_rotate

if [ "${1#*.}" = "pdb" ]; then input="${1%%.*}"; else echo "Check input file is in pdb format"; fi

gmx editconf -f ${input}.pdb -o ${input}_princ.pdb -princ -rotate ${2} -center 0 0 0 >& gmx_tmp<<EOF
0
EOF


rm gmx_tmp

echo ${input}_princ.pdb
