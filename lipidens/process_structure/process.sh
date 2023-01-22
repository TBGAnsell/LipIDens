#!/usr/bin/env bash

# Bash script to orient the protein along the principle axis and rotate by 90 degrees. Retruns orientated atomistic structure.

### Input order for script ###
#1 - protein_AT_full
#2 - protein_rotate


if [ "${1#*.}" = "pdb" ]; then input="${1%%.*}"; else echo "Check input file is in pdb format"; fi

pdb2pqr30 --ff CHARMM --keep-chain ${input}.pdb prot_pqr.pdb >& gmx_tmp

gmx pdb2gmx -f prot_pqr.pdb -ignh -ff charmm27 -water tip3p -o prot_sol.pdb >& gmx_tmp

gmx editconf -f prot_sol.pdb -o ${input}_princ.pdb -princ -d 8 -c -rotate ${2} >& gmx_tmp<<EOF
0
EOF

cat << EOF > em_short.mdp
integrator = steep
nsteps = 5000
emtol = 100
emstep = 0.001
EOF

gmx grompp -f em_short.mdp -c ${input}_princ.pdb -r ${input}_princ.pdb -o em.tpr -maxwarn 5 >& gmx_tmp

gmx mdrun -deffnm em -c prot_at_em_princ.pdb >& gmx_tmp

gmx trjconv -f prot_at_em_princ.pdb -o ${input}_princ.pdb -s em.tpr >& gmx_tmp<<EOF
0
EOF


rm gmx_tmp

echo ${input}_princ.pdb
