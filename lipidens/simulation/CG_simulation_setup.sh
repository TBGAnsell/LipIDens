#!/usr/bin/env bash

# Setup coarse-grained simulations of membrane protein system using orientated atomsitic input.

### Input order for CG setup script ###
# 1-protein_AT_full
# 2-protein_shift
# 3-bilayer
# 4-boxsize
# 5-replicates
# 6-python3_path
# 7-dssp_path
# 8-str(n_cores)
# 9-path
# 10-CG_simulation_time
# 11-martinize2 path
# 12-forcefield
# 13-martini_maxwarn 

cd ${9}
echo ""
echo "Setting up coarse-grained simulations:"
for i in `seq 1 ${5}`
do
    echo ""
    echo "Setting up replicate: ${i}"
    mkdir -p run$i
    cd run$i
    mkdir -p output_files

    MDP=../mdp_files
    simtime=$((${10}*50000000))
    sed -i "s/XXX/${simtime%.*}/g" ${MDP}/md.mdp

    # Change defaut lincs_iter and lincs_order params if CHOL with virtual sites is present in the membrane
    if [[ ${3} == *"CHOL"* ]] && [ ${i} == 1 ]; then lin=$'lincs_iter   =     2 \nlincs_order   =     12' ; echo "$lin" >> ${MDP}/md.mdp ; fi

    # cp itp files to the dir
    cp ../itp_files/*.itp .

    # OLD Martinize for martini 2 #
    #${6} ../python_files/martinize_gmx2019.py -f ${1} -o system.top -x protein_cg.gro -dssp ${7} -p backbone -ff elnedyn22 -cys auto -ef 1000 -el 0.5 -eu 0.9 -ea 0 -ep 0 >& output_files/martinize

    if [ ${12} != "martini_v3.0.0" ]; then
      ${11} -f ${1} -o system.top -x protein_cg.gro -dssp ${7} -p backbone -ff elnedyn22 -cys auto -ef 1000 -el 0.5 -eu 0.9 -ea 0 -ep 0 -maxwarn ${13} >& output_files/martinize
    else
      ${11} -f ${1} -o system.top -x protein_cg.gro -dssp ${7} -p backbone -ff martini3001 -cys auto -ef 1000 -el 0.5 -eu 0.9 -ea 0 -ep 0 -maxwarn ${13} >& output_files/martinize
    fi

    ${6} ../python_files/insane_mod.py -f protein_cg.gro -o system.gro  -pbc square -box ${4} -center -dm ${2} ${3} -sol W -p tmp.top

    sed -i '2d' system.top
    cat ../itp_files/top_header.txt system.top >out && mv out system.top
    echo "" >> system.top

    sed -n '10,$p' tmp.top >> system.top

    gmx grompp -f ${MDP}/em_steep_1.mdp -c system.gro -r system.gro -p system.top -o ions.tpr -maxwarn 1 >& output_files/ion_grompp
    gmx genion -s ions.tpr -o system_ions.gro -p system.top -pname NA -nname CL -conc 0.15 -neutral >& output_files/genion <<EOF
W
EOF
    if [ ${12} != "martini_v3.0.0" ]; then
      sed 's/\(NA \)\([ \t]*\)\( NA\)/NA+\2\NA+/g' system_ions.gro | sed "s/\(CL \)\([ \t]*\)\( CL\)/CL-\2\CL-/g" > system_ions_corrected.gro
      sed 's/NA  /NA+/g' system.top | sed 's/CL  /CL-/g' > system_corrected.top
    else
      cp system_ions.gro system_ions_corrected.gro
      cp system.top system_corrected.top
    fi
    if [[ -f system_ions_corrected.gro ]]; then echo "IMPORTANT: Please check you are happy with the protein orientation and bilayer position in system_ions_corrected.gro. If not then adjust the protein_shift and protein_rotate inputs and run again."; fi
    echo ""
    lip_list=`echo $3 | grep -Eoi "[a-z/]{2,}" | sort -u | perl -ne 'print "r${_}*|"' | awk '{print}' ORS=''`

    ### make index file ###
    gmx make_ndx -f system_ions_corrected.gro -o index.ndx >& output_files/index <<EOF
del 0
del 1-40
0|${lip_list}
1&!0
!1
name 1 Protein_Lipids
name 2 Lipids
name 3 Solvent
q
EOF

    # EM
    gmx grompp -f ${MDP}/em_steep_1.mdp -c system_ions_corrected.gro -r system_ions_corrected.gro -p system_corrected.top -n index.ndx -o em_steep_1.tpr >& output_files/em_grompp
    if [[ ! -f em_steep_1.tpr ]]; then echo "Failed to generate em_steep_1.tpr - check system!"; exit 0; fi
    gmx mdrun -deffnm em_steep_1 -v >& output_files/em
    if [[ ! -f em_steep_1.gro ]]; then echo "Energy minimisation failed!"; exit 0; fi

    # EQ
    echo "Running equilibration-1"
    gmx grompp -f ${MDP}/eq.1.mdp -c em_steep_1.gro -r em_steep_1.gro -p system_corrected.top -n index.ndx -o eq.1.tpr -maxwarn 1 >& output_files/eq1_grompp
    gmx mdrun -deffnm eq.1 -ntmpi 1 -ntomp ${8} -pin on -pinoffset 0 >& output_files/eq1
    if [[ ! -f eq.1.gro ]]; then echo "Equilibration-1 failed! Try increasing the box dimensions (using the 'boxsize' variable) or adjusting the protein orientation within the bilayer ('protein_shift' and 'protein_rotate' variables)."; exit 0; fi

    echo "Running equilibration-2"
    gmx grompp -f ${MDP}/eq.2.mdp -c eq.1.gro -r eq.1.gro -p system_corrected.top -n index.ndx -o eq.2.tpr -maxwarn 1 >& output_files/eq2_grompp
    gmx mdrun -deffnm eq.2 -ntmpi 1 -ntomp ${8} -pin on -pinoffset 0 >& output_files/eq2
    if [[ ! -f eq.2.gro ]]; then echo "Equilibration-2 failed! Try increasing the box dimensions (using the 'boxsize' variable) or adjusting the protein orientation within the bilayer ('protein_shift' and 'protein_rotate' variables)."; exit 0; fi

    # MD
    gmx grompp -f ${MDP}/md.mdp -c eq.2.gro -r eq.2.gro -p system_corrected.top -n index.ndx -o md.tpr -maxwarn 1 >& output_files/md_grompp
    if [[ -f md.tpr ]]; then echo "md.tpr file generated"; fi

    rm *#
    rm step*.pdb

    cd ../
done
