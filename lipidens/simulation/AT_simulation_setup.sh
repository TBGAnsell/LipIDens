#!/usr/bin/env bash

### Input order for CG setup script ###
# 1-replicates_AT
# 2-AT_path
# 3-AT_simulation_time

cd ${2}
echo ""
echo "Setting up atomisitc simulations:"
for i in `seq 1 ${1}`
do
    echo ""
    echo "Setting up replicate: ${i}"
    mkdir -p run$i
    cd run$i
    mkdir -p output_files

    MDP=../mdp_files_AT
    simtime=$((${3}*500000))
    sed -i "s/XXX/${simtime%.*}/g" ${MDP}/md*.mdp

    # cp itp and top files to the dir
    cp ../itp_files_AT/*.itp .
    cp ../*.top .
    cp ../*.pdb system.pdb

    all_mols=`sed -ne '/molecules/,$ p' *.top | grep "^[^#;[]" | sort -u | awk '{print $1}' | grep -v '^..$'`
    echo $all_mols

    rem=( 'PROT' 'Prot' 'prot' 'SOL' 'Sol' 'sol' 'ION' 'Ion' 'ion' 'WAT' 'Wat' 'wat' 'SPC' 'TIP' )

    for i in "${rem[@]}"
    do
      all_mols="$( echo "$all_mols" | sed -e "s/${i}.*//g" )"
    done

    lip_list=`echo $all_mols | tr " " "\n" | perl -ne 'print "r${_}*|"' | awk '{print}' ORS=''`


    ### make index file ###
    gmx make_ndx -f system.pdb -o index.ndx >& output_files/index <<EOF
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


    # MD
    gmx grompp -f ${MDP}/md*.mdp -c system.pdb -r system.pdb -p *.top -n index.ndx -o md.tpr -maxwarn 1 >& output_files/md_grompp
    if [[ -f md.tpr ]]; then echo "md.tpr file generated"; fi
    if [[ ! -f md.tpr ]]; then echo "Failed to generate md.tpr, check system!"; exit 0; fi

    rm *#
    rm step*.pdb

    cd ../
done
