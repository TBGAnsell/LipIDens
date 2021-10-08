#!/usr/bin/env bash

# Perform some basic post-processing of CG trajectories (make PBC whole and skip some frames)

### Input order for CG processing script ###
# 1-stride
# 2-replicates
# 3-path

cd ${3}
for i in `seq 1 ${2}`
do
    cd run$i
    echo ""
    if [[ ! -f md.gro ]]; then echo "Cannot locate md.gro file for run${i} - please run the simulation!"; exit 0; fi
    echo "Processing replicate: ${i}"
    mkdir -p output_files
    if [ -f md.part0001.xtc ]; then gmx trjcat -f md.part*.xtc -o md.xtc; fi
    if [ -f index.ndx ]; then cp index.ndx index_old.ndx; fi

    gmx make_ndx -f md.tpr -o index.ndx >& output_files/index_new <<EOF
q
EOF
    gmx trjconv -s md.tpr -f md.xtc -n index.ndx -o md_stride.xtc -pbc whole -dt ${1} -tu ns >& output_files/trjconv_whole <<EOF
0
EOF
    gmx trjconv -s md.tpr -f md_stride.xtc -o md_stride_firstframe.gro -dump 0 -n index.ndx >& output_files/trjconv_dump <<EOF
0
EOF
    if [[ -f md_stride_firstframe.gro ]]; then echo "Trajectory $i processed"; fi
    echo ""
    cd ../
done
