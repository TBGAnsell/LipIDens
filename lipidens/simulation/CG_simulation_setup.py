#!/usr/bin/env python3

import gromacs as gmx
import os

os.chdir(path)
print("Setting up coarse-grained simulations:")
for i in range(1, replicates+1):
    print("\n"+"Setting up replicate:"+i)
    os.mkdirs(path+"/run"+i)
    os.mkdirs(path+"/run"+i+"/output_files")
    os.chdir(path+"/run"+i)
    MDP="../mdp_files"
    simtime=CG_simulation_time*50000000




