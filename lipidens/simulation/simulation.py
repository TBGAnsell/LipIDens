#!/usr/bin/env python3

import os
import subprocess
import re
import shutil
import sys
import urllib.request
import zipfile
import numpy as np

###############################
### Setup of CG simulations ###
###############################

def get_py_paths(protocol_path):
    """
    Establish python3 and dssp paths used as input for CG_simulation_setup.sh under run_CG(). Calls the bash script simulation/get_paths.sh.
    
    Params:
    -------
    protocol_path: str
        protocol path

    Returns:
    --------
    python3_path: str
        path to python3 environment 
    dssp_path: str
        path to dssp package
    martinize2_path: str
        path to martinize
    """
    python3_path=subprocess.check_output(['{}/simulation/get_paths.sh'.format(protocol_path), 'python']).decode(sys.stdout.encoding).strip()
    dssp_path=subprocess.check_output(['{}/simulation/get_paths.sh'.format(protocol_path), 'mkdssp']).decode(sys.stdout.encoding).strip()
    martinize2_path=subprocess.check_output(['{}/simulation/get_paths.sh'.format(protocol_path), 'martinize2']).decode(sys.stdout.encoding).strip()

    return python3_path, dssp_path, martinize2_path

def system_setup(protocol_path, path):
    """
    Establish directories and files needed for CG simulations.

    Params:
    -------
    protocol_path: str
        protocol path
    path: str
        path
    """
    print("Location of CG simulations:", path)

    os.makedirs(os.path.join(path, "itp_files"), exist_ok=True)
    shutil.copytree("{}/simulation/mdp_files".format(protocol_path), "{}/mdp_files".format(path), dirs_exist_ok=True)
    shutil.copytree("{}/simulation/python_files".format(protocol_path), "{}/python_files".format(path), dirs_exist_ok=True)

    return

def fetch_CG_itp(forcefield, path):
    """
    Get CG martini .itp files from martini website for Martini2.2.

    Params:
    -------
    forcefield: str
        select CG forcefield type
    path: str
        path
    """
    if forcefield in ["martini_v2.0", "martini_v2.1", "martini_v2.2"]:
        urllib.request.urlretrieve("https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini2/particle-definitions/{}.itp".format(forcefield), "{}/itp_files/{}.itp".format(path, forcefield))
        #urllib.request.urlretrieve("http://cgmartini.nl/images/parameters/lipids/Collections/martini_v2.0_lipids_all_201506.itp", "{}/itp_files/martini_v2.0_lipids.itp".format(path))
        shutil.copy("{}/simulation/martini_v2.0_lipids_all_201506.it".format(protocol_path),"{}/itp_files/martini_v2.0_lipids.itp".format(path))
        urllib.request.urlretrieve("https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini2/ions/martini_v2.0_ions.itp", "{}/itp_files/martini_v2.0_ions.itp".format(path))

    elif forcefield in ["martini_v3.0.0"]:
        urllib.request.urlretrieve("https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip", "{}/itp_files/martini_3.0.zip".format(path))
        with zipfile.ZipFile("{}/itp_files/martini_3.0.zip".format(path), 'r') as zip_dir:
            zip_dir.extractall("{}/itp_files/".format(path))

    else:
        print("Error with topology file download for the specified CG forcefield - please consider a different forcefield or download the required .itp files")


def bilayer_select(bilayer):
    """
    Select bilayer composition from either predefined selection (bilayer_type_dict) or inputted by the user in insane.py format.
    
    Params:
    -------
    bilayer: str
        input bilayer composition

    Returns:
    --------
    bilayer: str
        bilayer composition in insane format    
    """
    bilayer_type = {"Gram neg. inner membrane": "-u POPE:67 -u POPG:23 -u CDL2:10 -l POPE:67 -l POPG:23 -l CDL2:10",
    "Gram neg. outer membrane": "-u PGIN:100 -l POPE:90 -l POPG:5 -l CDL2:5",
    "Plasma membrane": "-u POPC:20 -u DOPC:20 -u POPE:5 -u DOPE:5 -u DPSM:15 -u DPG3:10 -u CHOL:25 -l POPC:5 -l DOPC:5 -l POPE:20 -l DOPE:20 -l POPS:8 -l DOPS:7 -l POP2:10 -l CHOL:25",
    "ER membrane": "-u POPC:37 -u DOPC:37 -u POPE:8 -u DOPE:8 -u CHOL:10 -l POPC:15 -l DOPC:15 -l POPE:20 -l DOPE:20 -l POPS:10 -l POP2:10 -l CHOL:10",
    "Raft-like microdomain": "-u DPPC:27 -u DPPE:8 -u DPSM:15 -u DPG3:10 -u CHOL:40 -l DPPC:15 -l DPPE:35 -l DPPS:10 -l CHOL:40",
    "Simple": "-u POPC:70 -u CHOL:30 -l POPC:70 -l CHOL:30"}

    if bilayer in bilayer_type.keys():
        print("\nPredefined bilayer selected: {}".format(bilayer))
        bilayer=bilayer_type[bilayer]
    elif bilayer is None:
        print("\nNo bilayer selected: please choose a membrane composition")
        quit()
    else:
        print("\nCustom bilayer selected: {}".format(bilayer))
    return bilayer

def top_header(forcefield, path):
    """
    Make topology file header in preparation for running the CG simulations.

    Params:
    -------
    forcefield: str
        CG forcefield type
    path: str
        path
    """
    itp_list=[i for i in os.listdir('{}/itp_files'.format(path)) if i.endswith('.itp')]
    itp_list.remove("{}.itp".format(forcefield))

    with open("{}/itp_files/top_header.txt".format(path), 'w+') as fn:
        fn.write('#include "{}.itp"\n'.format(forcefield))
        for itp in itp_list:
            fn.write('#include "{}"\n'.format(itp))


def run_CG(protocol_path, protein_AT_full, protein_shift, bilayer, boxsize, replicates, python3_path, dssp_path, n_cores, path, CG_simulation_time, martinize2_path, forcefield, martini_maxwarn, ring_lipids):
    """
    Run the setup of CG simulations calling the CG_simulation_setup.sh script.
    
    Params:
    -------
    protocol_path: str
        path to LipIDens
    protein_AT_full: str
        protein coordinate file
    protein_shift: float
        value to shift protein position in Z
    bilayer: str
        bilayer composition in insane format
    boxsize: str
        box size in X,Y,Z
    replicates: int
        number of replicates
    python3_path: str
        path to python3 environment
    dssp_path: str
        path to dssp package
    n_cores: int
        number of cores
    path: str
        path
    CG_simulation_time: int
        CG simulation time in us
    martinize2_path: str
        path to martinize2
    forcefield: str
        forcefield type
    martini_maxwarn: int
        maxwarn number for martinize2
    ring_lipids: bool
        ring lipids around/within the protein     
    """
    protein_AT_full=os.path.join(os.getcwd(), protein_AT_full)
    repstr=np.arange(1, replicates+1)
    repstr=' '.join(map(str, repstr))
    subprocess.check_call(["{}/simulation/CG_simulation_setup.sh".format(protocol_path), protein_AT_full, str(protein_shift), bilayer, boxsize, repstr, python3_path, dssp_path, str(n_cores), path, str(CG_simulation_time), martinize2_path, str(forcefield),str(martini_maxwarn), ring_lipids])
    for i in range(1, replicates+1):
        print(f"\n After initial setup attempt for run{i}:")
        if os.path.isfile(f"{path}/run{i}/protein_cg.gro"):
            print("\nmartinize (coarse-graining) was successful")
            if not os.path.isfile(f"{path}/run{i}/md.tpr"):
                print(f"run{i}/md.tpr not generated - retrying replicate simulation setup")
                subprocess.check_call(["{}/simulation/CG_simulation_setup.sh".format(protocol_path), protein_AT_full, str(protein_shift), bilayer, boxsize, str(i), python3_path, dssp_path, str(n_cores), path, str(CG_simulation_time), martinize2_path, str(forcefield),str(martini_maxwarn), ring_lipids])
                if not os.path.isfile(f"{path}/run{i}/md.tpr"):
                    print(f"\nsecond setup attempt for run{i} was not successful, please check the outputs")
                
        else:
            print(f"\nmartinize2 (coarse-graining) was not successful - please check the outputs for run{i}")
            print("\nMost likely there is an issue with your dssp path or version")


def trjconv_CG(protocol_path, stride, replicates, path):
    """
    Basic processing of CG simulations using trjconv commands in CG_simulation_process.sh.

    Params:
    -------
    protocol_path: str
        path to LipIDens
    stride: int
        skip X number of frames during trajectory processing
    replicates: int
        number of replicates
    path: str
        path
    """
    subprocess.check_call(["{}/simulation/CG_simulation_process.sh".format(protocol_path), str(stride), str(replicates), path])

###############################
### Setup of AT simulations ###
###############################

def CG2AT(protocol_path, protein_AT_full, input_CG_frame, path):
    """
    Run CG2AT using atomistic structure and provided frame from CG simulations.

    Params:
    -------
    protocol_path: str
        path to LipIDens
    protein_AT_full: str
        protein coordinate file
    input_CG_frame: str
        CG coordinate file for backmapping
    path: str
        path
    """
    print("\nRunning CG2AT\n")

    # check file locations for atomistic and CG inputs
    if os.path.isfile(f"{path}/{protein_AT_full}"):
        protein_AT_full=f"{path}/{protein_AT_full}"
    else:
        print(f"\n{path}/{protein_AT_full} not found.")
        a_name=str(input("Define absolute path to atomistic protein coordinate file: "))
        if os.path.isfile(f"{a_name}"):
            protein_AT_full=a_name
        else:
            print(f"\n{a_name} not found.")
            exit()
    if os.path.isfile(f"{path}/{input_CG_frame}"):
        input_CG_frame=f"{path}/{input_CG_frame}"
    else:
        print(f"\n{path}/{input_CG_frame} not found.")
        c_name=str(input("Define absolute path to CG input frame: "))
        if os.path.isfile(f"{c_name}"):
            input_CG_frame=c_name
        else:
            print(f"\n{c_name} not found.")
            exit()
    subprocess.check_call(["{}/simulation/run_CG2AT.sh".format(protocol_path), protein_AT_full, input_CG_frame, path])


def system_setup_AT(protocol_path, path, model_type):
    """
    Check whether CG2AT produced all files needed for AT simulations (.itp and .pdb). Establish directories and files needed for AT simulations.
    
    Params:
    -------
    protocol_path: str
        path to LipIDens
    path: str
        path
    model_type: str
        de nove or aligned structure

    Returns:
    --------
    AT_path: str
        path to atomistic simulations
    """
    AT_path=os.path.join(path, "Atomistic_Sims")
    os.makedirs(AT_path, exist_ok=True)

    print("Location of AT simulations:", AT_path)

    if os.path.isdir("{}/CG2AT".format(path)):
        #Check whether input PDB file exists
        if os.path.isfile("{}/CG2AT/FINAL/final_cg2at_{}.pdb".format(path, model_type)):
            print("CG2AT .pdb file found - {} structure".format(model_type))
            struc="{}/CG2AT/FINAL/final_cg2at_{}.pdb".format(path, model_type)
            shutil.copy(struc, AT_path)
        else:
            print("CG2AT .pdb file not found - check whether CG2AT was successful")
            exit()

        #Check whether itp files exist
        if any(fle.endswith('.itp') for fle in os.listdir("{}/CG2AT/FINAL".format(path))):
            print("CG2AT .itp files found")
            os.makedirs(os.path.join(AT_path, "itp_files_AT"), exist_ok=True)
            itp_list=[i for i in os.listdir("{}/CG2AT/FINAL".format(path)) if i.endswith('.itp')]
            for itp in itp_list:
                shutil.copy("{}/CG2AT/FINAL/{}".format(path,itp), "{}/itp_files_AT".format(AT_path))
        else:
            print("CG2AT .itp files not found - check whether CG2AT was successful")
            exit()

        #Check whether top file exist
        if any(fle.endswith(".top") for fle in os.listdir("{}/CG2AT/FINAL".format(path))):
            print("CG2AT .top file found")
            top_list=[t for t in os.listdir("{}/CG2AT/FINAL".format(path)) if t.endswith('.top')]
            for top in top_list:
                shutil.copy("{}/CG2AT/FINAL/{}".format(path, top), "{}".format(AT_path))
        else:
            print("CG2AT .top file not found - check whether CG2AT was successful")
            exit()
    else:
        print("CG2AT directory not found - was CG2AT successful")
        exit()


    os.makedirs(os.path.join(AT_path, "mdp_files_AT"), exist_ok=True)

    # Check which forcefield used with CG2AT and get corresponding md.mdp file
    forcefield_list=["amber", "charmm", "opls"]
    with open("{}/{}".format(AT_path, top_list[0]), "r+" ) as f:
        for line in f:
            for ff in forcefield_list:
                for match in re.finditer(ff, line):
                    mdp=("{}/simulation/mdp_files_AT/md_{}.mdp".format(protocol_path, ff))
                    break
    print("\nAtomistic .mdp file:", mdp)

    shutil.copy(mdp, "{}/mdp_files_AT".format(AT_path))

    return AT_path

def run_AT(AT_path, replicates_AT, protocol_path, AT_simulation_time):
    """
    Run the setup of AT simulations calling the AT_simulation_setup.sh script. Output is md.tpr file ready for simulation.
    
    Params:
    -------
    AT_path: str
        path to atomistic simulations
    replicates_AT: int
        number of atomistic replicates
    protocol_path: str
        path to LipIDens
    AT_simulation_time: int
        atomistic simulation time in ns
    """

    subprocess.check_call(["{}/simulation/AT_simulation_setup.sh".format(protocol_path), str(replicates_AT), AT_path, str(AT_simulation_time)])
    return
