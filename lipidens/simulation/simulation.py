#!/usr/bin/env python3

import os
import subprocess
import re
import shutil
import sys
import urllib.request
import zipfile

###############################
### Setup of CG simulations ###
###############################

def get_py_paths(protocol_path):
    """
    Establish python3 and dssp paths used as input for CG_simulation_setup.sh under run_CG(). Calls the bash script simulation/get_paths.sh.
    """
    python3_path=subprocess.check_output(['{}/simulation/get_paths.sh'.format(protocol_path), 'python']).decode(sys.stdout.encoding).strip()
    dssp_path=subprocess.check_output(['{}/simulation/get_paths.sh'.format(protocol_path), 'dssp']).decode(sys.stdout.encoding).strip()
    martinize2_path=subprocess.check_output(['{}/simulation/get_paths.sh'.format(protocol_path), 'martinize2']).decode(sys.stdout.encoding).strip()

    return python3_path, dssp_path, martinize2_path

def system_setup(protocol_path, path):
    """
    Establish directories and files needed for CG simulations.
    """
    print("Location of CG simulations:", path)

    os.makedirs(os.path.join(path, "itp_files"), exist_ok=True)
    shutil.copytree("{}/simulation/mdp_files".format(protocol_path), "{}/mdp_files".format(path), dirs_exist_ok=True)
    shutil.copytree("{}/simulation/python_files".format(protocol_path), "{}/python_files".format(path), dirs_exist_ok=True)

    return

def fetch_CG_itp(forcefield, path):
    """
    Get CG martini .itp files from martini website for Martini2.2.
    """
    if forcefield in ["martini_v2.0", "martini_v2.1", "martini_v2.2"]:
        urllib.request.urlretrieve("http://cgmartini.nl/images/parameters/ITP/{}.itp".format(forcefield), "{}/itp_files/{}.itp".format(path, forcefield))
        urllib.request.urlretrieve("http://cgmartini.nl/images/parameters/lipids/Collections/martini_v2.0_lipids_all_201506.itp", "{}/itp_files/martini_v2.0_lipids.itp".format(path))
        urllib.request.urlretrieve("http://cgmartini.nl/images/parameters/ITP/martini_v2.0_ions.itp", "{}/itp_files/martini_v2.0_ions.itp".format(path))

    elif forcefield in ["martini_v3.0.0"]:
        urllib.request.urlretrieve("http://cgmartini.nl/images/martini_v300.zip", "{}/itp_files/martini_3.0.zip".format(path))
        with zipfile.ZipFile("{}/itp_files/martini_3.0.zip".format(path), 'r') as zip_dir:
            zip_dir.extractall("{}/itp_files/".format(path))

    else:
        print("Error with topology file download for the specified CG forcefield - please consider a different forcefield or download the required .itp files")


def bilayer_select(bilayer):
    """
    Select bilayer composition from either predefined selection (bilayer_type_dict) or inputted by the user in insane.py format.
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
    Make topolopy file header in preperation for running the CG simulations.
    """
    itp_list=[i for i in os.listdir('{}/itp_files'.format(path)) if i.endswith('.itp')]
    itp_list.remove("{}.itp".format(forcefield))

    with open("{}/itp_files/top_header.txt".format(path), 'w+') as fn:
        fn.write('#include "{}.itp"\n'.format(forcefield))
        for itp in itp_list:
            fn.write('#include "{}"\n'.format(itp))


def run_CG(protocol_path, protein_AT_full, protein_shift, bilayer, boxsize, replicates, python3_path, dssp_path, n_cores, path, CG_simulation_time, martinize2_path, forcefield):
    """
    Run the setup of CG simulations calling the CG_simulation_setup.sh script.
    """
    protein_AT_full=os.path.join(os.getcwd(), protein_AT_full)
    subprocess.check_call(["{}/simulation/CG_simulation_setup.sh".format(protocol_path), protein_AT_full, str(protein_shift), bilayer, boxsize, str(replicates), python3_path, dssp_path, str(n_cores), path, str(CG_simulation_time), martinize2_path, str(forcefield)])

def trjconv_CG(protocol_path, stride, replicates, path):
    """
    Basic processing of CG simulations using trjconv commands in CG_simulation_process.sh.
    """
    subprocess.check_call(["{}/simulation/CG_simulation_process.sh".format(protocol_path), str(stride), str(replicates), path])

###############################
### Setup of AT simulations ###
###############################

def CG2AT(protocol_path, protein_AT_full, input_CG_frame, path):
    """
    Run CG2AT using atomsitic structure and provided frame from CG simulations.
    """
    print("\nRunning CG2AT\n")
    subprocess.check_call(["{}/simulation/run_CG2AT.sh".format(protocol_path), protein_AT_full, input_CG_frame, path])


def system_setup_AT(protocol_path, path, model_type):
    """
    Check whether CG2AT produced all files needed for AT simulations (.itp and .pdb). Establish directories and files needed for AT simulations.
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
    """

    subprocess.check_call(["{}/simulation/AT_simulation_setup.sh".format(protocol_path), str(replicates_AT), AT_path, str(AT_simulation_time)])
    return
