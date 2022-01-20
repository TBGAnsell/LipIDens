#!/usr/bin/env python3

"""
Script to plot the probability distribution of protein-lipid distances and test a range of lower and upper cutoff combinations.

Author: Wanling Song, Modified: T. Bertie Ansell

"""
#%matplotlib inline
import seaborn as sns
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pylipid
from pylipid.api import LipidInteraction
print(pylipid.__version__)

from itertools import product
import mdtraj as md
from pylipid.util import get_traj_info, check_dir
import matplotlib.ticker as ticker
from itertools import product
import shutil
import pickle

def get_lipids(bilayer):
    """
    Extract lipid names from inputted bilayer.
    """
    if bilayer !=None:
        lip_list=list(set(re.findall(r'\w[A-Z0-9]{2,}', bilayer)))
    elif bilayer==None:
        bilayer=str(input("\nEnter lipid names for analysis seperated by space e.g. POPC DOPC: "))
        lip_list=list(set(re.findall(r'\w[A-Z0-9]{2,}', bilayer)))
    else:
        print("Bilayer not found, please define the bilayer composition")
        exit()
    print("Lipids to test:", lip_list)
    return lip_list

def load_traj(path):
    """
    Load single trajectory used for obtaining lipid probability distribution.
    """
    xtc_def="md_stride.xtc"
    top_def="md_stride_firstframe.gro"

    print("\nLoading processed trajectories for cutoff testing:")
    try:
        if os.path.isfile("{}/run1/{}".format(path, xtc_def)):
            trajfile = "{}/run1/{}".format(path, xtc_def)
        else:
            print("\n{}/run1/{} not found.\nHave you processed the CG trajectories?".format(path, xtc_def))
            r_alt=str(input("Do you wish to define alternative 'xtc' file name? (y/n): "))
            if r_alt=="y":
                xtc_def=str(input("\nDefine alternative 'xtc' file name: "))
                if os.path.isfile("{}/run1/{}".format(path, xtc_def)):
                    trajfile="{}/run1/{}".format(path, xtc_def)
                else:
                    print("\n{}/run1/{} not found.".format(path, xtc_def))
                    exit()
            elif r_alt=="n":
                print("Re-run protocol to process CG trajectories")
                exit()
            else:
                print("INVALID: must enter y/n")
                exit()
        if os.path.isfile("{}/run1/{}".format(path, top_def)):
            topfile =  "{}/run1/{}".format(path, top_def)
        else:
            print("\n{}/run1/{} not found.\nHave you processed the CG trajectories?".format(path, top_def))
            r_alt=str(input("Do you wish to define alternative topology file? (y/n): "))
            if r_alt=="y":
                top_def=str(input("\nDefine alternative topology file name: "))
                if os.path.isfile("{}/run1/{}".format(path, top_def)):
                    topfile =  "{}/run1/{}".format(path, top_def)
                else:
                    print("\n{}/run1/{} not found.".format(path, top_def))
                    exit()
            elif r_alt=="n":
                print("Re-run protocol to process CG trajectories")
                exit()
            else:
                print("INVALID: must enter y/n")
                exit()
        traj = md.load(trajfile, top=topfile)
    except Exception as e:
        print(e)
        exit()
    return traj

def set_lipid(path, lipid):
    """
    Establish save directory of cutoff test data.
    """
    save_dir = "{}/PyLipID_cutoff_test_{}".format(path, lipid)
    fig_dir = check_dir(save_dir, "Figures", print_info=False)
    return fig_dir

def plot_minimum_distances(distances, times, title, fn):
    """
    Plot the per residue minimum distance to individual lipids.
    """
    fig, ax = plt.subplots(1, 1, figsize=(3, 2.5))
    ax.plot(times, distances)
    ax.set_xlabel(r"Time ($\mu$s)")
    ax.set_ylabel("Minimum distances (nm)")
    ax.set_title(title)
    ax.set_ylim(0, 1.0)
    sns.despine(top=True, right=True)
    plt.tight_layout()
    plt.savefig(fn,format='pdf', dpi=200)
    plt.close()
    return

def compute_minimum_distance(traj, lipid, fig_dir, nprot, lipid_atoms=None,
                            contact_frames=10, distance_threshold=0.65):
    """
    Obtain minimum distances of specified lipid to each residue in the protein if contact comes within distance_threshold
    for longer than the number of contact_frames.
    """
    DIST_CONTACT_ALL = []
    traj_info, _, _ = get_traj_info(traj, lipid, lipid_atoms=lipid_atoms)
    for protein_idx in np.arange(nprot, dtype=int):
        for residue_idx, residue_atom_indices in enumerate(
            traj_info["protein_residue_atomid_list"][protein_idx]):
            dist_matrix = np.array([np.min(
                            md.compute_distances(traj, np.array(list(product(residue_atom_indices, lipid_atom_indices)))),
                            axis=1) for lipid_atom_indices in traj_info["lipid_residue_atomid_list"]])
            # plot distances
            for lipid_idx in np.arange(len(dist_matrix)):
                if sum(dist_matrix[lipid_idx] < distance_threshold) >= contact_frames:
                    DIST_CONTACT_ALL.append(dist_matrix[lipid_idx])
                    plot_minimum_distances(dist_matrix[lipid_idx], traj.time/1000000.0,
                                           "{}-{}{}".format(traj_info["residue_list"][residue_idx], lipid, lipid_idx),
                                           "{}/dist_{}_{}{}.pdf".format(fig_dir, traj_info["residue_list"][residue_idx], lipid, lipid_idx))

    distance_set = np.concatenate(DIST_CONTACT_ALL)
    return distance_set

def plot_PDF(distance_set, num_of_bins, fn, lipid):
    """
    Plot the probability distribution of minimum lipid distances.
    """
    fig, ax = plt.subplots(1,1, figsize=(4,3))
    ax.hist(distance_set, bins=num_of_bins, density=True, color='lightcoral')
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Minimum distance (nm)")
    ax.set_ylabel("Probablity Density")
    ax.set_title(lipid)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    sns.despine(top=True, right=True)
    plt.tight_layout()
    plt.savefig(fn,format='pdf', dpi=200)
    return

def test_cutoffs(cutoff_list, trajfile_list, topfile_list, lipid, lipid_atoms, nprot=1,
                 stride=1, save_dir=None, timeunit="us"):
    """
    Perform exhastive cutoff testing by caluclating the number of binding sites, average duration and number of contacting residues
    for each pair of cutoffs in the cutoff list.
    """
    num_of_binding_sites = {}
    duration_avgs = {}
    num_of_contacting_residues = {}
    for cutoffs in cutoff_list:
        print("\n Testing cutoff pair:", cutoffs, "\n")
        li = LipidInteraction(trajfile_list, topfile_list=topfile_list, cutoffs=cutoffs, lipid=lipid,
                              lipid_atoms=lipid_atoms, nprot=1, timeunit=timeunit,
                              save_dir=save_dir, stride=stride)
        li.collect_residue_contacts()
        li.compute_residue_duration()
        li.compute_binding_nodes(print_data=False) # switch print and write to False for cleaner output.
        num_of_binding_sites[cutoffs] = len(li.node_list)
        duration_avgs[cutoffs] = li.dataset["Duration"].mean()
        num_of_contacting_residues[cutoffs] = sum(li.dataset["Duration"]>0)
    shutil.rmtree(li.save_dir)
    return num_of_binding_sites, duration_avgs, num_of_contacting_residues


def exhaustive_search_setup(path, lower_cutoff, upper_cutoff, replicates):
    """
    Obtain list of cutoff pairs to use for exhastive cutoff testing using user specified lower and upper cutoff lists.
    Load all coarse-grain trajectories to test.
    """
    print("\nInitiating exhastive cutoff search:\n")
    print("Lower cutoffs to test:", lower_cutoff)
    print("Upper cutoffs to test:", upper_cutoff, "\n")
    cutoff_list = list(product(lower_cutoff, upper_cutoff))
    trajfile_list=[]
    topfile_list=[]

    xtc_def="md_stride.xtc"
    top_def="md_stride_firstframe.gro"

    for n in range(1,replicates+1):
        try:
            if os.path.isfile("{}/run{}/{}".format(path, n, xtc_def)):
                trajfile="{}/run{}/{}".format(path, n, xtc_def)
            else:
                print("\n{}/run{}/{} not found.".format(path, n, xtc_def))
                xtc_def=str(input("Define alternative 'xtc' file name: "))
                if os.path.isfile("{}/run{}/{}".format(path, n, xtc_def)):
                    trajfile="{}/run{}/{}".format(path, n, xtc_def)
                else:
                    print("\n{}/run{}/{} not found.".format(path, n, xtc_def))
                    exit()
            if os.path.isfile("{}/run{}/{}".format(path, n, top_def)):
                topfile="{}/run{}/{}".format(path, n, top_def)
            else:
                print("\n{}/run{}/{} not found.".format(path, n, top_def))
                top_def=str(input("Define alternative topology file name: "))
                if os.path.isfile("{}/run{}/{}".format(path, n, top_def)):
                    topfile="{}/run{}/{}".format(path, n, top_def)
                else:
                    print("\n{}/run{}/{} not found.".format(path, n, top_def))
                    exit()
        
            trajfile_list.append(trajfile)
            topfile_list.append(topfile)
        except Exception as e:
            print(e)
            exit()

    print("List of trajectories to test:", trajfile_list, "\n")
    return cutoff_list, trajfile_list, topfile_list

def ex_data_process(path, lipid, num_of_binding_sites, duration_avgs, num_of_contacting_residues, cutoff_list):
    """
    Process output data from PyLipID (num_binding_sites, duration_avgs, num_of_contacting_residues). Save in pickle format.
    """
    test_data = {"num_of_binding_sites": num_of_binding_sites,
             "duration_avgs": duration_avgs,
             "num_of_contacting_residues": num_of_contacting_residues,
             "test_cutoff_list": cutoff_list}
    with open(f"{path}/PyLipID_cutoff_test_{lipid}/test_cutoff_data_{lipid}.pickle", "wb") as f:
        pickle.dump(test_data, f, 2)

def graph(cutoff_list, metric_values, ylabel, title, fn):
    """
    Plot the data from the exhastive cutoff testing.
    """
    fig, ax = plt.subplots(1, 1, figsize=(len(cutoff_list)*0.42, 3.6))
    ax.scatter(np.arange(len(cutoff_list)), metric_values, s=50, color='lightcoral')
    ax.set_xticks(np.arange(len(cutoff_list)))
    ax.set_xticklabels(cutoff_list, rotation=45, ha='right')
    ax.set_xlabel("Dual cutoff")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    sns.despine(top=True, right=True)
    plt.tight_layout()
    plt.savefig(fn, format='pdf',dpi=200)
    return
