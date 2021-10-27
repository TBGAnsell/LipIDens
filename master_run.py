#!/usr/bin/env python3

import os
import protocol
import protocol.system_setup
import protocol.process_structure as pps
import protocol.simulation as ps
import protocol.test_PyLipID_cutoffs as lip_test
import protocol.run_PyLipID as lip_run
import protocol.screen_PyLipID_data as lip_screen

"""
Master script for step wise implementation of the protocol.

Author: T. Bertie Ansell

"""

### Full mdp files - not shortened ###

### Section 1: USER DEFINED VARIABLES ###
protein_AT_full='MscS_prot.pdb' # Atomistic .pdb file of protein
nprot = 1 # Number of homomeric protein chains, for heteromers use nprot=1
protein_shift=0 # Shift protein position in membrane to align TM region within the bilayer (value can also be negative).
protein_rotate='0 90 0' # Rotate the protein position to align TM region within the bilayer (angle in x y z)
boxsize='15,15,15' # CG simulation box size

save_dir="Test" # Save directory name

forcefield='martini_v2.2' # Currently compatible with martini_v2.0, martini_v2.1, martini_v2.2
membrane_composition='Gram neg. inner membrane' # Define membrane composition.
                                                # Can either select from predefined membrane compositions or
                                                # define using insane.py syntax.
                                                # Current predefined membrane compositions available:
                                                # ['Gram neg. inner membrane', 'Gram neg. outer membrane',
                                                # 'Plasma membrane', 'ER membrane', 'Raft-like microdomain', 'Simple']
                                                #
                                                # Example of custom bilayer composition (-u=upper leafelt, -l=lower leaflet):
                                                #membrane_composition='-u POPC:50 -u DOPC:50 -l POPE:30 -l CHOL:10 -l DOPE:60'

#membrane_composition='-u DPG1:30 -l LPPA:30 -l FPMG:10'
CG_simulation_time=4 # time in us - recomended to simulate for at least 5us per replicate, ideally 10-15us
replicates=2 # number of CG replicates
stride=10 # Skip every X no. frame during trajectory processing and running PyLipID
n_cores=16 # Number of CPU to use to run gromacs mdrun commands

##################
### Code Below ###
##################
protocol_path=os.path.dirname(protocol.__file__)
path=protocol.system_setup.setup(protocol_path, save_dir)

### Structure processing ###
protein_AT_full=pps.process_structure(protocol_path, protein_AT_full, protein_rotate)

### Setting up and performing CG simulations ###
python3_path, dssp_path = ps.get_py_paths(protocol_path)
ps.system_setup(protocol_path, path)
ps.fetch_CG_itp(forcefield, path)
bilayer=ps.bilayer_select(membrane_composition)
ps.top_header(forcefield, path)
ps.run_CG(protocol_path, protein_AT_full, protein_shift, bilayer, boxsize, replicates, python3_path, dssp_path, n_cores, path, CG_simulation_time)
### PAUSE POINT - run CG trajectories ###

## Process CG trajectories ###
ps.trjconv_CG(protocol_path, stride, replicates, path)

#################
#################
#################

### Section 2: USER DEFINED VARIABLES ###
lipid_atoms = None # all lipid atom/bead will be considered
contact_frames = 30  # will only plot data if the contact was formed over ${contact_frames} frames.
distance_threshold = 0.65

lower_cutoff = [0.4, 0.425, 0.45, 0.475, 0.5, 0.55] # list of lower cutoffs to test in exhasutive search
upper_cutoff = [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9] # list of upper cutoffs to test in exhaustive search
timeunit = "us"

##################
### Code Below ###
##################
### Selecting PyLipID input parameters ###
traj=lip_test.load_traj(path)
lip_list=lip_test.get_lipids(bilayer) ### Make None default - and also for save_dir (add as bilayer=None in function??/ check pylipid)
for lipid in lip_list:
    print(lipid)
    fig_dir=lip_test.set_lipid(path, lipid)
    distance_set = lip_test.compute_minimum_distance(traj, lipid, fig_dir, 1, lipid_atoms=lipid_atoms,
                                               contact_frames=contact_frames, distance_threshold=distance_threshold)
    lip_test.plot_PDF(distance_set, 1000, "{}/PyLipID_cutoff_test_{}/dist_distribut_{}.pdf".format(path, lipid, lipid), lipid)

    cutoff_list, trajfile_list, topfile_list = lip_test.exhaustive_search_setup(path, lower_cutoff, upper_cutoff, replicates)
    num_of_binding_sites, duration_avgs, num_of_contacting_residues = lip_test.test_cutoffs(
                                     cutoff_list, trajfile_list, topfile_list, lipid, lipid_atoms,
                                     nprot=nprot, stride=stride, save_dir="{}/PyLipID_cutoff_test_{}".format(path, lipid), timeunit=timeunit)
    lip_test.ex_data_process(path, lipid, num_of_binding_sites, duration_avgs, num_of_contacting_residues, cutoff_list)
    lip_test.graph(cutoff_list, [num_of_binding_sites[cutoffs] for cutoffs in cutoff_list],
          "num. of binding sites", lipid, f"{path}/PyLipID_cutoff_test_{lipid}/test_cutoff_num_of_bs_{lipid}.pdf")
    lip_test.graph(cutoff_list, [duration_avgs[cutoffs] for cutoffs in cutoff_list],
          f"Durations ({timeunit})", lipid, f"{path}/PyLipID_cutoff_test_{lipid}/test_cutoff_durations_{lipid}.pdf")
    lip_test.graph(cutoff_list, [num_of_contacting_residues[cutoffs] for cutoffs in cutoff_list],
          "num. of contacting residues", lipid,
          f"{path}/PyLipID_cutoff_test_{lipid}/test_cutoff_num_of_contacting_residues_{lipid}.pdf")

#################
#################
#################

### Section 3: USER DEFINED VARIABLES ###
cutoffs = [0.5, 0.7]  # dual-cutoff scheme for coarse-grained simulations. Single-cutoff scheme can be
                      # achieved by using the same value for two cutoffs.
dt_traj = None  # the timestep of trajectories. Need to use this param when trajectories are in a format
                # with no timestep information. Not necessary for trajectory formats of e.g. xtc, trr.

binding_site_size = 4  # binding site should contain at least four residues.
n_top_poses = 3     # write out num. of representative bound poses for each binding site.
n_clusters = "auto"  # cluster the bound poses for a binding site into num. of clusters. PyLipID
                     # will write out a pose conformation for each of the cluster. By default, i.e.
                     # "auto", PyLipID will use a density based clusterer to find possible clusters.
save_pose_format = "gro"  # format that poses are written in
save_pose_traj = True  # save all the bound poses in a trajectory for each binding site. The generated
                       # trajectories can take some disk space (up to a couple GB depending on your system).
save_pose_traj_format = "xtc"  # The format for the saved pose trajectories. Can take any format that is supported
                               # by mdtraj.

timeunit = "us"  # micro-sec. "ns" is nanosecond. Time unit used for reporting the results.
resi_offset = 0  # shift the residue index, useful for MARTINI models.

radii = None  # Radii of protein atoms/beads. In the format of python dictionary {atom_name: radius}
              # Used for calculation of binding site surface area. The van der waals radii of common atoms were
              # defined by mdtraj (https://github.com/mdtraj/mdtraj/blob/master/mdtraj/geometry/sasa.py#L56).
              # The radii of MARTINI 2.2 beads were included in PyLipID.

pdb_file_to_map = None   # if a pdb coordinate of the receptor is provided, a python script
                         # "show_binding_site_info.py" will be generated which maps the binding
                         # site information to the structure in PyMol. As PyMol cannot recognize
                         # coarse-grained structures, an atomistic structure of the receptor is needed.

fig_format = "pdf"  # format for all pylipid produced figures. Allow for formats that are supported by
                    # matplotlib.pyplot.savefig().

num_cpus = None  # the number of cpu to use when functions are using multiprocessing. By default,
                 # i.e. None, the functions will use up all the cpus available. This can use up all the memory in
                 # some cases.
##################
### Code Below ###
##################
### Running PyLipID analysis ###
trajfile_list, topfile_list=lip_run.get_trajectories(path, replicates)
for lipid in lip_list:
    lip_run.run_pylipid(trajfile_list, topfile_list, dt_traj, stride, lipid, lipid_atoms, cutoffs, nprot, binding_site_size,
        n_top_poses, n_clusters, save_dir, save_pose_format, save_pose_traj, save_pose_traj_format, timeunit, resi_offset,
         radii, pdb_file_to_map, fig_format, num_cpus)

#################
#################
#################

### Section 4: ###
##################
### Code Below ###
##################
### Screen PyLipID data ###
for lipid in lip_list:
    data = lip_screen.get_data(path, lipid)
    if data is not None:
        lip_screen.plot_screen_data(data, path, lipid)

#################
#################
#################

### Section 5: ###
### Rank sites ###

### Section 6 - User defined variables ###
input_CG_frame="md_fit_firstframe.gro"
model_type="de_novo" #option to specify whether to use 'aligned' or 'de_novo' structure from CG2AT for input of the AT simulations
replicates_AT=1
AT_simulation_time=100 #time in ns

##################
### Code Below ###
##################
### Setting up and running atomistic simulations ###
ps.CG2AT(protocol_path, protein_AT_full, input_CG_frame, save_dir)
AT_path=ps.system_setup_AT(protocol_path, path, model_type)
ps.run_AT(AT_path, replicates_AT, protocol_path, AT_simulation_time)
