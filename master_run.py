#!/usr/bin/env python3

import os
import lipidens
import lipidens.system_setup
import lipidens.process_structure as pps
import lipidens.simulation as ps
import lipidens.test_PyLipID_cutoffs as lip_test
import lipidens.run_PyLipID as lip_run
import lipidens.screen_PyLipID_data as lip_screen
import lipidens.rank_sites as rs

"""
Master python script for step wise implementation of the LipIDens protocol. The protocol embeds the PyLipID analysis toolkit for calculation of lipid binding sites
and their associated kinetics within a broader pipeline to assist interpretation of cryo-EM densities.
The LipIDens protocol uses the GROMACS (> version 5) simulation software.

The protocol code can be decomposed into numbered stages corresponding to:

1) Setting up and performing coarse-grained simulations
2) Testing PyLipID cut-off values
3) Selecting PyLipID input parameters and running PyLipID analysis
4) Screening PyLipID data
5) Ranking site lipids
6) Setting up atomistic simulations

Each subsection should be run independently by e.g. commenting out latter stages of the protocol. This is essential because, while the protocol establishes inputs needed
for coarse-grained or atomistic simulations it will not run the simulations. The user run simulations locally or offload to a high performance computing facility.
Sections corresponding to simulation runs are marked by ### PAUSE POINT ###. Once simulations have reached completion the pipeline can be resumed to analyse and process data.

Within each subsection the code can be decomposed into a list of USER DEFINED VARIABLES and the associated CODE for that section.

Author: T. Bertie Ansell, Wanling Song

"""

### Full mdp files - not shortened ###
##########################################################
### Section 1:Setting up and performing CG simulations ###
### USER DEFINED VARIABLES ###############################
##########################################################
protein_AT_full='TMD.pdb' # Atomistic .pdb file of protein. Protein residues should not have missing atoms however missing segments/loops are permitted.
nprot = 1 # Number of homomeric protein chains, for heteromers use nprot=1
protein_shift=0 # Shift protein position in membrane to align TM region within the bilayer (value can also be negative).
protein_rotate='0 90 0' # Rotate the protein position to align TM region within the bilayer (angle in x y z)
boxsize='15,15,15' # CG simulation box size

save_dir="Test" # Save directory name

forcefield='martini_v2.2' # Currently compatible with martini_v2.0, martini_v2.1, martini_v2.2, martini_v3.0.0
membrane_composition='Simple' # Define membrane composition.
                                                # Can either select from predefined membrane compositions or
                                                # define using insane.py syntax.
                                                # Current predefined membrane compositions available:
                                                # ['Gram neg. inner membrane', 'Gram neg. outer membrane',
                                                # 'Plasma membrane', 'ER membrane', 'Raft-like microdomain', 'Simple']
                                                #
                                                # Example of custom bilayer composition (-u=upper leaflet, -l=lower leaflet):
                                                #membrane_composition='-u POPC:50 -u DOPC:50 -l POPE:30 -l CHOL:10 -l DOPE:60'

#membrane_composition='-u POPC:100 -l DOPC:100'
CG_simulation_time=15 # time in us - recommended to simulate for at least 5us per replicate, ideally 10-15us
replicates=3 # number of CG replicates
stride=10 # Skip every X no. frame during trajectory processing and running PyLipID
n_cores=16 # Number of CPU to use to run gromacs mdrun commands

#############################
### Section 1: CODE Below ###
#############################
protocol_path=os.path.dirname(lipidens.__file__)
path=lipidens.system_setup.setup(protocol_path, save_dir)

### Structure processing ###
protein_AT_full=pps.process_structure(protocol_path, protein_AT_full, protein_rotate)

### Setting up and CG simulations ###
python3_path, dssp_path, martinize2_path = ps.get_py_paths(protocol_path)
bilayer=ps.bilayer_select(membrane_composition)

ps.system_setup(protocol_path, path)
ps.fetch_CG_itp(forcefield, path)
ps.top_header(forcefield, path)
ps.run_CG(protocol_path, protein_AT_full, protein_shift, bilayer, boxsize, replicates, python3_path, dssp_path, n_cores, path, CG_simulation_time, martinize2_path, forcefield)
#############################################
### PAUSE POINT - run the CG trajectories ###
#############################################

### Process CG trajectories ###
ps.trjconv_CG(protocol_path, stride, replicates, path)

################################################
### Section 2:Testing PyLipID cut-off values ###
### USER DEFINED VARIABLES #####################
################################################
lipid_atoms = None # all lipid atom/bead will be considered
contact_frames = 30  # will only plot data if the contact was formed over X number of frames where X=contact_frames.
distance_threshold = 0.65 # plot data only if lipid comes within distance_treshold of the protein

lower_cutoff = [0.4, 0.425, 0.45, 0.475, 0.5, 0.55] # list of lower cutoffs to test in exhaustive search
upper_cutoff = [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9] # list of upper cutoffs to test in exhaustive search
timeunit = "us"

#############################
### Section 2: CODE Below ###
#############################
### Testing PyLipID cut-offs  ###
lip_list=lip_test.get_lipids(bilayer)
traj=lip_test.load_traj(path)
for lipid in lip_list:
   print("\n Testing:", lipid)
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


##################################################################################
### Section 3: Selecting PyLipID input parameters and running PyLipID analysis ###
### USER DEFINED VARIABLES #######################################################
##################################################################################
cutoffs = [0.5, 0.7]  # dual-cutoff scheme for coarse-grained simulations. Single-cutoff scheme can be
                      # achieved by using the same value for two cutoffs.
lipid_atoms = None # all lipid atom/bead will be considered
dt_traj = None  # the timestep of trajectories. Need to use this param when trajectories are in a format
                # with no timestep information. Not necessary for trajectory formats of e.g. xtc, trr.

binding_site_size = 4  # binding site should contain at least four residues.
n_top_poses = 3     # write out num. of representative bound poses for each binding site.
n_clusters = "auto"  # cluster the bound poses for a binding site into num. of clusters. PyLipID
                    # will write out a pose conformation for each of the cluster. By default, i.e.
                    # "auto", PyLipID will use a density based clusterer to find possible clusters.
save_pose_format = "gro"  # format that poses are written in
save_pose_traj = True  # save all the bound poses in a trajectory for each binding site. The generated
                     #  trajectories can take some disk space (up to a couple GB depending on your system).
save_pose_traj_format = "xtc"  # The format for the saved pose trajectories. Can take any format that is supported
                      #         by mdtraj.

timeunit = "us"  # micro-sec. "ns" is nanosecond. Time unit used for reporting the results.
resi_offset = 0  # shift the residue index, useful for MARTINI models.

radii = None  # Radii of protein atoms/beads. In the format of python dictionary {atom_name: radius}
             #  Used for calculation of binding site surface area. The van der waals radii of common atoms were
             #  defined by mdtraj (https://github.com/mdtraj/mdtraj/blob/master/mdtraj/geometry/sasa.py#L56).
             # The radii of MARTINI 2.2 beads were included in PyLipID.

pdb_file_to_map = "TMD.pdb"   # if a pdb coordinate of the receptor is provided, a python script
                        # "show_binding_site_info.py" will be generated which maps the binding
                        # site information to the structure in PyMol. As PyMol cannot recognize
                        # coarse-grained structures, an atomistic structure of the receptor is needed.

fig_format = "pdf"  # format for all pylipid produced figures. Allow for formats that are supported by
                   # matplotlib.pyplot.savefig().

num_cpus = None  # the number of cpu to use when functions are using multiprocessing. By default,
                # i.e. None, the functions will use up all the cpus available. This can use up all the memory in
                # some cases.

#############################
### Section 3: CODE Below ###
#############################
### Running PyLipID analysis ###
trajfile_list, topfile_list=lip_run.get_trajectories(path, replicates)
for lipid in lip_list:
   lip_run.run_pylipid(trajfile_list, topfile_list, dt_traj, stride, lipid, lipid_atoms, cutoffs, nprot, binding_site_size,
       n_top_poses, n_clusters, save_dir, save_pose_format, save_pose_traj, save_pose_traj_format, timeunit, resi_offset,
        radii, pdb_file_to_map, fig_format, num_cpus)

#########################################
### Section 4: Screening PyLipID data ###
#########################################
### Section 4: CODE Below ###############
#########################################
### Screen PyLipID data ###
for lipid in lip_list:
   data = lip_screen.get_data(path, lipid)
   if data is not None:
       lip_screen.plot_screen_data(data, path, lipid)


######################################
### Section 5: Ranking site lipids ###
### USER DEFINED VARIABLES ###########
######################################
BindingSite_ID_dict={"POPC": [1, 2 , 5, 3],     # Dictionary of lipids (keys) and a list of binding site indices (values) to compare. The residence time of
          "POPE": [2, 3, 6, 1],                # sites are compared in the listed order e.g. POPC site 1 is compared with POPE site 2 and CHOL site 1.
          "CHOL":[1, 3, 5, "X"]}

#############################
### Section 5: CODE Below ###
#############################
### Ranking site lipids ###
rank_data=rs.get_BSstat(path, BindingSite_ID_dict)
rs.plot_site_rank(path, BindingSite_ID_dict, rank_data)

###################################################
### Section 6: Setting up atomistic simulations ###
### USER DEFINED VARIABLES ########################
###################################################
input_CG_frame="md_stride_t500.gro" # Input CG frame for backmapping to atomistic detail
model_type="aligned" # Option to specify whether to use 'aligned' or 'de_novo' structure from CG2AT for input of the AT simulations
replicates_AT=1 # Number of atomistic simulation replicates
AT_simulation_time=100 # Time in ns

#############################
### Section 6: CODE Below ###
#############################
### Setting up and running atomistic simulations ###
ps.CG2AT(protocol_path, protein_AT_full, input_CG_frame, save_dir)
AT_path=ps.system_setup_AT(protocol_path, path, model_type)
ps.run_AT(AT_path, replicates_AT, protocol_path, AT_simulation_time)

#############################################
### PAUSE POINT - run the AT trajectories ###
#############################################
