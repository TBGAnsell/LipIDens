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

Each subsection can be run independently by e.g. selected the appropriate stage of the protocol. This is essential because, while the protocol establishes inputs needed
for coarse-grained or atomistic simulations it will not run the simulations. The user can run simulations locally or offload to a high performance computing facility.
Sections corresponding to simulation runs are marked by ### PAUSE POINT ###. Once simulations have reached completion the pipeline can be resumed to analyse and process data.

To analyse simulations only users should store the trajectories under: "<Directory_name>/runX" for X number of replicates.

Within each subsection the code can be decomposed into a list of USER DEFINED VARIABLES and the associated CODE for that section.

Author: T. Bertie Ansell, Wanling Song

"""
##############################################
### Establish paths to be used in protocol ###
##############################################
protocol_path=os.path.dirname(lipidens.__file__)
save_dir=str(input("\nDirectory name (default: LipIDens_data): ") or "LipIDens_data")
path=lipidens.system_setup.setup(protocol_path, save_dir)
input_step=lipidens.system_setup.run_type()

### Full mdp files - not shortened ###
##########################################################
### Section 1:Setting up and performing CG simulations ###
### USER DEFINED VARIABLES ###############################
##########################################################
if input_step=="1a":
    print("\nInputs:\nAtomistic .pdb file of protein. Protein residues should not have missing atoms however missing segments/loops are permitted.")
    protein_AT_full=str(input("Protein .pdb file name: "))
    print("\nPRESS ENTER FOR DEFAULTS")
    nprot=int(input("\nNumber of homomeric protein chains, for heteromers input '1' (default: 1): ") or 1)
    print("\nShift protein position in membrane to align TM region within the bilayer (value can also be negative).")
    protein_shift=float(input("Protein shift e.g. 1.3 (default: 0): ") or 0)
    print("\nRotate the protein position to align TM region within the bilayer (angle in x y z)")
    protein_rotate=str(input("Protein rotate (default: 0 90 0): ") or '0 90 0')
    boxsize=str(input("Simulation box size (nm) (default: 15,15,15): ") or '15,15,15')

    forcefield_dict={"1": "martini_v2.0",
                     "2": "martini_v2.1",
                     "3": "martini_v2.2",
                     "4": "martini_v3.0.0"}

    print(*[f"\n{key}: {val}" for key, val in forcefield_dict.items()])
    forcefield=forcefield_dict[str(input("\nSelect forcefield:"))]

    print("\nDefine membrane composition.\nCan either select from predefined membrane compositions or define using insane.py syntax.")
    print("\nCurrent predefined membrane compositions available:\n['Gram neg. inner membrane', 'Gram neg. outer membrane','Plasma membrane', 'ER membrane', 'Raft-like microdomain', 'Simple']\n\nExample of custom bilayer composition (-u=upper leaflet, -l=lower leaflet):\n'-u POPC:50 -u DOPC:50 -l POPE:30 -l CHOL:10 -l DOPE:60'")
    membrane_composition=str(input("\nMembrane composition: "))
    martini_maxwarn=int(input("Add maxwarn number for martinize2 (default: 0): ") or 0)
    ring_lipids=str(input("Ring lipids around/within the protein, True/False (default=False): ") or False)
    CG_simulation_time=int(input("\nCoarse grain simulation time (in us) (default: 15): ") or 15)
    replicates=int(input("Number of replicates (default: 10): ") or 10)

    n_cores=int(input("Number of CPU to use (default: 16): ") or 16)


    #############################
    ### Section 1: CODE Below ###
    #############################
    ### Structure processing ###
    protein_AT_full=pps.process_structure(protocol_path, protein_AT_full, protein_rotate, path)
    bilayer=ps.bilayer_select(membrane_composition)

    ### Setting up and CG simulations ###
    python3_path, dssp_path, martinize2_path = ps.get_py_paths(protocol_path)
    ps.system_setup(protocol_path, path)
    ps.fetch_CG_itp(forcefield, path)
    ps.top_header(forcefield, path)
    ps.run_CG(protocol_path, protein_AT_full, protein_shift, bilayer, boxsize, replicates, python3_path, dssp_path, n_cores, path, CG_simulation_time, martinize2_path, forcefield, martini_maxwarn, ring_lipids)
#############################################
### PAUSE POINT - run the CG trajectories ###
#############################################

### Process CG trajectories ###
if input_step=="1b":
    replicates=int(input("Number of replicates to process: "))
    stride=int(input("Skip every X no. frames during trajectory processing (default: 10): ") or 10)
    ps.trjconv_CG(protocol_path, stride, replicates, path)

################################################
### Section 2:Testing PyLipID cut-off values ###
### USER DEFINED VARIABLES #####################
################################################
if input_step=="2":
    nprot=int(input("\nNumber of homomeric protein chains, for heteromers input '1' (default: 1): ") or 1)
    replicates=int(input("Number of replicates to analyse: "))
    #param_dict={1 :[variable, value, description, type]}
    param_dict={"1": ["lipid_atoms", None, "Lipid atoms to test (default: None uses all atoms)", "list_str"],
                "2": ["contact_frames", 30, "Only plot data if contact formed over <contact_frames> number of frames", "int"],
                "3": ["distance_threshold", 0.65, "Only plot data if lipid comes within <distance_treshold> of the protein (nm)", "float"],
                "4": ["lower_cutoff", [0.4, 0.425, 0.45, 0.475, 0.5, 0.55], "List of lower cut-offs to test in exhaustive search (nm)", "list_float"],
                "5": ["upper_cutoff", [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9], "List of upper cut-offs to test in exhaustive search (nm)", "list_float"],
                "6": ["timeunit", "us", "Time unit to plot", "str"],
                "7": ["stride", 10, "Skip every X no. frames", "int"]}
    param_dict=lipidens.system_setup.param_check(param_dict, "testing PyLipID cut-offs")
    lipid_atoms, contact_frames, distance_threshold, lower_cutoff, upper_cutoff, timeunit, stride=(param_dict[f"{x}"][1] for x in range(1, len(param_dict)+1))

    #############################
    ### Section 2: CODE Below ###
    #############################
    ### Testing PyLipID cut-offs  ###
    lip_list=lip_test.get_lipids(bilayer=None)
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
if input_step=="3":
    nprot=int(input("\nNumber of homomeric protein chains, for heteromers input '1' (default: 1): ") or 1)
    replicates=int(input("\nNumber of replicates to analyse: "))
    param_dict={"1": ["lipid_atoms", None, "Lipid atoms to test (default: None uses all atoms).", "list_str"],
                "2": ["dt_traj", None, "Need to change this when trajectories are in a format with no timestep information.\nNot necessary for trajectory formats of e.g. xtc, trr.", "str"],
                "3": ["binding_site_size", 4, "Binding site should contain at least <binding_site_size> residues.", "int"],
                "4": ["n_top_poses", 3, "Write out num. of representative bound poses for each binding site.", "int"],
                "5": ["n_clusters", "auto", "Cluster the bound poses for a binding site into <n_clusters> clusters.\nBy default, i.e. auto, PyLipID will use a density based clusterer to find all possible clusters.", "int"],
                "6": ["save_pose_format", "gro", "Pose coordinate file format - 'gro' or 'pdb'.", "str"],
                "7": ["save_pose_traj", True, "Save the bound poses to a trajectory for each binding site.", "bool"],
                "8": ["save_pose_traj_format", "xtc", "The format for the saved pose trajectories - any format supported by mdtraj.", "str"],
                "9": ["timeunit", "us", "Time unit for plots.", "str"],
                "10": ["resi_offset", 0, "Shift the residue index by <resi_offset>.", "int"],
                "11": ["radii", None, "Radii of protein atoms/beads in the format of python dictionary {atom_name: radius}. \nThe default (None) includes the van der waals radii of common atoms were defined by mdtraj (https://github.com/mdtraj/mdtraj/blob/master/mdtraj/geometry/sasa.py#L56). \nThe radii of MARTINI 2.2 beads were included in PyLipID.", "dict"],
                "12": ["fig_format", "pdf", "Save format for figures - all formats supported by matplotlib.pyplot.savefig().", "str"],
                "13": ["num_cpus", None, "The number of cpu to use when functions are using multiprocessing.\nThe default (None) will use all available cpus.", "int"],
                "14": ["stride", 10, "Skip every X no. frames", "int"]}
    param_dict=lipidens.system_setup.param_check(param_dict, "running PyLipID analysis")
    lipid_atoms, dt_traj, binding_site_size, n_top_poses, n_clusters, save_pose_format, save_pose_traj, save_pose_traj_format, timeunit, resi_offset, radii, fig_format, num_cpus, stride=(param_dict[f"{x}"][1] for x in range(1, len(param_dict)+1))

    pdb_file_to_map = str(input("\nAtomistic pdb coordinate file of protein. If provided, binding site information will be mapped to structure: "))

    #############################
    ### Section 3: CODE Below ###
    #############################
    ### Running PyLipID analysis ###
    trajfile_list, topfile_list=lip_run.get_trajectories(path, replicates)
    lip_list=lip_test.get_lipids(bilayer=None)
    for lipid in lip_list:
        print(f"\nIMPORTANT: Select cut-off scheme used to caluclate protein-{lipid} interactions.\nA single cut-off can be achieved using the same value for the lower and upper cutoff.")
        cutoffs= input("\nDual cut-off scheme in format lower, upper cut-off (nm) seperated by space e.g. '0.5 0.7': ").split()
        cutoffs=[float(x) for x in cutoffs]
        lip_run.run_pylipid(trajfile_list, topfile_list, dt_traj, stride, lipid, lipid_atoms, cutoffs, nprot, binding_site_size,
       n_top_poses, n_clusters, save_dir, save_pose_format, save_pose_traj, save_pose_traj_format, timeunit, resi_offset,
        radii, pdb_file_to_map, fig_format, num_cpus)

#########################################
### Section 4: Screening PyLipID data ###
#########################################
### Section 4: CODE Below ###############
#########################################
### Screen PyLipID data ###
if input_step=="4":
    lip_list=lip_test.get_lipids(bilayer=None)
    for lipid in lip_list:
        data = lip_screen.get_data(path, lipid)
        if data is not None:
            lip_screen.plot_screen_data(data, path, lipid)

######################################
### Section 5: Ranking site lipids ###
### USER DEFINED VARIABLES ###########
######################################
if input_step=="5":
    #############################
    ### Section 5: CODE Below ###
    #############################
    ### Ranking site lipids ###
    lip_list=lip_test.get_lipids(bilayer=None)
    BS_predict_dict=rs.compare_sites(path, lip_list)
    BindingSite_ID_dict=rs.get_site_compare(lip_list, BS_ID_dict=BS_predict_dict)
    rank_data=rs.get_BSstat(path, BindingSite_ID_dict)
    rs.plot_site_rank(path, BindingSite_ID_dict, rank_data)
    Pcorr, Dmap, Sfac=rs.locate_density_path(path)
    dens_path=rs.backmap_poses(path, protocol_path, BindingSite_ID_dict, save_dir)
    rs.pymol_density_compare(path, Pcorr, Dmap, Sfac, dens_path, BindingSite_ID_dict)




###################################################
### Section 6: Setting up atomistic simulations ###
### USER DEFINED VARIABLES ########################
###################################################
if input_step=="6":
    input_CG_frame=str(input("\nInput coarse-grained frame for backmapping to atomistic detail: "))
    protein_AT_full=str(input("Atomistic protein .pdb file name: "))
    print("\nProtein coordinates are backmapped either to the structure conformation (aligned) or those of the CG frame (de_novo).")
    model_type_dict={"1": "de_novo",
                     "2": "aligned"}

    print(*[f"\n{key}: {val}" for key, val in model_type_dict.items()])
    model_type=model_type_dict[str(input("\nSelect protein backmapping model type (default: de_novo):") or 1)]
    replicates_AT=int(input("Number of atomistic simulation replicates (default: 3): ") or 3)
    AT_simulation_time=int(input("\nAtomistic simulation time (in ns) (default: 100): ") or 100)

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
