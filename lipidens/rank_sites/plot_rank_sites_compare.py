#!/usr/bin/env python3
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import operator
import subprocess
import shutil
import numpy as np
import re
import warnings
warnings.filterwarnings( "ignore")

"""

Script to rank the residence times of different lipids binding to the same site using outputs from PyLipID.

Author: TBG Ansell, Wanling Song

"""

def compare_sites(path, lip_list):
    """
    Create dictionary of predicted comparable binding sites for each lipid, against the reference lipid (first entry).
    Comparable sites selected based on best match between comprising residues.

    Params
    ------
    path: str
        path
    lip_list: list
        lipid to compare

    Returns
    -------
    BS_ID_dict: dict
        predicted sites to compare as keys=lipid: values=binding sites in comparable order
    """
    if len(lip_list)>1:
        try:
            # Compare lipid sites based on comprising residues. Comparison made for each lipid against the reference lipid (first lipid in list).
            match_sites=pd.DataFrame()

            ref_lipid=lip_list[0]
            print("\nIMPORTANT: First lipid entry used as reference - not reccomended to use sterols as the reference lipid as sites likely to differ from phospholipids.")
            print('\nReference lipid '+ref_lipid)
            ref_csv=pd.read_csv(f"{path}/Interaction_{ref_lipid}/Dataset_{ref_lipid}/Dataset.csv")
            BS_ref=ref_csv["Binding Site ID"].unique()

            for lipid in lip_list[1:]:
                tmp_df=pd.DataFrame()
                print('Comparing '+lipid)
                target_csv=pd.read_csv(f"{path}/Interaction_{lipid}/Dataset_{lipid}/Dataset.csv")
                BS_target=target_csv["Binding Site ID"].unique()

                for ref_site in range(0, BS_ref.max()+1):
                    ref_resid_list=ref_csv[ref_csv["Binding Site ID"]==ref_site]["Residue"].tolist()

                    target_BS_count={}
                    for target_site in range(0, BS_target.max()+1):
                        target_resid_list=target_csv[target_csv["Binding Site ID"]==target_site]["Residue"].tolist()
                        res_count=len(set([x for x in ref_resid_list+target_resid_list if x in ref_resid_list and x in target_resid_list]))

                        target_BS_count[int(target_site)]=res_count/len(ref_resid_list+target_resid_list)

                    t=max(target_BS_count.values())
                    # Threshold set to 0.2 based on comparsion to manually assigned matching sites.
                    if t>0.2:
                         target_max_sites=[key for key, value in target_BS_count.items() if value==t]
                    else:
                         target_max_sites=["X"]

                    df=pd.DataFrame({f"{ref_lipid}": ref_site,
                                    f"{lipid}": target_max_sites})

                    tmp_df=tmp_df.append(df, ignore_index=True)
                match_sites=pd.concat([match_sites, tmp_df], axis=1)

            match_sites=match_sites.loc[:,~match_sites.columns.duplicated()]
            print("\nAutomatically predicted matching binding sites for each lipid: \n", match_sites.to_string(index=False))
            print("\nIMPORTANT: check you are happy with the site match by examining lipid poses.\nYou may wish to remove replicate or poorly defined sites.")
            BS_ID_dict=match_sites.to_dict('list')

            # Check for replicate matching sites for each lipid
            for col in match_sites:
                repeat=match_sites[col].replace('X', np.nan).dropna().duplicated().any()
                if repeat:
                    print(f"\nWARNING: A binding site was included twice for {col} - you may wish to check which site matches best.")

        except Exception as e:
            print(e)
            exit()
    elif len(lip_list)==1:
        try:
            ref_lipid=lip_list[0]
            print('\nReference lipid '+ref_lipid)
            ref_csv=pd.read_csv(f"{path}/Interaction_{ref_lipid}/Dataset_{ref_lipid}/Dataset.csv")
            BS_ref=ref_csv["Binding Site ID"].unique().tolist()
            BS_ref.remove(-1)
            BS_ID_dict={ref_lipid: BS_ref}
        except Exception as e:
            print(e)
            exit()
    else:
        print("Must enter more than one lipid to compare sites between lipid species.")
        exit()

    return BS_ID_dict


def get_site_compare(lip_list, BS_ID_dict):
    """
    Generate dictionary of lipid sites to compare in format: dict={"lipid":[list of site IDs]}.
    Users can specify whether to accept the predicted matching sites or modify the prediction.

    Params
    -------
    lip_list: list
        lipids to compare
    BS_ID_dict: dict
        predicted comparable sites for each lipid

    Returns
    -------
    BS_ID_dict: dict
        unmodified or modified sites to compare as keys=lipids, values=comparable binidng sites (in order)
    """
    try:
        if BS_ID_dict:
            print(*[f"\n\t{lip}: {site}" for lip, site in BS_ID_dict.items()])
            res=input("\nAccept predicted matching sites (y/n): ") or 'y'
            if res=='y':
                print("\nCorresponding sites accepted.")
            elif res=='n':
                BS_ID_dict=None
            else:
                print("Must enter y/n")
                exit()
        if not BS_ID_dict:
            print('\n'*2+"Generate dictionary of lipids (keys) and a list of binding site indices (values) to compare.")
            print("Example:\nBindingSite_ID_dict={'POPC': [1, 2, 5, 3],\n\t'POPE': [2, 3, 6, 1],\n\t'CHOL':[1, 3, 5, 'X']}")
            print("\nThe residence time of sites are compared in the listed order \ne.g. POPC site 1 is compared with POPE site 2 and CHOL site 1.")
            print("'X' denotes when an equivilent site does not exist.")
            BS_ID_dict={}
            for idx, lipid in enumerate(lip_list):
                try:
                    print('\n'*2+'-'*10+f'{lipid}'+'-'*10)
                    site_ID_list=input(f"List {lipid} site IDs in comparsion order e.g '0 2 3 X 5': ")
                    remove_char=["[", "]", ","]
                    for rc in remove_char:
                        site_ID_list=site_ID_list.replace(rc,"")
                    site_ID_list=list(site_ID_list.split())
                    site_ID_list=[int(sd) if sd.isnumeric() else 'X' for sd in site_ID_list]
                    BS_ID_dict[lipid]=site_ID_list
                    print(f"\nBindingSite_ID_dict={BS_ID_dict}")
                except Exception as e:
                    print(e)

        print("\nLipid binding sites to compare:")
        print(*[f"\n\t{lip}: {site}" for lip, site in BS_ID_dict.items()])
    except Exception as e:
        print(e)
        exit()
    return BS_ID_dict

def get_BSstat(path, site_dict):
    """
    Load binding site data for lipid sites in site_dict.

    Params
    ------
    path: str
        path
    site_dict: dict
        binding sites to compare

    Returns
    -------
    data: pd.DaraFrame
        kinetic parameters to plot in comparison
    """
    os.makedirs(os.path.join(path, "Lipid_compare"), exist_ok=True)
    data=pd.DataFrame()
    for lipid in site_dict:
        try:
            print('\n'+lipid)
            interactions_csv=pd.read_csv(f"{path}/Interaction_{lipid}/Dataset_{lipid}/Dataset.csv")
            for idx, site in enumerate(site_dict[lipid]):
                if isinstance(site, int):
                    site_csv=interactions_csv[interactions_csv["Binding Site ID"] == site]
                    #site_csv.loc[:,"koff_diff"]=site_csv.loc[:,"Binding Site Koff"] - site_csv.loc[:,"Binding Site Koff Bootstrap avg"]
                    site_csv["koff_diff"]=site_csv["Binding Site Koff"] - site_csv["Binding Site Koff Bootstrap avg"]
                    df=pd.DataFrame({"Lipid": lipid,
                        "Site_Idx": [idx],
                        "BS": site,
                    "koff_diff": site_csv["koff_diff"].unique(),
                    "BS R Squared": site_csv["Binding Site R Squared"].unique(),
                    "BS Residence Time": site_csv["Binding Site Residence Time"].unique(),
                    "BS Residence Time boot error": 1/site_csv["Binding Site Koff Bootstrap avg"].unique() - site_csv["Binding Site Residence Time"].unique(),
                    "BS Occupancy": site_csv["Binding Site Occupancy"].unique()})
                elif site =="X":
                    df=pd.DataFrame({"Lipid": lipid,
                        "Site_Idx": [idx],
                        "BS": site,
                    "koff_diff": 0,
                    "BS R squared": 0,
                    "BS Residence Time": 0,
                    "BS Residence Time boot error": 0,
                    "BS Occupancy": 0})
                else:
                    print(f"Binding Site ID {site} not recognised. Should be an integer or string 'X' for null site")
                    exit()
                data=data.append(df, ignore_index=True)
        except Exception as e:
            print(e)
            pass
    return data


def plot_site_rank(path, site_dict, data):
    """
    Plot lipid residence times for corresponding sites by site index order in the site_dict[lipid] list.

    Params
    ------
    path: str
        path
    site_dict: dict
        binding sites to compare
    data: pd.DataFrame
        kinetic params for lipid sites
    """
    for i in range(0, max(len(v) for v  in site_dict.values())):
        data_idx=data[data["Site_Idx"]==i].reset_index()
        if not (data_idx["BS"] =="X").all():
            fig, ax = plt.subplots(2, 1, figsize=(4.2, 4), sharex=True,  gridspec_kw={'height_ratios': [1, 3]})

            ax[0].scatter(x=data_idx["Lipid"], y=data_idx["BS R Squared"], c='k')
            ax[0].set_ylabel(r"R$^2$")
            ax[0].set_xlabel("")


            #sns.barplot(x="Lipid", y="BS Residence Time", data=data_idx, ax=ax[1], palette="Set2")
            ax[1].bar(x=data_idx["Lipid"],height=data_idx["BS Residence Time"], yerr=[np.zeros(len(data_idx["BS Residence Time boot error"])), data_idx["BS Residence Time boot error"]], color=sns.color_palette("Set2"))
            ax[1].set_ylabel(r"Residence Time ($\mu$s)")
            ax[1].set_xlabel("")

            ax[0].axhline(y=1, ls=":", c='grey')
            ax[0].set_ylim(data_idx["BS R Squared"].min()-0.02, 1.0)
            #ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=45, ha="right")
            plt.xticks(rotation=45, ha="right")

            sns.despine(top=True, right=True)
            plt.tight_layout()
            #plt.show()
            plt.savefig("{}/Lipid_compare/Lipid_compare_BSstats_PyLipID_Site_idx{}_ref_{}_BS{}.pdf".format(path, data_idx.at[0, "Site_Idx"], data_idx.at[0, "Lipid"], data_idx.at[0, "BS"]), format='pdf')
            plt.close()
    print("\nSite comparison complete:", f"{path}/Lipid_compare")
    return

def locate_density_path(path):
    """
    Check the input structure and density file exists and get input sigma factor value.

    Params
    -------
    path: str
        path

    Returns
    -------
    density map location and sigma factor
    """
    print("Enter protein coordinate file and corresponding (aligned) density map.")
    a_name=str(input("Protein coordinate file: "))
    if os.path.isfile(f"{path}/{a_name}"):
        protein_AT_full=f"{path}/{a_name}"
    else:
        print(f"\n{path}/{a_name} not found.")
        a_name=str(input("Define absolute path to protein coordinate file: "))
        if os.path.isfile(f"{a_name}"):
            protein_AT_full=a_name
        else:
            print("\n{a_name} not found.")
            exit()

    d_name=str(input("Density map file name: "))
    if os.path.isfile(f"{path}/{d_name}"):
        density_map=f"{path}/{d_name}"
    else:
        print(f"\n{path}/{d_name} not found.")
        d_name=str(input("Define absolute path to density map file: "))
        if os.path.isfile(f"{d_name}"):
            density_map=d_name
        else:
            print("\n{d_name} not found.")
            exit()

    sigma_factor=float(input("Sigma factor value to display map at (default: 10): ") or 10)
    return protein_AT_full, density_map, sigma_factor

def backmap_poses(path, protocol_path, BS_ID_dict, save_dir):
    """
    For each lipid and each pose in BS_ID_dict, backmap the lipid pose to atomistic detail (using CG2AT before equilibration i.e. lipid pose position unaltered comapred to CG pose). Put comparable poses in same location (within Density_Pose_Compare) where BS ID referens to the refernece lipid.

    Params
    ------
    path: str
        path
    protocol_path: str
        LipIDens path
    BS_ID_dict: dict
        Comparible binding poses for each lipid
    save_dir: str
        save directory

    Returns
    -------
    dens_path: str
        location of backmapped poses
    """
    dens_path=os.path.join(path, "Density_Pose_Compare")
    os.makedirs(dens_path, exist_ok=True)
    print("Location of lipid poses:", dens_path)

    #establish reference lipid path to store comparable poses
    for site in list(BS_ID_dict.values())[0]:
        os.makedirs(dens_path+f"/BS_ID_{site}", exist_ok=True)

    if os.path.isfile(f"{path}/run1/martini_v2.2.itp"):
        ff_type="martini_2-2_charmm36 martini_2-2_charmm36-VS"
    elif os.path.isfile(f"{path}/run1/martini_v2.0.itp"):
        ff_type="martini_2-0_charmm36"
    else:
        ff_type="martini_3-0_charmm36 martini_3-0_charmm36-VS"
    print("Selected forcefield for backmapping: ",ff_type)

    for lipid in BS_ID_dict:
        for idx, site in enumerate(BS_ID_dict[lipid]):
            if site != 'X':
                input_CG_frame=f"{path}/Interaction_{lipid}/Bound_Poses_{lipid}/BSid{site}_rank/BSid{site}_top0.gro"
                save_CG_frame_path=f"{save_dir}/Interaction_{lipid}/Bound_Poses_{lipid}/BSid{site}_rank"
                print(f"\nRunning CG2AT: {lipid} site pose {site}\n")
                try:
                    subprocess.check_call(["{}/rank_sites/run_CG2AT.sh".format(protocol_path), input_CG_frame, save_CG_frame_path, ff_type])
                except subprocess.CalledProcessError as e:
                    print(e.output)
                    pass

                # get direct lipid backmapped pose (before equilibration) and move to directory for comparison
                #if os.path.isfile(f"{save_CG_frame_path}/CG2AT/{lipid}/{lipid}_merged.pdb"):
                if os.path.isfile(f"{save_CG_frame_path}/CG2AT/MERGED/merged_cg2at_de_novo.pdb"):
                    print(f"CG2AT {lipid} pose {site} backmapped")
                    #struc=f"{save_CG_frame_path}/CG2AT/{lipid}/{lipid}_merged.pdb"
                    struc=f"{save_CG_frame_path}/CG2AT/MERGED/merged_cg2at_de_novo.pdb"
                    shutil.copy(struc, dens_path+"/BS_ID_{}".format(list(BS_ID_dict.values())[0][idx]))
                    shutil.move(dens_path+"/BS_ID_{}/merged_cg2at_de_novo.pdb".format(list(BS_ID_dict.values())[0][idx]), dens_path+"/BS_ID_{}/BS_ID_{}_{}_merged.pdb".format(list(BS_ID_dict.values())[0][idx], list(BS_ID_dict.values())[0][idx],lipid))
                else:
                    print(f"CG2AT {lipid} {site} .pdb file not found - check whether backmap was successful")
                    pass

    return dens_path

def pymol_density_compare(path, protein_AT_full, density_map, sigma_factor, dens_path, BS_ID_dict):
    """
    Generate a pymol script comparing top ranked lipid binding poses for each site with densities in proximity to predicted sites. 
    For each site, backmapped top ranked lipid binding poses are loaded and aligned. Residues for each site are shown as spheres scaled
    by the predicted residence times. The density map is segmented in proximity to each predicted lipid binding site (using the given
    sigma factor) for comparison with lipid poses. 
      
    Elements of this function were adapted from PyLipID.

    Params
    -------
    path: str
        path
    protein_AT_full: str
        protein coordinate file
    density_map: str
        density map
    sigma_factor:
        density sigma factor to carve map
    dens_path: str
        path to backmapped atomistic lipid biniding poses for each site
    BS_ID_dict: dict
        Comparible binding poses for each lipid
   
    """
    
    r_lip=list(BS_ID_dict.keys())[0]
    r_num_of_sites=max([i for i in BS_ID_dict[r_lip] if isinstance(i, int)])

    if os.path.isfile(f"{path}/Interaction_{r_lip}/Dataset_{r_lip}/Dataset.csv"):
        ref_lip_csv=f"{path}/Interaction_{r_lip}/Dataset_{r_lip}/Dataset.csv"
    else:
        print("Unable to locate interaction dataset for reference lipid.")
        rlipcsv=str(input("Provide absolute path to Dataset.csv file for {r_lip}: "))
        if os.path.isfile(rlipcsv):
            ref_lip_csv=rlipcsv
        else:
            print("\n{rlipcsv} not found.")
            exit()


    t=f"""
#!/usr/bin/env python3

import os
import numpy as np
import re
from pymol import cmd
from collections import Counter
#cmd.reinitialize()

p_name="Protein"
cmd.set("retain_order",1)
#load and map sites
cmd.load("{protein_AT_full}", p_name)
cmd.load("{density_map}", "Density_map")
cmd.hide("everything")
cmd.show("cartoon", p_name)
cmd.color("white", p_name)
cmd.set("stick_radius", 0.5)
ref_lip="{r_lip}"
ref_num_of_sites={r_num_of_sites}
dens_path="{dens_path}"
sigma_factor={sigma_factor}


# load residence time data for sites and scale spheres
colours= np.array([np.random.choice(np.arange(256, dtype=float), size=3) for dummy in range(ref_num_of_sites+1)])
print("Colours", colours, len(colours))



with open("{ref_lip_csv}", "r") as f:
    data_lines = f.readlines()

column_names = data_lines[0].strip().split(",")
for column_idx, column_name in enumerate(column_names):
    if column_name == "Residue":
        column_id_residue_list = column_idx
    elif column_name == "Residue ID":
        column_id_residue_index = column_idx
    elif column_name == "Binding Site ID":
        column_id_BS = column_idx
    elif column_name == "Residence Time":
        column_id_value_to_show = column_idx

residue_list = []
residue_rank_set = []
binding_site_identifiers = []
values_to_show = []
for line in data_lines[1:]:
    data_list = line.strip().split(",")
    residue_list.append(data_list[column_id_residue_list])
    residue_rank_set.append(data_list[column_id_residue_index])
    binding_site_identifiers.append(float(data_list[column_id_BS]))
    values_to_show.append(data_list[column_id_value_to_show])

with open("{protein_AT_full}", "r") as f:
    pdb_lines = f.readlines()
residue_identifiers = []
for line in pdb_lines:
    line_stripped = line.strip()
    if line_stripped[:4] == "ATOM":
        identifier = (line_stripped[22:26].strip(), line_stripped[17:20].strip(), line_stripped[21].strip())
#                           residue index,              resname,                     chain id
        if len(residue_identifiers) == 0:
            residue_identifiers.append(identifier)
        elif identifier != residue_identifiers[-1]:
            residue_identifiers.append(identifier)

# set sphere scales 
values_to_show = np.array(values_to_show, dtype=float)
MIN = np.percentile(np.unique(values_to_show), 5)
MAX = np.percentile(np.unique(values_to_show), 100)
X = (values_to_show - np.percentile(np.unique(values_to_show), 50))/(MAX - MIN)
SCALES = 1/(0.5 + np.exp(-X * 5))

    """

    t+=r"""

residue_list = np.array(residue_list, dtype=str)
residue_rank_set = np.array(residue_rank_set, dtype=int)
binding_site_identifiers = np.array(binding_site_identifiers, dtype=int)
residue_identifiers = list(residue_identifiers)

for bs_id in np.arange(ref_num_of_sites+1):
    cmd.set_color(f"tmp_{bs_id}", list(colours[bs_id]))
    res_sel_list=[]
    id_pdb_unq=[]

    for entry_id in np.where(binding_site_identifiers == bs_id)[0]:
        selected_residue = residue_list[entry_id]
        selected_residue_rank = residue_rank_set[entry_id]
        identifier_from_pdb = residue_identifiers[selected_residue_rank]
        id_pdb_unq.append(identifier_from_pdb[2])
        if re.findall("[a-zA-Z]+$", selected_residue)[0] != identifier_from_pdb[1]:
            raise IndexError(
            "The {}th residue in the provided pdb file ({}{}) is different from that in the simulations ({})!".format(
                                                                                            entry_id+1,
                                                                                            identifier_from_pdb[0],
                                                                                            identifier_from_pdb[1],
                                                                                            selected_residue)
                                                                                            )
        if identifier_from_pdb[2] != " ":
            cmd.select(f"BSid{bs_id}_{selected_residue}",
            f"{p_name} and chain {identifier_from_pdb[2]} and resid {identifier_from_pdb[0]} and (not name C+O+N)")
            res_sel_list.append(identifier_from_pdb[0])
        else:
            cmd.select(f"BSid{bs_id}_{selected_residue}",
            f"{p_name} and resid {identifier_from_pdb[0]} and (not name C+O+N)")
            res_sel_list.append(identifier_from_pdb[0])
        cmd.show("spheres", f"BSid{bs_id}_{selected_residue}")
        cmd.set("sphere_scale", SCALES[entry_id], selection=f"BSid{bs_id}_{selected_residue}")  
        cmd.color(f"tmp_{bs_id}", f"BSid{bs_id}_{selected_residue}")
    cmd.group(f"BSid{bs_id}", f"BSid{bs_id}_*")
    
    # binding site residues for density selection
    res_sel_list="+".join(res_sel_list)

    # load and align top ranked lipid binding poses
    if os.path.isdir(f"{dens_path}/BS_ID_{bs_id}"):
        fle_lst=os.listdir(f"{dens_path}/BS_ID_{bs_id}")
        for fle in fle_lst:
            print(fle[:-4])
            cmd.load(f"{dens_path}/BS_ID_{bs_id}/{fle}")
            if len(Counter(id_pdb_unq))==1:
                chain_unq=Counter(id_pdb_unq).most_common(1)[0][0]
                cmd.cealign(target=f"{p_name} and chain {chain_unq}", mobile=fle[:-4])
            else:
                cmd.cealign(target=p_name, mobile=fle[:-4])
    
        # generate density around site
        #cmd.isomesh(f"BS_ID_{bs_id}_map", "Density_map", level=sigma_factor, selection=f"{p_name} and resid {res_sel_list} around 5")
        cmd.isomesh(f"BS_ID_{bs_id}_map", "Density_map", level=sigma_factor, selection=f"{p_name} and BSid{bs_id}* and sidechain", carve=6)
        cmd.group(f"BS_ID_{bs_id}", f"BS_ID_{bs_id}*")
        cmd.hide("cartoon", f"BS_ID_{bs_id}")
        cmd.color(f"tmp_{bs_id}", f"BS_ID_{bs_id}")       
cmd.center(p_name)

# save session
cmd.save(f"{dens_path}/Lipid_poses_density_compare.pse")
    """
    with open(f"{dens_path}/Lipid_poses_density_compare.py", "w+") as f:
        f.write(t)

    return
