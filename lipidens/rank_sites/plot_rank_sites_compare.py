#!/usr/bin/env python3
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import operator
import numpy as np

"""

Script to rank the residence times of different lipids binding to the same site using outputs from PyLipID.

Author: TBG Ansell

"""

def compare_sites(path, lip_list):
    if len(lip_list)>1:
        try:
            match_sites=pd.DataFrame()

            ref_lipid=lip_list[0]
            print('Analysing '+ref_lipid)
            ref_csv=pd.read_csv(f"{path}/Interaction_{ref_lipid}/Dataset_{ref_lipid}/Dataset.csv")
            BS_ref=ref_csv["Binding Site ID"].unique()

            for lipid in lip_list[1:]:
                tmp_df=pd.DataFrame()
                print('Analysing '+lipid)
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
                    if t>0.2:
                         target_max_sites=[key for key, value in target_BS_count.items() if value==t]
                    else:
                         target_max_sites=["X"]
                    #target_max_sites=[value for key, value in target_BS_count.items() if value==t]

                    df=pd.DataFrame({f"{ref_lipid}": ref_site,
                                    f"{lipid}": target_max_sites})

                    tmp_df=tmp_df.append(df, ignore_index=True)

                match_sites=pd.concat([match_sites, tmp_df], axis=1)

            match_sites=match_sites.loc[:,~match_sites.columns.duplicated()]
            print(match_sites)
            BS_ID_dict=match_sites.to_dict('list')


            print(BS_ID_dict)
        except Exception as e:
            print(e)
            exit()
    else:
        print("Must enter more than one lipid to compare sites between lipid species.")
        exit()

    return BS_ID_dict


def get_site_compare(lip_list):
    """
    Generate dictionary of lipid sites to compare in format: dict={"lipid":[list of site IDs]}.
    """
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
    return BS_ID_dict

def get_BSstat(path, site_dict):
    """
    Load binding site data for lipid sites in site_dict.
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
                    "BS Occupancy": site_csv["Binding Site Occupancy"].unique()})
                elif site =="X":
                    df=pd.DataFrame({"Lipid": lipid,
                        "Site_Idx": [idx],
                        "BS": site,
                    "koff_diff": 0,
                    "BS R squared": 0,
                    "BS Residence Time": 0,
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
    """
    for i in range(0, max(len(v) for v  in site_dict.values())):
        data_idx=data[data["Site_Idx"]==i]

        fig, ax = plt.subplots(2, 1, figsize=(4.2, 4), sharex=True,  gridspec_kw={'height_ratios': [1, 3]})

        ax[0].scatter(x=data_idx["Lipid"], y=data_idx["BS R Squared"], c='k')
        ax[0].set_ylabel(r"R$^2$")
        ax[0].set_xlabel("")


        sns.barplot(x="Lipid", y="BS Residence Time", data=data_idx, ax=ax[1], palette="Set2")
        ax[1].set_ylabel(r"Residence Time ($\mu$s)")
        ax[1].set_xlabel("")

        ax[0].axhline(y=1, ls=":", c='grey')
        ax[0].set_ylim(data_idx["BS R Squared"].min()-0.02, 1.0)
        ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=45, ha="right")

        sns.despine(top=True, right=True)
        plt.tight_layout()
        #plt.show()
        plt.savefig(f"{path}/Lipid_compare/Lipid_compare_BSstats_PyLipID_Site_idx_{i}.pdf", format='pdf')
        plt.close()
    print("\nSite comparison complete:", "{}/Lipid_compare".format(path))
    return
