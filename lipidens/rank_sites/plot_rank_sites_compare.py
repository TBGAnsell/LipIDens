#!/usr/bin/env python3
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

"""

Script to rank the residence times of different lipids binding to the same site using outputs from PyLipID.

Author: TBG Ansell

"""

def get_BSstat(path, site_dict):
    """
    Load binding site data for lipid sites in site_dict.
    """
    os.makedirs(os.path.join(path, "Lipid_compare"), exist_ok=True)
    data=pd.DataFrame()
    for lipid in site_dict:
        try:
            print(lipid)
            interactions_csv=pd.read_csv(f"{path}/Interaction_{lipid}/Dataset_{lipid}/Dataset.csv")
            for idx, site in enumerate(site_dict[lipid]):
                print(site)
                if isinstance(site, int):
                    site_csv=interactions_csv[interactions_csv["Binding Site ID"] == site]
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
