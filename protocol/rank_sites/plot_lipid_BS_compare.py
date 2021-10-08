#!/usr/bin/env python3

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
"""

Script to compare lipid binding parames at the same site in PyLipID.

Author: TBG Ansell

"""
### User input ###

path="XXX/Interaction_Lipids_all_0.475_0.7"
#lipid: [Binding sites] - order of binding sites correpsonds to 1) 2) 3)
site_dict={"POPE":[5, 1, 6],#
        "POPG":[5, 0, 6],#
        "CDL2":[6, 1, 9] #
        }

lipid_names = {"POPE": "POPE",
        "POPG": "POPG",
        "CDL2": "CDL"
        }


lipid_labels = {
        "POPE": '#CA00FE',
        "POPG": 'lightgrey',
        "CDL2": 'mediumturquoise'
        }

BS_cats={"BS R Squared": r"R$^2$",
        "BS Residence Time":r"Residence Time ($\mu$s)" ,
 # "BS Occupancy": "Occupancy (%)",
 # "koff_diff": "dKoff"
  }

### Code below ###
def get_BSstat():
    data=pd.DataFrame()
    for lipid in site_dict:
        print(lipid)
        interactions_csv=pd.read_csv("{}/Interaction_{}/Dataset_{}/Dataset.csv".format(path,lipid,lipid))
        for idx, site in enumerate(site_dict[lipid]):
            if site != "X":
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
            data=data.append(df, ignore_index=True)

    return data

def plot_data(data):
    for i in range(0, len(site_dict["POPE"])):
        data_idx=data[data["Site_Idx"]==i]

        fig, ax = plt.subplots(2, 1, figsize=(4.2, 4), sharex=True,  gridspec_kw={'height_ratios': [1, 3]})
        for idx, cat in enumerate(BS_cats):
            if cat=="BS R Squared":
                ax[idx].scatter(x=data_idx["Lipid"], y=data_idx["{}".format(cat)], c='k')
                ax[idx].set_ylabel(BS_cats[cat])
                ax[idx].set_xlabel("")
            else:
                #sns.stripplot(x="Lipid", y=cat, data=data, ax=ax[idx], palette=lipid_labels, s=7)
                sns.barplot(x="Lipid", y=cat, data=data_idx,ax=ax[idx], palette=lipid_labels)
                #sns.barplot(x="Lipid", y=cat, data=data_idx, hue="Site_Idx",ax=ax[idx])
                ax[idx].set_ylabel(BS_cats[cat])
                ax[idx].set_xlabel("")

        ax[0].axhline(y=1, ls=":", c='grey')
        #ax[2].axhline(y=100,  ls=":", c='grey')
        ax[0].set_ylim(0.7,1.0)
        #ax[1].set_ylim(0, 15)
        ax[1].set_xticklabels(lipid_names.values())

        sns.despine(top=True, right=True)
        plt.tight_layout()
        #plt.show()
        plt.savefig("{}/Lipid_compare_BSstats_PyLipID_Site_idx_{}.pdf".format(path, i), format='pdf')
        plt.close()
    return




data=get_BSstat()
print(data)
plot_data(data)
