#!/usr/bin/env python3

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

"""

Script to compare PyLipID outputs across binding sites and rank them.

Author: TBG Ansell

"""

def get_data(path, lipid):
    """
    Generate summary of binding site data as csv.

    Params:
    -------
    path: str
        path
    lipid: str
        lipid
    Returns:
    --------
    data: pandas.DataFrame()
        summary of lipid interaction data per binding site
    """
    ### Load the basic data ###
    # All data #
    try:
        print("\nObtaining PyLipID outputs for screening:", lipid, "\n")
        interactions_csv=pd.read_csv("{}/Interaction_{}/Dataset_{}/Dataset.csv".format(path, lipid, lipid))
        interactions_csv=interactions_csv.filter(regex='Bind', axis=1).drop_duplicates()
        interactions_csv=interactions_csv[interactions_csv["Binding Site ID"] >=0]
        # Store and load summary of binding site data #
        interactions_csv.to_csv("{}/Interaction_{}/BS_summary_{}.csv".format(path, lipid, lipid), index=False)
        data=pd.read_csv("{}/Interaction_{}/BS_summary_{}.csv".format(path,lipid,lipid), delimiter=",")
        print("Binding site data:\n", "{}/Interaction_{}/BS_summary_{}.csv".format(path,lipid,lipid))
    except Exception as e:
        print("PyLipID Data (Dataset.csv) not found for {} - check PyLipID run was successful".format(lipid))
        print(e)
        data=None
        pass
    return data

### Plot the data ###
def plot_screen_data(data, path, lipid):
    """
    Generate plot of screened site data in site order.

    Params:
    -------
    data: pandas.DataFrame()
        summary of lipid interaction data per binding site
    path: str
        path
    lipid: str
        lipid
    """
    # BS categories to plot #
    BS_cats={"koff_diff": r"$\Delta$$k_{off}$",
      "Binding Site Residence Time":r"Residence Time ($\mu$s)" ,
      "Binding Site Occupancy": "Occupancy (%)",
      "Binding Site Surface Area": r"Surface Area ($nm^2$)"}
    
    num_IDs=data["Binding Site ID"].max()
    if num_IDs < 15:
        x_scale=1
    else:
        x_scale=num_IDs/15

    fig, ax = plt.subplots(4, 1, figsize=(4*x_scale,6))
    pal=['#577590', '#90BE6D', '#F8961E', '#F94144']

    data["koff_diff"]=data["Binding Site Koff"] - data["Binding Site Koff Bootstrap avg"]

    for idx, cat in enumerate(BS_cats):
        if cat=="koff_diff":
            sns.stripplot(x="Binding Site ID", y=cat, data=data, ax=ax[idx], order=data.sort_values(by=cat, ascending=False, key=lambda k: abs(k))["Binding Site ID"].tolist(), color=pal[idx], s=7)
        else:
            sns.stripplot(x="Binding Site ID", y=cat, data=data, ax=ax[idx], order=data.sort_values(by=cat)["Binding Site ID"].tolist(), color=pal[idx], s=7)
        ax[idx].set_ylabel(BS_cats[cat])
        ax[idx].set_xlabel("")
        #plt.close()

    ax[3].set_xlabel("Binding Site")
    ax[0].axhline(y=0, ls=":", c='grey')
    ax[2].axhline(y=100,  ls=":", c='grey')

    sns.despine(top=True, right=True)
    plt.tight_layout()
    #plt.show()
    plt.savefig("{}/Interaction_{}/Site_stats_rank_compare.pdf".format(path, lipid), format='pdf')
    plt.close()
    print("Ranked comparison of binding sites complete:", "{}/Interaction_{}/Site_stats_rank_compare.pdf".format(path, lipid), "\n")
    return
