# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 16:18:42 2023

@author: dboateng
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pyClimat.plot_utils import apply_style


data_path = "C:/Users/dboateng/Desktop/Datasets/Alps_d18op/lapse_rate_summary.csv"
data_west_path = "C:/Users/dboateng/Desktop/Datasets/Alps_d18op/lapse_rate_summary_west.csv"
data_north_path = "C:/Users/dboateng/Desktop/Datasets/Alps_d18op/lapse_rate_summary_north.csv"

global_river = "C:/Users/dboateng/Desktop/Datasets/Alps_d18op/global_rivers.csv"

campani_path = "C:/Users/dboateng/Desktop/Datasets/Alps_d18op/campani_estimates.csv"

path_to_store = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

data = pd.read_csv(data_path)
data_west = pd.read_csv(data_west_path)
data_north = pd.read_csv(data_north_path)
data_campani = pd.read_csv(campani_path)
data_global = pd.read_csv(global_river)


# plot 
apply_style(fontsize=28, style="seaborn-paper", linewidth=3,)

fig, ax = plt.subplots(figsize=(20, 18))

boxpot = sns.boxplot(data=data, saturation=0.8, ax=ax, linewidth=3)
stripplot = sns.stripplot(data=data_west, size=23, ax=ax, edgecolor="black", linewidth=3)
stripplot = sns.stripplot(data=data_north, size=23, ax=ax, edgecolor="white", linewidth=3)
sns.pointplot(data=data_campani, x="SMB", y="elevation", ax=ax,
                      markers="x", linestyles="", scale=4,
                      dodge=True, color="black")

sns.pointplot(data=data_global, x="SMB", y="elevation", ax=ax,
                      markers="+", linestyles="", scale=4,
                      dodge=True, color="black")

ax.set_ylabel("Elevation [m]", fontweight="bold", fontsize=28)
ax.grid(True, linestyle="--", color="grey")

plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
plt.savefig(os.path.join(path_to_store, "summary_west_north_diff.svg"), bbox_inches="tight", format= "svg")
