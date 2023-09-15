# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 10:29:41 2023

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl 
import cartopy.crs as ccrs
from pyClimat.plot_utils import apply_style
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#path to data 

path_to_low_elev_data = "C:/Users/dboateng/Desktop/Datasets/Alps_d18op/low_elevation_distribution.csv"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"
# read data 

df = pd.read_csv(path_to_low_elev_data)



apply_style(fontsize=28, style="seaborn-paper", linewidth=3,)

fig, ax = plt.subplots(figsize=(20, 18))

boxpot = sns.boxplot(data=df, saturation=0.8, ax=ax, linewidth=3, y="d18op", x="source")
stripplot = sns.stripplot(data=df, size=23, ax=ax, edgecolor="black", linewidth=3, y="d18op", x="source",
                          hue="age", palette=sns.mpl_palette("winter_r", as_cmap=True))

ax.set_ylim(-22.5, 0)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', which='major', length=8, width=3)
ax.tick_params(axis='y', which='minor', length=4, width=1.5)

plt.savefig(os.path.join(path_to_plots, "low_elevation_distribution_age.svg"), format= "svg", bbox_inches="tight", dpi=600)