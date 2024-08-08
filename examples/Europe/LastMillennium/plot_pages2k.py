# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 15:15:34 2024

@author: dboateng
"""

import os
from io import StringIO
import sys
import numpy as np
from scipy.stats import zscore
import pandas as pd
import lipd
import pyleoclim as pyleo
import tempfile
import pooch

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates
import cartopy.feature as cfeature
import seaborn as sns

from pyClimat.plot_utils import *
from help_functions import *

path_to_lipd_iso2k = "C:/Users/dboateng/Desktop/Datasets/iso2k/iso2k1_0_1"
path_to_pickle_iso2k = "C:/Users/dboateng/Desktop/Datasets/iso2k/iso2k1_0_1.pkl"
path_to_pickle_pages2k = "C:/Users/dboateng/Desktop/Datasets/iso2k/pages2k.pkl"
path_to_lipd_pages2k = "C:/Users/dboateng/Desktop/Datasets/iso2k/pages2k"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"


# with SupressOutputs():
#     df = lipd2df(lipd_dirpath=path_to_lipd_pages2k, pkl_filepath=path_to_pickle_pages2k,
#                  crit="paleoData_useInGlobalTemperatureAnalysis")
    
    
pages2k_data = pd.read_pickle(path_to_pickle_pages2k)

# list of markers and colors for the different archive_type
markers = ["p", "p", "o", "v", "d", "*", "s", "s", "8", "D", "^"]
colors = [
    np.array([1.0, 0.83984375, 0.0]),
    np.array([0.73828125, 0.71484375, 0.41796875]),
    np.array([1.0, 0.546875, 0.0]),
    np.array([0.41015625, 0.41015625, 0.41015625]),
    np.array([0.52734375, 0.8046875, 0.97916667]),
    np.array([0.0, 0.74609375, 1.0]),
    np.array([0.25390625, 0.41015625, 0.87890625]),
    np.array([0.54296875, 0.26953125, 0.07421875]),
    np.array([1, 0, 0]),
    np.array([1.0, 0.078125, 0.57421875]),
    np.array([0.1953125, 0.80078125, 0.1953125]),
]

# create the plot

apply_style(fontsize=25, style=None, linewidth=2) 

projection = ccrs.Robinson(central_longitude=0, globe=None)
#projection = ccrs.PlateCarree()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 14), subplot_kw={"projection":  
                                                                        projection})

# add plot title
ax.set_title(
    f"PAGES2k Network (n={len(pages2k_data)})")

# set the base map
# ----------------
ax.set_global()

# add coast lines
ax.coastlines(resolution = "50m", linewidth=1.5, color="grey")

# add land fratures using gray color
ax.add_feature(cfeature.LAND, facecolor="gray", alpha=0.3)
ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth = 0.3)

# add gridlines for latitude and longitude
ax.gridlines(edgecolor="gray", linestyle=":")


# plot the different archive types
# -------------------------------

# extract the name of the different archive types
archive_types = pages2k_data.archiveType.unique()

# plot the archive_type using a forloop
for i, type_i in enumerate(archive_types):
    df = pages2k_data[pages2k_data["archiveType"] == type_i]
    # count the number of appearances of the same archive_type
    count = df["archiveType"].count()
    # generate the plot
    ax.scatter(
        df["geo_meanLon"],
        df["geo_meanLat"],
        marker=markers[i],
        color=colors[i],
        edgecolor="k",
        s=120,
        transform=ccrs.Geodetic(),
        label=f"{type_i} (n = {count})",
    )
# add legend to the plot
ax.legend(
    scatterpoints=1,
    bbox_to_anchor=(0, -0.6),
    loc="lower left",
    ncol=2,
    #fontsize=15,
)

plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.06)

plt.savefig(os.path.join(path_to_plots, "pages2k_overview.png"), format= "png", bbox_inches="tight")
