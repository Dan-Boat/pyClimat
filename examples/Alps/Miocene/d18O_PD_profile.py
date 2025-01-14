# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:04:35 2025

@author: dboateng
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 14:26:16 2024

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed, read_from_path
from pyClimat.analysis import compute_lterm_mean, extract_profile
from pyClimat.variables import extract_var


main_path = "D:/Datasets/Model_output_pst/PD"
gnip_path = "D:/Datasets/GNIP_data/world/scratch/station_world_overview_5years.csv" 


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/armelle"




df_gnip = pd.read_csv(gnip_path)
df_gnip = df_gnip[(df_gnip["lat"] >= 43) & (df_gnip["lat"] <= 49.5)]
df_gnip = df_gnip[(df_gnip["lon"] >= 3) & (df_gnip["lon"] <= 17)]


#load datasets 
PD_data = read_from_path(main_path, "PD_1980_2014_monthly.nc", decode=True)
PD_wiso = read_from_path(main_path, "PD_1980_2014_monthly_wiso.nc", decode=True)


# # extract values and compute long-term means

def extract_vars_and_analysis(data, wiso):
    
    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)

    d18op_alt = compute_lterm_mean(data=d18op, time="annual")
    
    prof_mean = extract_profile(data=d18op_alt, maxlon=25, minlon=-5, maxlat=47, minlat=45, dim="lon",
                                method="mean")
    
    prof_std = extract_profile(data=d18op_alt, maxlon=25, minlon=-5, maxlat=47, minlat=45, dim="lon",
                                method="std")
    
    return prof_mean, prof_std
    

PD_d18Op_mean, PD_d18Op_std = extract_vars_and_analysis(data=PD_data, wiso=PD_wiso)


apply_style2(fontsize=25, style="seaborn-talk", linewidth=3, usetex=True)
    
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 16,))


ax.fill_between(PD_d18Op_mean.index, y1=PD_d18Op_mean - PD_d18Op_std, 
                y2=PD_d18Op_mean + PD_d18Op_std, color="black", alpha=0.1)

ax.plot(PD_d18Op_mean.index, PD_d18Op_mean, linewidth=3, color="black", label="ECHAM5-wiso (1979-2014)")

ax.set_ylabel("$\delta^{18}$Op VSMOW (‰)", fontsize=24, fontweight="bold")   
ax.set_xlabel("Longitude (°E)", fontsize=24, fontweight="bold")
    
ax.xaxis.tick_top()   
ax.xaxis.set_label_position('top')


ax.scatter(y=df_gnip["d18op"], x=df_gnip["lon"], s=350, marker="o", color="black", 
                alpha=0.8, label="GNIP stations (>5 years of data)")

ax.legend(bbox_to_anchor=(0.01, -0.15, 1., 0.102), loc=3, ncol=3, borderaxespad=0., frameon = True, 
                  fontsize=20)


plt.savefig(os.path.join(path_to_plots, "d18O_PD_profile.pdf"), format= "pdf", bbox_inches="tight", dpi=600)
