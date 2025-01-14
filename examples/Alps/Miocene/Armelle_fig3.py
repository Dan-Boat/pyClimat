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

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"


df_gnip = pd.read_csv(gnip_path)
df_gnip = df_gnip[(df_gnip["lat"] >= 35) & (df_gnip["lat"] <= 50)]
df_gnip = df_gnip[(df_gnip["lon"] >= 0) & (df_gnip["lon"] <= 12)]


#load datasets 
PD_data = read_from_path(main_path, "PD_1980_2014_monthly.nc", decode=True)
PD_wiso = read_from_path(main_path, "PD_1980_2014_monthly_wiso.nc", decode=True)


proxy_path = "C:/Users/dboateng/Desktop/Datasets/Miocene_proxydata_neww.csv"

df = pd.read_csv(proxy_path, usecols=["Latitude", "Longitude", "Age_Ma", "d18Ow_VSMOW", "d18Ow_err_VSMOW"])

# select age range
df_mco = df[(df["Age_Ma"] >= 14.7) & (df["Age_Ma"] <= 16.9)]

df_mco = df_mco.groupby(["Longitude", "Latitude"]).agg(d18Op=("d18Ow_VSMOW", "mean"), d18Op_err=("d18Ow_err_VSMOW", "mean")).reset_index()

df_mmct = df[(df["Age_Ma"] >= 13.8) & (df["Age_Ma"] <= 14.7)]

df_mmct = df_mmct.groupby(["Longitude", "Latitude"]).agg(d18Op=("d18Ow_VSMOW", "mean"), d18Op_err=("d18Ow_err_VSMOW", "mean")).reset_index()




mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"
#W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

# reading data 
# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)




# # extract values and compute long-term means

def extract_vars_and_analysis(data, wiso):
    
    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)

    d18op_alt = compute_lterm_mean(data=d18op, time="annual")
    
    prof_mean = extract_profile(data=d18op_alt, maxlon=10, minlon=6, maxlat=49, minlat=35, dim="lat",
                                method="mean")
    
    prof_std = extract_profile(data=d18op_alt, maxlon=10, minlon=6, maxlat=49, minlat=35, dim="lat",
                                method="std")
    
    return prof_mean, prof_std
    
# read data 
Mio278_d18Op_mean, Mio278_d18Op_std = extract_vars_and_analysis(data=W1E1_278_data, wiso=W1E1_278_wiso)
Mio450_d18Op_mean, Mio450_d18Op_std = extract_vars_and_analysis(data=W1E1_450_data, wiso=W1E1_450_wiso)

PD_d18Op_mean, PD_d18Op_std = extract_vars_and_analysis(data=PD_data, wiso=PD_wiso)


apply_style2(fontsize=25, style="seaborn-talk", linewidth=3, usetex=True)
    
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 16,))

ax.fill_betweenx(Mio278_d18Op_mean.index, x1=Mio278_d18Op_mean - Mio278_d18Op_std, 
                x2=Mio278_d18Op_mean + Mio278_d18Op_std, color="blue", alpha=0.1)
ax.plot(Mio278_d18Op_mean, Mio278_d18Op_mean.index, linewidth=3, color="blue", label="MIO 278 ppm")

ax.fill_betweenx(Mio450_d18Op_mean.index, x1=Mio450_d18Op_mean - Mio450_d18Op_std, 
                x2=Mio450_d18Op_mean + Mio450_d18Op_std, color="red", alpha=0.1)
ax.plot(Mio450_d18Op_mean, Mio450_d18Op_mean.index, linewidth=3, color="red", label="MIO 450 ppm")

ax.fill_betweenx(PD_d18Op_mean.index, x1=PD_d18Op_mean - PD_d18Op_std, 
                x2=PD_d18Op_mean + PD_d18Op_std, color="black", alpha=0.1)
ax.plot(PD_d18Op_mean, PD_d18Op_mean.index, linewidth=3, color="black", label="PD")

ax.set_xlabel("$\delta^{18}$Op VSMOW (‰)", fontsize=24, fontweight="bold")   
ax.set_ylabel("Latitude (°N)", fontsize=24, fontweight="bold")
    
ax.xaxis.tick_top()   
ax.xaxis.set_label_position('top')

ax.scatter(x=df_mmct["d18Op"], y=df_mmct["Latitude"], s=450, marker="*", color="blue", 
               alpha=0.8, label="MMCT")

# ax.errorbar(x=df_mmct["d18Op"], y=df_mmct["Latitude"], xerr=df_mmct["d18Op_err"],
#                  yerr=None, ls="none", color="black", capthick=2, capsize=3)


ax.scatter(x=df_mco["d18Op"], y=df_mco["Latitude"], s=350, marker="s", color="red", 
               alpha=0.8, label="MCO")
# ax.errorbar(x=df_mco["d18Op"], y=df_mco["Latitude"], xerr=df_mco["d18Op_err"],
#                  yerr=None, ls="none", color="black", capthick=2, capsize=3)



ax.scatter(x=df_gnip["d18op"], y=df_gnip["lat"], s=350, marker="o", color="black", 
                alpha=0.8, label="PD")

ax.legend(bbox_to_anchor=(0.01, -0.15, 1., 0.102), loc=3, ncol=3, borderaxespad=0., frameon = True, 
                  fontsize=20)


plt.savefig(os.path.join(path_to_plots, "d18Ow_model_proxy_profile.pdf"), format= "pdf", bbox_inches="tight", dpi=600)
