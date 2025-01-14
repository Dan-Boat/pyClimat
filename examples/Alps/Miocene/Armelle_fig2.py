# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 17:42:36 2024

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
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean
from pyClimat.variables import extract_var


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/armelle"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"



proxy_path = "C:/Users/dboateng/Desktop/Datasets/Miocene_proxydata_new.csv"

df = pd.read_csv(proxy_path, usecols=["Latitude", "Longitude", "Age_Ma", "d18Ow_VSMOW", "symbol"])

# select age range
df_mco = df[(df["Age_Ma"] >= 14.7) & (df["Age_Ma"] <= 16.9)]

df_mco = df_mco.groupby(["Longitude", "Latitude"]).agg(d18Op=("d18Ow_VSMOW", "mean")).reset_index()

df_mmct = df[(df["Age_Ma"] >= 13.8) & (df["Age_Ma"] <= 14.7)]

df_mmct = df_mmct.groupby(["Longitude", "Latitude"]).agg(d18Op=("d18Ow_VSMOW", "mean")).reset_index()




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
    
    return d18op_alt
    
# read data 
Mio278_d18Op = extract_vars_and_analysis(data=W1E1_278_data, wiso=W1E1_278_wiso)
Mio450_d18Op = extract_vars_and_analysis(data=W1E1_450_data, wiso=W1E1_450_wiso)


# plot comments (reduce the resolution of the coastlines)
# plot for isotopes


apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)

fig,(ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20,12), subplot_kw={"projection":projection})

plot_annual_mean(variable="$\delta^{18}$Op VSMOW", data_alt=Mio278_d18Op, ax=ax1,
                  cmap="Spectral_r", units="‰", vmax=0, vmin=-15, 
                levels=22, level_ticks=9, add_colorbar=True, cbar_pos= [0.35, 0.05, 0.35, 0.02], 
                orientation="horizontal", plot_coastlines=True, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=True, 
                plot_projection=projection, title="(a) MIO 278 ppm", center=False, domain="Alps")

ax1.scatter(x=df_mmct["Longitude"], y=df_mmct["Latitude"], c=df_mmct["d18Op"], cmap="Spectral_r",
            vmax=0, vmin=-15, edgecolor="black", s= 300, transform=ccrs.PlateCarree(),
            linewidth=2)

plot_annual_mean(variable="$\delta^{18}$Op VSMOW", data_alt=Mio450_d18Op, ax=ax2,
                  cmap="Spectral_r", units="‰", vmax=0, vmin=-15, 
                levels=22, level_ticks=9, add_colorbar=False, 
                orientation="horizontal", plot_coastlines=True, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=True, 
                plot_projection=projection, title="(b) MIO 450 ppm", center=False,
                domain="Alps")

ax2.scatter(x=df_mco["Longitude"], y=df_mco["Latitude"], c=df_mco["d18Op"], cmap="Spectral_r",
            vmax=0, vmin=-15, edgecolor="black", s= 300, transform=ccrs.PlateCarree(),
            linewidth=2)
            
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, wspace=0.1)

plt.savefig(os.path.join(path_to_plots, "d18Ow_model_proxy_mean.pdf"), format= "pdf", bbox_inches="tight", dpi=600)