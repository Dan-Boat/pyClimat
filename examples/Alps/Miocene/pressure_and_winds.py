# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 11:10:26 2023

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_from_path
from pyClimat.analysis import extract_var, compute_lterm_mean, compute_lterm_diff


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_path = "D:/Datasets/Model_output_pst/a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h/output_processed/"
W1E1_450_path = "D:/Datasets/Model_output_pst/a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h/output_processed/"
W1E1_PI_path = "D:/Datasets/Model_output_pst/a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h/output_processed/"

# reading data 
# read data (monthly)

# main data
filename_lterm = "1003_1017_monthly.nc"

W1E1_278_data = read_from_path(os.path.join(W1E1_278_path, "MONTHLY_MEANS"), filename_lterm, decode=True)
W1E1_450_data = read_from_path(os.path.join(W1E1_450_path, "MONTHLY_MEANS"), filename_lterm, decode=True)
W1E1_PI_data = read_from_path(os.path.join(W1E1_PI_path, "MONTHLY_MEANS"), filename_lterm, decode=True)


# extract values and compute long-term means

def extract_vars_and_analysis(data, pi_data):
    
    slp = extract_var(Dataset=data , varname="slp", units="hPa")
    u = extract_var(Dataset= data , varname="u", lev=850, lev_units="hPa")
    v = extract_var(Dataset= data , varname="v", lev=850, lev_units="hPa")
    
    
    slp_pi = extract_var(Dataset=pi_data , varname="slp", units="hPa")
    u_pi = extract_var(Dataset= pi_data , varname="u", lev_units="hPa", lev=850)
    v_pi = extract_var(Dataset= pi_data , varname="v", lev_units="hPa", lev=850)
   
    
    #compute climatologies difference
    
    slp_diff_djf = compute_lterm_diff(data_control=slp_pi, data_main=slp, time="season",
                                      season="DJF")
    u_diff_djf = compute_lterm_diff(data_control=u_pi, data_main=u, time="season",
                                      season="DJF")
    v_diff_djf = compute_lterm_diff(data_control=v_pi, data_main=v, time="season",
                                      season="DJF")
    
    slp_diff_jja = compute_lterm_diff(data_control=slp_pi, data_main=slp, time="season",
                                      season="JJA")
    u_diff_jja = compute_lterm_diff(data_control=u_pi, data_main=u, time="season",
                                      season="JJA")
    v_diff_jja = compute_lterm_diff(data_control=v_pi, data_main=v, time="season",
                                      season="JJA")
   
    
   
    
    
    return_data = {"slp_DJF":slp_diff_djf, "u_DJF":u_diff_djf, "v_DJF":v_diff_djf,
                   "slp_JJA":slp_diff_jja, "u_JJA":u_diff_jja, "v_JJA":v_diff_jja}
    
    return return_data
    
# read data 
Mio278_data = extract_vars_and_analysis(data=W1E1_278_data, pi_data=W1E1_PI_data)
Mio450_data = extract_vars_and_analysis(data=W1E1_450_data, pi_data=W1E1_PI_data)



apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)
#projection = ccrs.PlateCarree()

fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(25,20), 
                                                      subplot_kw={"projection":projection})

plot_annual_mean(variable="Mean sea level pressure anomalies", data_alt=Mio278_data.get("slp_DJF"), ax=ax1,
                 cmap="RdBu_r", units="hPa", vmax=20, vmin=-20, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(a) MIO 278ppm - PI (DJF)", orientation="horizontal",
                cbar_pos= [0.35, 0.05, 0.25, 0.02], plot_winds=True, data_v=Mio278_data.get("v_DJF"), 
                data_u=Mio278_data.get("u_DJF"), domain="Europe")

plot_annual_mean(variable="Mean sea level pressure anomalies", data_alt=Mio450_data.get("slp_DJF"), ax=ax2,
                 cmap="RdBu_r", units="hPa", vmax=20, vmin=-20, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(b) MIO 450ppm - PI (DJF)", orientation="horizontal"
                , plot_winds=True, data_v=Mio450_data.get("v_DJF"), 
                data_u=Mio450_data.get("u_DJF"), domain="Europe")

plot_annual_mean(variable="Mean sea level pressure anomalies", data_alt=Mio278_data.get("slp_JJA"), ax=ax3,
                 cmap="RdBu_r", units="hPa", vmax=20, vmin=-20, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(c) MIO 278ppm - PI (JJA)", orientation="horizontal"
                , plot_winds=True, data_v=Mio278_data.get("v_JJA"), 
                data_u=Mio278_data.get("u_JJA"), domain="Europe")

plot_annual_mean(variable="Mean sea level pressure anomalies", data_alt=Mio450_data.get("slp_JJA"), ax=ax4,
                 cmap="RdBu_r", units="hPa", vmax=20, vmin=-20, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(d) MIO 450ppm - PI (JJA)", orientation="horizontal",
                plot_winds=True, data_v=Mio450_data.get("v_JJA"), 
                data_u=Mio450_data.get("u_JJA"), domain="Europe")


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
plt.savefig(os.path.join(path_to_plots, "pressure and winds.svg"), format= "svg", bbox_inches="tight", dpi=600)



