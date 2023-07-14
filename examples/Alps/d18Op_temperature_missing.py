# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 17:04:31 2023

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
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var


module_output_main_path = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/plots"


# Path to experiments

exp_name_aw100e100 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_aw100e0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_aw100e200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"

# west set up 
exp_name_aw200e100 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
exp_name_aw200e0 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
exp_name_aw200e200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"
exp_name_aw100e150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"


# for supplementary (same but for annual)
years= "1003_1017"
period = "1m"


# reading dataset
aw100e100_data, aw100e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e100, years=years,
                                                  period=period)
aw100e0_data, aw100e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e0, years=years,
                                                  period=period)
aw100e200_data, aw100e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e200, years=years,
                                                  period=period)
aw200e100_data, aw200e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e100, years=years,
                                                  period=period)
aw200e0_data, aw200e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e0, years=years,
                                                  period=period)
aw200e200_data, aw200e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e200, years=years,
                                                  period=period)
aw100e150_data, aw100e150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e150, years=years,
                                                  period=period)


def extract_and_analysis(exp_data, exp_wiso, W1E1_data, W1E1_wiso, diff="missing",
                                    d18Op_lapse_rate=-0.002, temp_lapse_rate=-0.0056):
    
    # extract d18op
    exp_d18op_data = extract_var(Dataset=exp_data, varname="d18op", units="per mil", Dataset_wiso=exp_wiso)
    w1e1_d18op_data = extract_var(Dataset=W1E1_data, varname="d18op", units="per mil", Dataset_wiso=W1E1_wiso)
    
    exp_temp_data = extract_var(Dataset=exp_data, varname="temp2", units="°C")
    w1e1_temp_data = extract_var(Dataset=W1E1_data, varname="temp2", units="°C")
    
    
    
    #extract elevation 
    exp_elev_data = extract_var(Dataset=exp_data, varname="elev", units="m",)
    w1e1_elev_data = extract_var(Dataset=W1E1_data, varname="elev", units="m")
    
    # changes in elevation 
    elev_change = exp_elev_data - w1e1_elev_data
    
    # lat_range = (elev_change.lat > 50) & (elev_change.lat < 40)
    # lon_range = (elev_change.lon > 20) & (elev_change.lon < 0)
    
    #elev_change_modified = elev_change.where((lat_range & lon_range), 0,)
    
    expected_d18op_change = elev_change * d18Op_lapse_rate # per mil 
    
    simulated_d18op_change = exp_d18op_data - w1e1_d18op_data
    
    d18Op_missing = simulated_d18op_change - expected_d18op_change
    
    
    expected_temp_change = elev_change * temp_lapse_rate # per mil 
   
    simulated_temp_change = exp_temp_data - w1e1_temp_data
   
    temp_missing = simulated_temp_change - expected_temp_change
    
    # calculate annual means of missing
    d18Op_missing_alt = compute_lterm_mean(data=d18Op_missing, time="annual")
    
    temp_missing_alt = compute_lterm_mean(data=temp_missing, time="annual")
    
    temp_expected_alt = compute_lterm_mean(data=expected_temp_change, time="annual")
    
    return_data = {"d18Op": d18Op_missing_alt, "temp": temp_missing_alt, "temp_exp": temp_expected_alt}
    
    
    return return_data


labels = ["W2E1", "W1E0", "W2E0", "W1E2",]
data_all = [aw200e100_data,  aw100e0_data, aw200e0_data, aw100e200_data]
wiso_all = [aw200e100_wiso,  aw100e0_wiso, aw200e0_wiso, aw100e200_wiso]

extracts = {}


for i,topo in enumerate(labels):
    extracts[topo] = extract_and_analysis(exp_data=data_all[i], exp_wiso=wiso_all[i],
                                          W1E1_data=aw100e100_data, W1E1_wiso=aw100e100_wiso)
    
    
projection = ccrs.PlateCarree()


apply_style(fontsize=28, style=None, linewidth=2.5)

fig, ((ax1,ax2),(ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(20, 15),
                                                      subplot_kw={"projection": projection})
axes = [ax1, ax2, ax3, ax4]


for i,label in enumerate(labels):
    if i == 0:
        
        plot_annual_mean(variable="Temperature difference", data_alt=extracts[label].get("temp"), ax=axes[i],
                          cmap=RdBu_r, units="°C", vmax=2.5, vmin=-2.5, 
                        levels=18, level_ticks=6, add_colorbar=True, cbar_pos= [0.30, 0.05, 0.45, 0.02], 
                        orientation="horizontal", plot_coastlines=True, bottom_labels=True,
                        left_labels=True, fig=fig, plot_borders=False, plot_projection=projection, 
                        domain="Europe", title=label, label_format="%.1f")
        
    else:
        plot_annual_mean(variable="Temperature difference", data_alt=extracts[label].get("temp"), ax=axes[i],
                          cmap=RdBu_r, units="°C", vmax=2.5, vmin=-2.5, 
                        levels=18, level_ticks=6, add_colorbar=False, plot_coastlines=True, bottom_labels=True,
                        left_labels=True, fig=fig, plot_borders=False, plot_projection=projection, domain="Europe", 
                        max_pvalue=0.1, title=label)
        
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
plt.savefig(os.path.join(path_to_plots, "temp_correction.svg"), format= "svg", bbox_inches="tight", dpi=600)



apply_style(fontsize=28, style=None, linewidth=2.5)

fig, ((ax1,ax2),(ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(20, 15),
                                                      subplot_kw={"projection": projection})
axes = [ax1, ax2, ax3, ax4]


for i,label in enumerate(labels):
    if i == 0:
        
        plot_annual_mean(variable="$\delta^{18}$Op vs SMOW difference", data_alt=extracts[label].get("d18Op"), ax=axes[i],
                          cmap="RdBu", units="‰", vmax=4, vmin=-4, 
                        levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.30, 0.05, 0.45, 0.02], 
                        orientation="horizontal", plot_coastlines=True, bottom_labels=True,
                        left_labels=True, fig=fig, plot_borders=False, plot_projection=projection, 
                        domain="Europe", title=label)
        
    else:
        plot_annual_mean(variable="$\delta^{18}$Op vs SMOW difference", data_alt=extracts[label].get("d18Op"), ax=axes[i],
                          cmap="RdBu", units="‰", vmax=4, vmin=-4, 
                        levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=True,
                        left_labels=True, fig=fig, plot_borders=False, plot_projection=projection, domain="Europe", 
                        title=label)
        
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
plt.savefig(os.path.join(path_to_plots, "d18Op_correction.svg"), format= "svg", bbox_inches="tight", dpi=600)

