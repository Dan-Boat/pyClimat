# -*- coding: utf-8 -*-
"""
Created on Wed May 17 19:04:45 2023

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
from pyClimat.analysis import extract_var, compute_lterm_mean, compute_lterm_diff


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"
W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

# reading data 
# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_PI_data, W1E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)


# extract values and compute long-term means

def extract_vars_and_analysis(data, wiso, pi_data, pi_wiso):
    
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    
    temp2_pi = extract_var(Dataset=pi_data , varname="temp2", units="°C")
    prec_pi = extract_var(Dataset= pi_data , varname="prec", units="mm/month")
    d18op_pi = extract_var(Dataset=pi_data , varname="d18op", units="per mil", Dataset_wiso= pi_wiso)
    
    #compute climatologies difference
    
    temp2_diff = compute_lterm_diff(data_control=temp2_pi, data_main=temp2, time="annual")
    prec_diff = compute_lterm_diff(data_control=prec_pi, data_main=prec, time="annual")
    d18op_diff = compute_lterm_diff(data_control=d18op_pi, data_main=d18op, time="annual")
    
    
    return_data = {"temperature":temp2_diff, "precipitation":prec_diff, "d18Op":d18op_diff,}
    
    return return_data
    
# read data 
Mio278_data = extract_vars_and_analysis(data=W1E1_278_data, wiso=W1E1_278_wiso, pi_data=W1E1_PI_data, 
                                        pi_wiso=W1E1_PI_wiso)
Mio450_data = extract_vars_and_analysis(data=W1E1_450_data, wiso=W1E1_450_wiso, pi_data=W1E1_PI_data, 
                                        pi_wiso=W1E1_PI_wiso)



apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)

fig,((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(24,18), subplot_kw={"projection":projection})

plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=Mio278_data.get("d18Op"), ax=ax1,
                 cmap="PRGn", units="‰", vmax=5, vmin=-5, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(a) MIO 278ppm - PI", orientation="horizontal",
                cbar_pos= [0.05, 0.05, 0.25, 0.02])

plot_annual_mean(variable="Precipitation anomalies", data_alt=Mio278_data.get("precipitation"), ax=ax2,
                 cmap="BrBG", units="mm/month", vmax=150, vmin=-150, 
                levels=22, level_ticks=7, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(b) MIO 278ppm - PI", orientation="horizontal",
                cbar_pos= [0.35, 0.05, 0.25, 0.02])

plot_annual_mean(variable="Temperature anomalies", data_alt=Mio278_data.get("temperature"), ax=ax3,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(c) MIO 278ppm - PI", orientation="horizontal",
                cbar_pos= [0.65, 0.05, 0.25, 0.02])


plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=Mio450_data.get("d18Op"), ax=ax4,
                 cmap="PRGn", units="‰", vmax=5, vmin=-5, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(e) MIO 450ppm - PI", center=False,)

plot_annual_mean(variable="Precipitation anomalies", data_alt=Mio450_data.get("precipitation"), ax=ax5,
                 cmap="BrBG", units="mm/month", vmax=150, vmin=-150, 
                levels=22, level_ticks=7, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(f) MIO 450ppm - PI",)

plot_annual_mean(variable="Temperature anomalies", data_alt=Mio450_data.get("temperature"), ax=ax6,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(g) MIO 450ppm - PI",)


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
plt.savefig(os.path.join(path_to_plots, "global_anomalies.svg"), format= "svg", bbox_inches="tight", dpi=600)
plt.show()
