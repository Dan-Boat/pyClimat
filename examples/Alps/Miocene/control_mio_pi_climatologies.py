# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 11:27:11 2023

@author: dboateng

Plot the annual means of the control experiments (and summer means if required)
"""
# import models
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
from pyClimat.analysis import extract_var, compute_lterm_mean


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"


# read data

CTL_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

years = "1003_1017"
period = "1m"

CTL_data, CTL_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=CTL_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

def extract_relevant_vars(data, wiso):
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    u10 = extract_var(Dataset=data , varname="u10")
    v10 = extract_var(Dataset=data , varname="v10")
    
    # compute annual means
    temp2_alt = compute_lterm_mean(data=temp2, time="annual")
    prec_alt = compute_lterm_mean(data=prec, time="annual")
    d18op_alt = compute_lterm_mean(data=d18op, time="annual")
    u10_alt = compute_lterm_mean(data=u10, time="annual")
    v10_alt = compute_lterm_mean(data=v10, time="annual")
   
    
    
    return temp2_alt, prec_alt, d18op_alt, u10_alt, v10_alt


# extract slm for miocene

slm_mio = W1E1_278_data.slm[0]
#xtract variables

CTL_temp2, CTL_prec, CTL_d18op, CTL_u10, CTL_v10 = extract_relevant_vars(data=CTL_data, wiso=CTL_wiso)
W1E1_278_temp2, W1E1_278_prec, W1E1_278_d18op, W1E1_278_u10, W1E1_278_v10 = extract_relevant_vars(data=W1E1_278_data, wiso=W1E1_278_wiso)
W1E1_450_temp2, W1E1_450_prec, W1E1_450_d18op, W1E1_450_u10, W1E1_450_v10 = extract_relevant_vars(data=W1E1_450_data, wiso=W1E1_450_wiso)


# plot
apply_style(fontsize=22, style=None, linewidth=2) 
    
projection = ccrs.PlateCarree()
fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 13), subplot_kw={"projection":  projection})

plot_annual_mean(ax=ax1, variable="Precipitation", data_alt=CTL_prec, cmap=YlGnBu, units="mm/month", vmax=200, vmin=50, domain="Europe", 
                  levels=22, level_ticks=6, title="CTL (PI)", left_labels=True, bottom_labels=True, add_colorbar=True, cbar_pos = [0.35, 0.25, 0.25, 0.02],
                  plot_projection=projection, plot_winds=True, data_u=CTL_u10, data_v=CTL_v10,  orientation= "horizontal", fig=fig, plot_coastlines=True)

plot_annual_mean(ax=ax2, variable="Precipitation", data_alt=W1E1_278_prec, cmap=YlGnBu, units="mm/month", vmax=200, vmin=50, domain="Europe", 
                  levels=22, level_ticks=6, title="W1E1 (MIO 278ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                  plot_winds=True, data_u=W1E1_278_u10, data_v=W1E1_278_v10, fig=fig, plot_coastlines=False, 
                  sea_land_mask=slm_mio)




