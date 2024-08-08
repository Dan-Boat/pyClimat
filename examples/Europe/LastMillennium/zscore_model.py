# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 16:03:59 2024

@author: dboateng
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_transect, compute_lterm_diff, compute_lterm_mean
from pyClimat.variables import extract_var
from pyClimat.utils import extract_region


from path_to_data_lm import *

main_path = "D:/Datasets/iGCM_datasets/"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"

filenames = ["CESM", "ECHAM5", "GISS", "HADCM3", "CCSM"]



def extract_climatologies(filename, time="annual", season=None):
    
    filename_data = filename + "_wiso_vars.nc"
    
    
    d18O = read_from_path(path=main_path, filename=filename_data, decode=True,
                      varname="d18O")  
    
    data  = d18O.sel(time= slice("0950-01", "1950-01"))
    
    mean = data.mean(dim="time", skipna=True)
    std = data.std(dim="time", skipna=True)
    zscore = (data - mean) / std
    
    
    date_range_mca = d18O.sel(time= slice("0950-01", "1250-01")).time
    date_range_lia = d18O.sel(time= slice("1650-01", "1850-01")).time
    date_range_pi = d18O.sel(time= slice("1851-01", "1950-01")).time
        
        
    
    d18O_mca = compute_lterm_mean(data=zscore, time=time, time_range=date_range_mca, 
                                  season=season)
    d18O_lia = compute_lterm_mean(data=zscore, time=time, time_range=date_range_lia, 
                                  season=season)
    
    d18O_pi = compute_lterm_mean(data=zscore, time=time, time_range=date_range_pi, 
                                  season=season)
    
    return d18O_mca, d18O_lia, d18O_pi


# test function 
# cesm_mca, cesm_lia, cesm_diff = extract_climatologies(filename="CESM", time="annual")
# giss_mca, giss_lia, giss_diff = extract_climatologies(filename="GISS", time="annual")
data_mca, data_lia, data_pi = extract_climatologies(filename="HADCM3", time="annual")


def plot_d18Op_global(data_mca, data_lia, data_pi, axes=None, fig=None, axes_cbar=True, labels=None):
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    if axes is None:
        fig,(ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(30,8), subplot_kw={"projection":projection})
        axes = [ax1, ax2, ax3]
    
    
    
            
    plot_annual_mean(variable="Mean(Z-score)", data_alt=data_mca, ax=axes[0],
                     cmap=RdBu_r, units="‰", vmax=0.3, vmin=-0.3, 
                    levels=22, level_ticks=11, add_colorbar=False, bottom_labels=True,
                    left_labels=True, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[0], center=True, domain="NH Wide")
    
   
    plot_annual_mean(variable="Mean(Z-score)", data_alt=data_lia, ax=axes[1],
                     cmap=RdBu_r, units="‰", vmax=0.3, vmin=-0.3, 
                    levels=22, level_ticks=11, add_colorbar=False, 
                    bottom_labels=True,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[1], center=True, domain="NH Wide")
    
    plot_annual_mean(variable="Mean(Z-score)", data_alt=data_pi, ax=axes[2],
                     cmap=RdBu_r, units="‰", vmax=0.3, vmin=-0.3, 
                    levels=22, level_ticks=9, add_colorbar=axes_cbar, cbar_pos= [0.35, 0.02, 0.35, 0.02], 
                    orientation="horizontal", bottom_labels=True,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[2], center=True, domain="NH Wide",
                    label_format="%.1f")
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10)
    plt.savefig(os.path.join(path_to_plots, "hadcm3_mca_lia_pi.png"), format= "png", bbox_inches="tight", dpi=300)
    
    
    
plot_d18Op_global(data_mca, data_lia, data_pi, labels = ["(a) MCA   [iHadCM3]", "(b) LIA", "(c) PI"])