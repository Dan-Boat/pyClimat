# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 19:34:29 2023

@author: dboateng
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
import matplotlib.patches as patches

from pyClimat.analysis import extract_profile
from pyClimat.plot_utils import *
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

path_to_gtopo = "D:/Datasets/ECHAM5/Inputs/CTL_gtopo/global_gtopo30.nc"

data = xr.open_dataset(path_to_gtopo)

elev_data = data.z

def extract_topo_profile(maxlon=20, minlon=-5, maxlat=46.8, minlat=46, 
                          dim="lon"):
    
    extract_lon = extract_profile(data = elev_data, maxlon=25, minlon=-5, maxlat=46.5, minlat=46, 
                              dim="lon", to_pandas=True)
    
    extract_lat = extract_profile(data = elev_data, maxlon=10.5, minlon=10, maxlat=54, minlat=40, 
                              dim="lat", to_pandas=True)
    
    extract_lon = extract_lon/1000
    
    extract_lat = extract_lat/1000
    
    apply_style(fontsize=28, style="seaborn-talk", linewidth=3,)
    
    fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(30, 15), sharey=False)
    
    extract_lon.plot(ax=ax1, linestyle="-", color="black", linewidth=2, legend=False)
    extract_lon.plot(kind="area", color="black", alpha=0.1, legend=False, stacked=False,
                   ax=ax1) 
    
    ax1.set_ylim(0, 12)
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.set_ylabel("Elevation [km]", fontweight="bold", fontsize=28)
    
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.tick_params(axis='y', which='major', length=8, width=3)
    ax1.tick_params(axis='y', which='minor', length=4, width=1.5)
    
    
    extract_lat.plot(ax=ax2, linestyle="-", color="black", linewidth=2, legend=False)
    extract_lat.plot(kind="area", color="black", alpha=0.1, legend=False, stacked=False,
                   ax=ax2) 
    
    ax2.set_ylim(0, 12)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.tick_params(axis='y', which='major', length=8, width=3)
    ax2.tick_params(axis='y', which='minor', length=4, width=1.5)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.90, top=0.97, bottom=0.05, wspace=0.12)
    plt.savefig(os.path.join(path_to_plots, "elevation_gtopo_profile.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    

df = extract_topo_profile()