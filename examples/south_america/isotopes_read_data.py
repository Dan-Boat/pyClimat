# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 09:16:58 2024

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl 
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean

from pyClimat.utils import extract_indices_around


def calculate_regional_means(ds, lon_target, lat_target, radius_deg,):
    """
    Calculate regional means around a specific longitude and latitude location
    with a given radius for a NetCDF dataset using xarray.
    """
    # Find indices of the nearest grid point to the target location
    
    if hasattr(ds, "longitude"):
        ds = ds.rename({"longitude":"lon", "latitude":"lat"})
        
    ds = ds.assign_coords({"lon": (((ds.lon + 180) % 360) - 180)})
    
    indices = extract_indices_around(ds, lat_target, lon_target, radius_deg)
    
    regional_mean = ds.isel(lat=indices[0], lon=indices[1]).mean(dim=("lon", "lat")).data
        
    return np.float64(regional_mean)



path_to_data = "D:/Datasets/Model_output_pst/ECHAM5-wiso_7k"

filename = "Hol-T_echam5_wiso_wisoaprt_d_0004_to_7000_2_yearmean.nc"

filename_mon = "Hol-T_echam5_wiso_wisoaprt_d_0004_to_7000_2.nc"

#dates = xr.cftime_range(start="1850", periods=83957, freq="MS", calendar="noleap")

data = read_from_path(path=path_to_data, filename=filename_mon, varname= "wisoaprt_d", decode=True)

data["time"] = dates


# compute means

data_seasonal_means = compute_lterm_mean(data=data, time="season")

data_seasonal_means_JJA = data_seasonal_means.sel(season="JJA")

data_seasonal_means_JJA = data_seasonal_means_JJA[0, :, :]



# extract point

location_data = calculate_regional_means(ds=data, lon_target=-56, lat_target=-20, radius_deg=100)


# plotting 
# projection = ccrs.Robinson(central_longitude=0, globe=None)

# apply_style(fontsize=28, style=None, linewidth=2.5) 
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,12), subplot_kw={"projection":projection})

# plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=data_seasonal_means_JJA, ax=ax,
#               cmap=RdYlBu, units="‰", vmax=2, vmin=-20, 
#             levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
#             left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="Mid-Holocene", 
#             orientation="horizontal",
#             cbar_pos= [0.28, 0.05, 0.45, 0.02], center=False, domain="South America")


# fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
# plt.tight_layout() 
# plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)


