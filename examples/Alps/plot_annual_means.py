#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 12:14:06 2021

@author: dboateng

This is example script for using Climat for visualising ECHAM module output (long-term annual means)
The script contains directories of module outputs and path to save plot
Note: it is structured solely for the personal needs of the author, therefore, it must be adopted advisably.
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)

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


# Paths
module_output_main_path = "D:/Datasets/Model_output_pst"

exp_name = "t004_dkrz-mistral_e5w2.3_AMIP_t159l31.6h"    # simulation with present-day simulation (not different from PI simulations)

#years= "1003_1017"
years = "1980_2000"
period = "1m"

GNIP_path = "/home/dboateng/Datasets/GNIP/water_isotopes/ascii/GNIP_STATIONS_EUROPE"
filename = "GNIP_EU_lat_lon_iso.txt"

# reading dataset
control_data, control_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name, years=years,
                                                  period=period)
df_gnip = read_GNIP_data(path=GNIP_path , filename=filename)
#  extracting variables 

temp2 = extract_var(Dataset=control_data , varname="temp2", units="°C")
prec = extract_var(Dataset= control_data , varname="prec", units="mm/month")
d18op = extract_var(Dataset=control_data , varname="d18op", units="per mil", Dataset_wiso= control_wiso)
elev = extract_var(Dataset=control_data , varname="elev", units="m")
u10 = extract_var(Dataset=control_data , varname="u10")
v10 = extract_var(Dataset=control_data , varname="v10")

# contruct annual means 
temp2_alt = compute_lterm_mean(data=temp2, time="annual")
prec_alt = compute_lterm_mean(data=prec, time="annual")
d18op_alt = compute_lterm_mean(data=d18op, time="annual")
elev_alt = compute_lterm_mean(data=elev, time="annual")
u10_alt = compute_lterm_mean(data=u10, time="annual")
v10_alt = compute_lterm_mean(data=v10, time="annual")


#visualising variables and saving

# defining coordinates for gridbox

#west --> 44, 47 N 1, 8 E
lat_w, h_w = 44 , 3  # lat and height
lon_w , w_w = 1, 7  # long and width 

# #east --> 46, 48.5 N 10, 17 E
# lat_e, h_e = 46, 2.5
# lon_e, w_e = 10, 7

#north --> 46.5, 50 N 5, 16 E
lat_n, h_n = 47, 3.5
lon_n, w_n = 5, 11

# south--> 43, 47 N 7.5, 15 E
lat_s, h_s = 43, 4
lon_s, w_s = 8, 7.5

# #isotopic profiles
# lat_A, h_A = 46, 1
# lon_A, w_A = 0, 25

# lat_B, h_B = 40, 14
# lon_B, w_B = 10, 2

# fonts and ploting stlye 
apply_style(fontsize=22, style=None, linewidth=2) 

projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")
fig, ((ax1,ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(20, 15), subplot_kw={"projection":
                                                                                                                  projection})

# elevation
plot_annual_mean(ax=ax1, variable="Elevation", data_alt=elev_alt, cmap=Greys, units="m", vmax=3000, vmin=0, domain="Europe", 
                  levels=22, level_ticks=6, title="[A]", left_labels=True, bottom_labels=False, use_colorbar_default=True)
# adding transects 

# ax1.add_patch(patches.Rectangle(xy =(lon_w, lat_w), width= w_w, height=h_w, ls= "--", color= red, transform = projection, 
#                                 fc="None", lw=2.5,))

# ax1.add_patch(patches.Rectangle(xy =(lon_n, lat_n), width= w_n, height=h_n, ls= "--", color= black, transform = projection, 
#                                 fc="None", lw=2.5))

# ax1.add_patch(patches.Rectangle(xy =(lon_s, lat_s), width= w_s, height=h_s, ls= "--", color= green, transform = projection, 
#                                 fc="None", lw=2.5))

# d18Op
plot_annual_mean(ax=ax2, variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_alt, cmap=YlGnBu, units="‰", vmax=2, vmin=-16, domain="Europe", 
                  levels=22, level_ticks=10, GNIP_data=df_gnip , title="[B]", left_labels=False, bottom_labels=False, use_colorbar_default=True,
                  center =False)


plot_annual_mean(ax=ax3, variable="Temperature", data_alt=temp2_alt, cmap=RdBu_r, units="°C", vmax=25, vmin=-10, domain="Europe", 
                  levels=22, level_ticks=11, title="[C]", left_labels=True, bottom_labels=True, use_colorbar_default=True)

plot_annual_mean(ax=ax4, variable="Precipitation", data_alt=prec_alt, cmap=Blues, units="mm/month", vmax=250, vmin=0, domain="Europe", 
                  levels=22, level_ticks=6, data_u10=u10_alt, data_v10=v10_alt, title="[D]", left_labels=False,
                  bottom_labels=True, use_colorbar_default=True)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.06)

plt.savefig(os.path.join(path_to_store, "fig1.svg"), format= "svg", bbox_inches="tight", dpi=300)
#plt.savefig(os.path.join(path_to_store, "fig1.png"), format= "png", bbox_inches="tight", dpi=300)


# fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(13, 13), subplot_kw={"projection":  projection})

# # d18Op
# plot_annual_mean(ax=ax1, variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_alt, cmap=YlGnBu, units="‰", vmax=2, vmin=-16, domain="Europe", 
#                  levels=22, level_ticks=10, GNIP_data=df_gnip , title=None, left_labels=True, bottom_labels=True, use_colorbar_default=True)

# fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
# plt.tight_layout() 
# plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.06)

# plt.savefig(os.path.join(path_to_store, "egu3.svg"), format= "svg", bbox_inches="tight", dpi=300)
# plt.savefig(os.path.join(path_to_store, "egu3.png"), format= "png", bbox_inches="tight", dpi=300)
