#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 17:19:18 2022

@author: dboateng
"""


import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os 
from pyClimat.data import read_from_path, read_ERA_processed
from pyClimat.analysis import extract_var, compute_lterm_mean

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean

tp_path = "/home/dboateng/Datasets/ERA5/monthly_1950_2021/t2m_monthly.nc"
main_path = "/home/dboateng/Model_output_pst"
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")

filename_lterm = "1003_1017_1m_mlterm.nc"

path_to_store = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

# read long-term means

PI_data = read_from_path(pi_path, filename_lterm)
t2m = read_ERA_processed(path=tp_path, varname="t2m")   - 273.15 #°C

pi_temp2 = extract_var(Dataset=PI_data, varname="temp2", units="°C")


PI_temp2_alt = compute_lterm_mean(data=pi_temp2, time="month", month="JJAS")
t2m_alt = compute_lterm_mean(data=t2m, time="month", month="JJAS")

def plot_PI_ERA5_t2m():
    apply_style(fontsize=22, style=None, linewidth=2) 
    
    projection = ccrs.PlateCarree()
    fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(20, 13), subplot_kw={"projection":
                                                                                                                      projection})
    
    # add patch (Sahel: 10-20 N, 20W - 30E; coast of Guinea: 5-10N, 20W-30E, Sahara region: 20-30N, 20W-30E)
    
    #sahara --> 30, 20 N 20W, 30 E
    
    lat_sh, h_sh = 20 , 10  # lat and height
    lon_sh , w_sh = -20, 50  # long and width 
    
    #sahel --> 20, 10 N 20W, 30 E
    
    lat_sl, h_sl = 10 , 9.5  # lat and height
    lon_sl , w_sl = -20, 50  # long and width 
    
    #guinea -->5, 10 N 20W, 30 E
    
    lat_g, h_g = 5 , 4.5  # lat and height
    lon_g , w_g = -20, 50  # long and width 
    
    
    plot_annual_mean(ax=ax1, variable="Temperature", data_alt=PI_temp2_alt, cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                      levels=22, level_ticks=11, title="[A] PI-ECHAM5", left_labels=True, bottom_labels=True, use_colorbar_default=True,
                      )
    
    
    plot_annual_mean(ax=ax2, variable="Temperature", data_alt=t2m_alt, cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                      levels=22, level_ticks=11, title="[B] ERA5", left_labels=False, bottom_labels=True, 
                      add_colorbar=False, )
    
 
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.06)
    
    plt.savefig(os.path.join(path_to_store, "t2m_echam_era.svg"), format= "svg", bbox_inches="tight", dpi=300)
    
plot_PI_ERA5_t2m()