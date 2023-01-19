#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:45:48 2023

@author: dboateng

This script generates the plot of precipitation and temperature for the Monsoon region using both ERA5 and ECHAM5-wiso (1980-200)
"""
import os 
import pandas as pd 

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean 
from pyClimat.data import read_ERA_processed, read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff, extract_var

#useful functions
#define variables
t2m = "temp2"
prec = "prec"
v10 = "v10"
u10 = "u10"


def extract_only_prec_temp_winds(data_surface):
    
    data_t2m = extract_var(Dataset=data_surface, varname=t2m, units="°C")
    
    data_prec = extract_var(Dataset=data_surface, varname=prec, units="mm/month")
    
    data_v10 = extract_var(Dataset=data_surface, varname=v10) # default in m/s
    data_u10 = extract_var(Dataset=data_surface, varname=u10) # default in m/s
    
    
    return data_prec, data_t2m, data_v10, data_u10



# define paths 

path_to_plots = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"
ERA5_path = "/home/dboateng/Datasets/ERA5/monthly_1950_2021/"

ERA5_tp_path = os.path.join(ERA5_path, "tp_monthly.nc")
ERA5_t2m_path = os.path.join(ERA5_path, "t2m_monthly.nc")
ERA5_v10_path = os.path.join(ERA5_path, "v10_monthly.nc")
ERA5_u10_path = os.path.join(ERA5_path, "u10_monthly.nc")

#define time of amip 
from1980to2000 = pd.date_range(start="1980-01-01", end="2000-12-31", freq="MS")


#read in postprocessed and analysed data 
ERA5_t2m = read_ERA_processed(path=ERA5_t2m_path, varname="t2m")   - 273.15 #°C
ERA5_tp = read_ERA_processed(path=ERA5_tp_path, varname="tp") * 1000 * 30  #mm/month
ERA5_v10 = read_ERA_processed(path=ERA5_v10_path, varname="v10") #m/s
ERA5_u10 = read_ERA_processed(path=ERA5_u10_path, varname="u10") #m/s

main_path = "/home/dboateng/Model_output_pst"
exp_name = "t004_dkrz-mistral_e5w2.3_AMIP_t159l31.6h"    # simulation with present-day simulation (not different from PI simulations)

years = "1980_2000"
period = "1m"

PD_data = read_ECHAM_processed(main_path=main_path , exp_name= exp_name, years=years,
                                                  period=period, read_wiso=False)



PD_prec, PD_t2m, PD_v10, PD_u10 = extract_only_prec_temp_winds(PD_data)

#ERA
ERA5_t2m_alt = compute_lterm_mean(data=ERA5_t2m, time="month", month="JJAS", time_range=from1980to2000)
ERA5_tp_alt = compute_lterm_mean(data=ERA5_tp, time="month", month="JJAS", time_range=from1980to2000)
ERA5_v10_alt = compute_lterm_mean(data=ERA5_v10, time="month", month="JJAS", time_range=from1980to2000)
ERA5_u10_alt = compute_lterm_mean(data=ERA5_u10, time="month", month="JJAS", time_range=from1980to2000)


#PD
PD_t2m_alt = compute_lterm_mean(data=PD_t2m, time="month", month="JJAS")
PD_prec_alt = compute_lterm_mean(data=PD_prec, time="month", month="JJAS")
PD_v10_alt = compute_lterm_mean(data=PD_v10, time="month", month="JJAS")
PD_u10_alt = compute_lterm_mean(data=PD_u10, time="month", month="JJAS")

#ERA5-ECHAM5
# interpolate ERA to ECHAM resolution
ERA5_t2m_alt = ERA5_t2m_alt.rename({"longitude": "lon", "latitude":"lat"})
ERA5_tp_alt = ERA5_tp_alt.rename({"longitude": "lon", "latitude":"lat"})


ERA5_t2m_alt_interp = ERA5_t2m_alt.interp(lat=PD_t2m_alt.lat).interp(lon=PD_t2m_alt.lon)
ERA5_tp_alt_interp = ERA5_tp_alt.interp(lat=PD_prec_alt.lat).interp(lon=PD_prec_alt.lon)

t2m_diff = PD_t2m_alt - ERA5_t2m_alt_interp
prec_diff = PD_prec_alt - ERA5_tp_alt_interp



# set up plot function for specific variable

apply_style(fontsize=22, style=None, linewidth=2) 

projection = ccrs.PlateCarree()
fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(28, 15), subplot_kw={"projection": projection})


plot_annual_mean(ax=ax1, fig=fig, variable="Precipitation", data_alt=ERA5_tp_alt, cmap=YlGnBu, units="mm/month", vmax=400, vmin=50, domain="West Africa", 
                  levels=22, level_ticks=6, title="ERA5", left_labels=True, bottom_labels=False, 
                  add_colorbar=True, cbar_pos = [0.20, 0.52, 0.25, 0.02], orientation= "horizontal")

plot_annual_mean(ax=ax2, fig=fig, variable="Precipitation", data_alt=PD_prec_alt, cmap=YlGnBu, units="mm/month", vmax=400, vmin=50, domain="West Africa", 
                  levels=22, level_ticks=6, title="ECHAM5-wiso (PD)", left_labels=False, bottom_labels=False, 
                  add_colorbar=False)

plot_annual_mean(ax=ax3, fig=fig, variable="Precipitation anomalies", data_alt=prec_diff, cmap=BrBG, units="mm/month", vmax=150, vmin=-150, domain="West Africa", 
                  levels=22, level_ticks=11, title="ECHAM5-wiso (PD) - ERA5", left_labels=False, bottom_labels=False, 
                  add_colorbar=True, cbar_pos = [0.65, 0.52, 0.25, 0.02], orientation= "horizontal")


plot_annual_mean(ax=ax4, fig=fig, variable="Temperature", data_alt=ERA5_t2m_alt, cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                  levels=22, level_ticks=11, title="ERA5", left_labels=True, bottom_labels=True, 
                  add_colorbar=True, cbar_pos = [0.20, 0.05, 0.25, 0.02], orientation= "horizontal")

plot_annual_mean(ax=ax5, fig=fig, variable="Temperature", data_alt=PD_t2m_alt, cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                  levels=22, level_ticks=11, title="ECHAM5-wiso (PD)", left_labels=False, bottom_labels=True, 
                  add_colorbar=False)

plot_annual_mean(ax=ax6, fig=fig, variable="Temperature anomalies", data_alt=t2m_diff, cmap=RdBu_r, units="°C", vmax=10, vmin=-10, domain="West Africa", 
                  levels=22, level_ticks=11, title="ECHAM5-wiso (PD) - ERA5", left_labels=False, bottom_labels=True, 
                  add_colorbar=True, cbar_pos = [0.65, 0.05, 0.25, 0.02], orientation= "horizontal")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, "compare_era_and_echam.svg"), format= "svg", bbox_inches="tight", dpi=300)

plt.show()    
    