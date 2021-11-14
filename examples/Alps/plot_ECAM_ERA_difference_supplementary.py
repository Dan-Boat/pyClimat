#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 18:44:18 2021

@author: dboateng
This script plots the difference between ERA-Interim and PI control run
"""
# importimng package 

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

t2m_path = "/home/dboateng/Model_output_pst/era/t2m_monthly.nc"
tp_path = "/home/dboateng/Model_output_pst/era/tp_monthly.nc"

# ERA-interim 
t2m = read_ERA_processed(path=t2m_path, varname="t2m") - 273.15 #°C
tp = read_ERA_processed(path=tp_path, varname="tp") * 1000  #mm/month

t2m = t2m.rename({"longitude": "lon", "latitude":"lat"})

#setting paths 
module_output_main_path = "/home/dboateng/Model_output_pst"
years = "1003_1017"

a002 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

#reading 
a002_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name=a002, years=years, read_wiso=False, period="1m")

t2m = t2m.interp(lat=a002_data.lat).interp(lon=a002_data.lon)
tp = tp.interp(lat=a002_data.lat).interp(lon=a002_data.lon)

#extracting temp and prec
a002_temp2 = extract_var(Dataset=a002_data , varname="temp2", units="°C")
a002_prec = extract_var(Dataset= a002_data , varname="prec", units="mm/month")

# computing long-term difference
a002_temp2_diff = compute_lterm_diff(data_control=t2m , data_main=a002_temp2, time="season", season_calendar="standard")
a002_prec_diff = compute_lterm_diff(data_control=tp , data_main=a002_prec, time="season", season_calendar="standard")

#visualising variables and saving
projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")
fig, ((ax1,ax2),(ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(17, 10), subplot_kw={"projection":projection})

plot_seasonal_mean(variable="Temperature", data_slt=a002_temp2_diff , cmap=RdBu_r, units="°C", seasons=["JJA", "DJF"], 
                   axes=[ax1,ax2], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.56, 0.02, 0.30],title=True,
                   season_label = ["[A]              Summer (JJA)", "[B]               Winter (DJF)"], fig_title= "Alps 100% - ERA-Interim")

plot_seasonal_mean(variable="Precipitation", data_slt=a002_prec_diff , cmap=RdBu, units="mm/month", seasons=["JJA", "DJF"], 
                   axes=[ax3,ax4], fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=6, cbar_pos = [0.90, 0.13, 0.02, 0.3],
                   title=True, season_label=["[C]", "[D]"])
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "figS2.png"), format= "png", bbox_inches="tight", dpi=600)
plt.show()