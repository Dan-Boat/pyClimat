#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 16:20:22 2022

@author: dboateng

This script intend to analyse the simulated changes in the Miocene compared to the PI-industrial
Interested variables: sea level pressure, precipitation, winds, geopotential height, and the integrated vertical moisture or precipitable water

! The script heavily depend on the pyClimat package which can installed with pip install pyclimat (is you face issues with cartopy, please check how to install it 
                                                                                                   on your system or contact me for help)
"""

#loading of modules 
import os 
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


from pyClimat.data import read_from_path
from pyClimat.analysis import extract_var, compute_lterm_diff, compute_lterm_mean
from pyClimat.plot_utils import *
from pyClimat.plots import plot_seasonal_mean


path_to_store = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/south_america/plots"
# define paths 
main_path = "/home/dboateng/Model_output_pst"

# reading files 

mh_path = os.path.join(main_path, "MH", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")

filename_lterm = "1003_1017_1m_mlterm.nc"

MH_data = read_from_path(mh_path, filename_lterm)
PI_data = read_from_path(pi_path, filename_lterm)

filename_plev_lterm = "1003_1017_1m_mlterm_plev.nc"


PI_plev_data = read_from_path(pi_path, filename_plev_lterm)
MH_plev_data = read_from_path(mh_path, filename_plev_lterm)
 

#calculating the anomalies

# define variables
t2m = "temp2"
prec = "prec"
mslp = "slp"
geosp = "geopoth"



def extract_all_variables(data_surface, data_plev):
    
    data_t2m = extract_var(Dataset=data_surface, varname=t2m, units="°C")
    data_prec = extract_var(Dataset=data_surface, varname=prec, units="mm/month")
    data_qvi = extract_var(Dataset=data_surface, varname="qvi", units="kg/m²")
    data_slp = extract_var(Dataset=data_plev, varname=mslp, units="hPa")
    data_v = extract_var(Dataset=data_plev, varname="v", lev_units="hPa")
    data_u = extract_var(Dataset=data_plev, varname="u", lev_units="hPa")
    data_geosp = extract_var(Dataset=data_plev, varname=geosp, lev_units="hPa")
    data_v500 = data_v.sel(lev=500)
    data_u500 = data_u.sel(lev=500)
    data_geosp_500 = data_geosp.sel(lev=500)

    return data_t2m, data_prec, data_qvi, data_slp, data_v500, data_u500, data_geosp_500 

PI_t2m, PI_prec, PI_qvi, PI_slp, PI_v500, PI_u500, PI_geosp_500  = extract_all_variables(PI_data, PI_plev_data)
MH_t2m, MH_prec, MH_qvi, MH_slp, MH_v500, MH_u500, MH_geosp_500  = extract_all_variables(MH_data, MH_plev_data)

# differience

t2m_diff = compute_lterm_diff(data_control=PI_t2m, data_main=MH_t2m, time="season", season_calendar="standard")
prec_diff = compute_lterm_diff(data_control=PI_prec, data_main=MH_prec, time="season", season_calendar="standard")
slp_diff = compute_lterm_diff(data_control=PI_slp, data_main=MH_slp, time="season", season_calendar="standard")
geosp_diff = compute_lterm_diff(data_control=PI_geosp_500, data_main=MH_geosp_500, time="season", season_calendar="standard")
qvi_diff = compute_lterm_diff(data_control=PI_qvi, data_main=MH_qvi, time="season", season_calendar="standard")


#means

PI_t2m_alt = compute_lterm_mean(data=PI_t2m, time="season", season_calendar="standard")
PI_prec_alt =  compute_lterm_mean(data=PI_prec, time="season", season_calendar="standard")
PI_slp_alt =  compute_lterm_mean(data=PI_slp, time="season", season_calendar="standard")
PI_qvi_alt =  compute_lterm_mean(data=PI_qvi, time="season", season_calendar="standard")
PI_geop_500_alt =  compute_lterm_mean(data=PI_geosp_500, time="season", season_calendar="standard")
PI_v500_alt =  compute_lterm_mean(data=PI_v500, time="season", season_calendar="standard")
PI_u500_alt =  compute_lterm_mean(data=PI_u500, time="season", season_calendar="standard")

MH_v500_alt =  compute_lterm_mean(data=MH_v500, time="season", season_calendar="standard")
MH_u500_alt =  compute_lterm_mean(data=MH_u500, time="season", season_calendar="standard")


#plotting
projection = ccrs.PlateCarree()
apply_style(fontsize=22, style=None, linewidth=2) 

fig, ((ax1,ax2), (ax3, ax4) )= plt.subplots(nrows = 2, ncols = 2, figsize=(25, 23), subplot_kw={"projection":
                                                                                                                      projection})

plot_seasonal_mean(variable="Precipitation", data_slt=prec_diff, cmap=BrBG, units="mm/month", seasons=["JJA", "DJF"], 
                   axes=[ax3, ax4], vmin=-180, vmax=180, levels=22, level_ticks=11, add_colorbar=True, 
                   cbar_pos = [0.90, 0.30, 0.02, 0.25], fig=fig, domain="South America", orientation="vertical",
                   bottom_labels=True, season_label=["[C] MH - PI", "[D] MH - PI"], title=True)

plot_seasonal_mean(variable="Precipitation", data_slt=PI_prec_alt, cmap=YlGnBu, units="mm/month", seasons=["JJA", "DJF"], 
                   axes=[ax1, ax2], vmin=50, vmax=500, levels=22, level_ticks=11, add_colorbar=True, 
                   cbar_pos = [0.90, 0.65, 0.02, 0.25], fig=fig, domain="South America", orientation="vertical", bottom_labels=True,
                   season_label=["[A] PI (JJA)", "[B] PI (DJF)"], title=True)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.15, hspace=0.05)
plt.savefig(os.path.join(path_to_store, "PI_tp_t2m_slp.svg"), format= "svg", bbox_inches="tight", dpi=300)

plt.show()

