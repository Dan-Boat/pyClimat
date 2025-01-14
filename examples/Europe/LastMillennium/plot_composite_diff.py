# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:29:02 2024

@author: dboateng
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean

from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_transect
from pyClimat.variables import extract_var
from pyClimat.utils import extract_region


from path_to_data_lm import *

main_path = "D:/Datasets/iGCM_datasets/"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots/data"


def get_indices(df, nao_mode, ea_mode, nao_factor, ea_factor):
    
    
    # extract the indices for OP and EG
    
    df_eq_index = (df[str(nao_mode)] * nao_factor > 0) & (df[str(ea_mode)] * ea_factor > 0) | (df[str(nao_mode)] * nao_factor < 0) & (df[str(ea_mode)] * ea_factor < 0)
    df_op_index = (df[str(nao_mode)] * nao_factor > 0) & (df[str(ea_mode)] * ea_factor < 0) | (df[str(nao_mode)] * nao_factor < 0) & (df[str(ea_mode)] * ea_factor > 0)
    
    # extract them from the pcs 
    
    df_eq = df[df_eq_index]
    df_op = df[df_op_index]
    
    #create xarray for the pcs
    
    EQ_indices = xr.DataArray(df_eq[str(nao_mode)] * nao_factor, dims="time", coords={"time": df_eq["time"]})
    OP_indices = xr.DataArray(df_op[str(nao_mode)] * nao_factor, dims="time", coords={"time": df_op["time"]})
    
    return EQ_indices, OP_indices


def cal_diff_OP_EQ(path_to_data, filename, pcs_path, nao_mode=1, 
                nao_factor=-1, ea_mode=3, ea_factor=1):  
    
    df = pd.read_csv(pcs_path, parse_dates=["time"])

    atmos_index_eq, atmos_index_op = get_indices(df, nao_mode, ea_mode, nao_factor, ea_factor)
    
  
    filename_data = filename + "_wiso_vars.nc"
    
    temp = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="tsurf")  - 273.15 #°C
    
    d18O = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="d18O")  
    prec = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="prec")
    
    maxlon = 150
    minlon = -140
    maxlat = 86
    minlat = 24
    
    
    t2m_season = extract_region(data=temp, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF") 
    
    prec_season = extract_region(data=prec, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF") 
    
    d18O_season = extract_region(data=d18O, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF")
    
    
    t2m_season["time"] = df.time.values
    t2m_op = t2m_season.sel(time = t2m_season.time.isin(atmos_index_op.time.values)).mean(dim="time")
    t2m_eq = t2m_season.sel(time = t2m_season.time.isin(atmos_index_eq.time.values)).mean(dim="time")
    
    
    prec_season["time"] = df.time.values
    prec_op = prec_season.sel(time = prec_season.time.isin(atmos_index_op.time.values)).mean(dim="time")
    prec_eq = prec_season.sel(time = prec_season.time.isin(atmos_index_eq.time.values)).mean(dim="time")
    
    
    d18O_season["time"] = df.time.values
    d18O_op = d18O_season.sel(time = d18O_season.time.isin(atmos_index_op.time.values)).mean(dim="time")
    d18O_eq = d18O_season.sel(time = d18O_season.time.isin(atmos_index_eq.time.values)).mean(dim="time")
    
    t2m_diff = t2m_op - t2m_eq
    prec_diff = prec_op - prec_eq
    d18O_diff = d18O_op - d18O_eq
    
    result = xr.Dataset({"t2m_diff":t2m_diff.sortby("lon"), 
                         "prec_diff":prec_diff.sortby("lon"), 
                         "d18O_diff":d18O_diff.sortby("lon")})
    
    return result


cesm_data = cal_diff_OP_EQ(path_to_data=main_path, filename="CESM", pcs_path=CESM_DJF_PCs,
                           nao_mode=1, nao_factor=-1, ea_mode=3, ea_factor=1)

giss_data = cal_diff_OP_EQ(path_to_data=main_path, filename="GISS", pcs_path=GISS_DJF_PCs,
                           nao_mode=1, nao_factor=-1, ea_mode=3, ea_factor=1)



apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)

fig,((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(22,25), subplot_kw={"projection":projection})

plot_annual_mean(variable="Precipitation", data_alt=cesm_data["prec_diff"], ax=ax1,
                 cmap=BrBG, units="mm/month", vmax=10, vmin=-10, plot_coastlines=True, 
                 levels=22, level_ticks=7, add_colorbar=True, cbar_pos= [0.95, 0.65, 0.02, 0.20], 
                 orientation="vertical", bottom_labels=False,left_labels=False, fig=fig, plot_borders=True,
                 plot_projection=projection, title="(a) CESM  (OP-EQ)", center=True, domain="NH",
                 label_format="%.0f", )

plot_annual_mean(variable="Precipitation", data_alt=giss_data["prec_diff"], ax=ax2,
                 cmap=BrBG, units="mm/month", vmax=10, vmin=-10,
                 plot_coastlines=True, levels=22, level_ticks=7, add_colorbar=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=True,
                plot_projection=projection, title="(b) GISS", center=True, domain="NH",
                label_format="%.0f")

plot_annual_mean(variable="Temperature", data_alt=cesm_data["t2m_diff"], ax=ax3,
                 cmap=RdBu_r, units="°C", vmax=1.5, vmin=-1.5, 
                levels=22, level_ticks=7, add_colorbar=True, cbar_pos= [0.95, 0.35, 0.02, 0.20], 
                orientation="vertical", bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=True,
                plot_projection=projection, title="(c)", center=True, domain="NH",
                label_format="%.1f")

plot_annual_mean(variable="Temperature", data_alt=giss_data["t2m_diff"], ax=ax4,
                 cmap=RdBu_r, units="°C", vmax=1.5, vmin=-1.5, 
                levels=22, level_ticks=7, add_colorbar=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=True,
                plot_projection=projection, title="(d)", center=True, domain="NH",
                label_format="%.1f")

plot_annual_mean(variable="$\delta^{18}$Op VSMOW", data_alt=cesm_data["d18O_diff"], ax=ax5,
                 cmap="PiYG", units="‰", vmax=0.5, vmin=-0.5, 
                levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.95, 0.05, 0.02, 0.20], 
                orientation="vertical", bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=True,
                plot_projection=projection, title="(e)", center=True, domain="NH",
                label_format="%.1f")

plot_annual_mean(variable="$\delta^{18}$Op VSMOW", data_alt=giss_data["d18O_diff"], ax=ax6,
                 cmap="PiYG", units="‰", vmax=0.5, vmin=-0.5, 
                levels=22, level_ticks=11, add_colorbar=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=True,
                plot_projection=projection, title="(e)", center=True, domain="NH",
                label_format="%.1f")
    
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
plt.savefig(os.path.join(path_to_plots, "NAO_op_eq_diff.png"), format= "png", bbox_inches="tight", dpi=600)



