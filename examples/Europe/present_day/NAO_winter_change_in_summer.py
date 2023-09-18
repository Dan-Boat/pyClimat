# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 15:15:25 2023

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
from pyClimat.plots import plot_eofsAsCovariance
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_transect
from pyClimat.variables import extract_var
from pyClimat.utils import extract_region


from paths_to_data import *
echam_pd_data_path = "D:/Datasets/Model_output_pst/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021"

df_echam_pcs_DJF = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

df_echam_pcs_JJA = pd.read_csv(PCs_ECHAM_JJA_path, parse_dates=["time"])


df_djf = df_echam_pcs_DJF.rename(columns={"1":"NAO_DJF", "time": "time_winter"})
df_djf = df_djf[["time_winter", "NAO_DJF"]]
df_djf["NAO_DJF"] = df_djf["NAO_DJF"] * -1

df_jja = df_echam_pcs_JJA.rename(columns={"1":"NAO_JJA", "time": "time_summer"})
df_jja = df_jja[["time_summer", "NAO_JJA"]]
df_jja["NAO_JJA"] = df_jja["NAO_JJA"] 


df = df_djf.join(df_jja)

time_summer_nao_positive = df[df["NAO_DJF"] > 0]
time_summer_nao_negative = df[df["NAO_DJF"] < 0]

nao_positive = xr.DataArray(time_summer_nao_positive["NAO_DJF"], dims="time",
                            coords={"time": time_summer_nao_positive["time_summer"]})

nao_negative = xr.DataArray(time_summer_nao_negative["NAO_DJF"], dims="time",
                            coords={"time": time_summer_nao_negative["time_summer"]})


def cal_diff_nao_DJF_index_in_JJA(varname, index_positive, index_negative,
                          lev=None, lev_units=None, season="JJA", units=None):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data = extract_var(Dataset=echam_data, varname=varname, units=units, Dataset_wiso=echam_wiso,
                       lev=lev, lev_units=lev_units)
    
    data_season = extract_region(data=data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season=season)
    
    data_season["time"] = data_season.time.values.astype("datetime64[ns]")
    
    data_positive = data_season.sel(time = data_season.time.isin(index_positive.time.values)).mean(dim="time")
    data_negative = data_season.sel(time = data_season.time.isin(index_negative.time.values)).mean(dim="time")
    
    data_diff = data_positive - data_negative
    data_diff = data_diff.sortby("lon")
    
    return data_diff


d18Op_diff = cal_diff_nao_DJF_index_in_JJA(varname="d18op", index_positive=nao_positive, 
                                           index_negative=nao_negative, units="per mil")

temp_diff = cal_diff_nao_DJF_index_in_JJA(varname="temp2", index_positive=nao_positive, 
                                           index_negative=nao_negative, units="°C")

prec_diff = cal_diff_nao_DJF_index_in_JJA(varname="prec", index_positive=nao_positive, 
                                           index_negative=nao_negative, units="mm/month")

from pyClimat.plots import plot_annual_mean

apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
#apply_style(fontsize=22, style=None, linewidth=2)
projection = ccrs.PlateCarree()
fig, ((ax1, ax2, ax3)) = plt.subplots(nrows = 1, ncols=3, 
                                                 figsize=(28, 13), subplot_kw={"projection": projection})


plot_annual_mean(variable='Precipitation', data_alt=prec_diff, cmap=BrBG, units="mm/month",
                 ax=ax3, fig=fig, vmax=30, vmin=-30, levels=22, domain="Europe", level_ticks=9, 
                 cbar_pos = [0.65, 0.01, 0.30, 0.02], title=" Prec (NAO+ - NAO-(DJF)) [JJA]",
                 left_labels=True, bottom_labels=True)

plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=d18Op_diff, cmap="PiYG", units="‰",
                 ax=ax1, fig=fig, vmax=2, vmin=-2, levels=22, domain="Europe", level_ticks=11, 
                 cbar_pos = [0.05, 0.01, 0.30, 0.02], title="$\delta^{18}$Op (NAO+ - NAO-(DJF)) [JJA]",
                 left_labels=True, bottom_labels=True, label_format="%.1f")

plot_annual_mean(variable='Temperature', data_alt=temp_diff, cmap=RdBu_r, units="°C",
                 ax=ax2, fig=fig, vmax=2, vmin=-2, levels=22, domain="Europe", level_ticks=5, 
                 cbar_pos = [0.35, 0.01, 0.30, 0.02], title="Temp (NAO+ - NAO-(DJF)) [JJA]",
                 left_labels=True, bottom_labels=True, label_format="%.1f")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.10, hspace=0.001)
plt.savefig(os.path.join(path_to_plots, "composite_nao_djf_climete_jja.pdf"), format= "pdf", bbox_inches="tight", dpi=300)
plt.show()