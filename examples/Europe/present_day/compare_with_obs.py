# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 20:10:19 2023

@author: dboateng

This script aim to compare the NAO and EA from ERA5 and ECHAM to obs records
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

#defining paths
main_path_to_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/PD"
main_path_to_obs = "C:/Users/dboateng/Desktop/Datasets/NAO"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"


EA_valencia_DJF_path = os.path.join(main_path_to_obs, "EA_Valencia_DJF.csv")
EA_valencia_JJA_path = os.path.join(main_path_to_obs, "EA_Valencia_JJA.csv")

NAO_Gilbraltar_JJA_path = os.path.join(main_path_to_obs, "NAO_Gilbraltar_JJA.csv")
NAO_Gilbraltar_DJF_path = os.path.join(main_path_to_obs, "NAO_Gilbraltar_DJF.csv")

NAO_CDC_JJA_path = os.path.join(main_path_to_obs, "NAO_CDC_JJA.csv")
NAO_CDC_DJF_path = os.path.join(main_path_to_obs, "NAO_CDC_DJF.csv")

EOFs_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_eofsAsCovariance.nc") # standard method
EOFs_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_eofsAsCovariance.nc")

PCs_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_pcs.csv") # standard method
PCs_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_pcs.csv")

Vars_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_pcs_variance.csv") # standard method
Vars_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_pcs_variance.csv")


EOFs_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_eofsAsCovariance.nc") # standard method
EOFs_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_eofsAsCovariance.nc")

PCs_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_pcs.csv") # standard method
PCs_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_pcs.csv")

Vars_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_pcs_variance.csv") # standard method
Vars_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_pcs_variance.csv")



# read all the required datasets (start with winter)

# read EOFs and plot with thier variance
def plot_DJF_NAO_Covariance():
    eof_era_djf_data = xr.open_dataarray(EOFs_ERA_DJF_path)
    eof_echam_djf_data = xr.open_dataarray(EOFs_ECHAM_DJF_path)
    variance_era_djf = pd.read_csv(Vars_ERA_DJF_path, index_col=["mode"])
    variance_echam_djf = pd.read_csv(Vars_ECHAM_DJF_path, index_col=["mode"])
    
    # plot eofs
    
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.AlbersEqualArea(
        central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))
    
    fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols=2, 
                                                 figsize=(18, 14), subplot_kw={"projection": projection})
    
    modes = [1,2]
    units="hPa" 
    variable="Mean Sea Level Pressure"
    vmax=15
    vmin=-15
    figname ="NAO_era_echam_obs_DJF"
    plot_eofsAsCovariance(variable= variable, data=eof_era_djf_data.sel(mode=modes[0]) * -1, mode_var=variance_era_djf.loc[modes[0]], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                          ax=ax1, fig=fig, title="[A] ERA5 (1959-2021)", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=eof_echam_djf_data.sel(mode=modes[0]) * -1, mode_var=variance_echam_djf.loc[modes[0]], units=units, 
                          vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] ECHAM5-wiso (1979-2014)", bottom_labels=True, use_AlberEqualArea=True)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.88, top=0.90, bottom=0.06, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, figname + ".svg"), format= "svg", bbox_inches="tight", dpi=300)


def plot_DJF_EA_Covariance():
    eof_era_djf_data = xr.open_dataarray(EOFs_ERA_DJF_path)
    eof_echam_djf_data = xr.open_dataarray(EOFs_ECHAM_DJF_path)
    variance_era_djf = pd.read_csv(Vars_ERA_DJF_path, index_col=["mode"])
    variance_echam_djf = pd.read_csv(Vars_ECHAM_DJF_path, index_col=["mode"])
    
    # plot eofs
    
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.AlbersEqualArea(
        central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))
    
    fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols=2, 
                                                 figsize=(18, 14), subplot_kw={"projection": projection})
    
    modes = [1,2,3,4]
    units="hPa" 
    variable="Mean Sea Level Pressure"
    vmax=15
    vmin=-15
    figname ="EA_era_echam_obs_DJF"
    plot_eofsAsCovariance(variable= variable, data=eof_era_djf_data.sel(mode=modes[2]), mode_var=variance_era_djf.loc[modes[2]], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                          ax=ax1, fig=fig, title="[A] ERA5 (1959-2021)", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=eof_echam_djf_data.sel(mode=modes[1]), mode_var=variance_echam_djf.loc[modes[1]], units=units, 
                          vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] ECHAM5-wiso (1979-2014)", bottom_labels=True, use_AlberEqualArea=True)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.88, top=0.90, bottom=0.06, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, figname + ".svg"), format= "svg", bbox_inches="tight", dpi=300)

#plot_DJF_EA_Covariance()


#plot the NAO indices for all and vs obs
# read all the indices
pcs_era_djf = pd.read_csv(PCs_ERA_DJF_path, index_col=["time"])
pcs_echam_djf = pd.read_csv(PCs_ECHAM_DJF_path, index_col=["time"])
nao_gilbraltar_djf = pd.read_csv(NAO_Gilbraltar_DJF_path, index_col=["time"])
nao_cdc_djf = pd.read_csv(NAO_CDC_DJF_path, index_col=["time"])
ea_valencia_djf = pd.read_csv(EA_valencia_DJF_path, index_col=["time"])

#change index to pd.datetime
nao_gilbraltar_djf.index = pd.to_datetime(nao_gilbraltar_djf.index)
nao_gilbraltar_djf = nao_gilbraltar_djf.loc["1950-01-12":"2021-01-02"]

ea_valencia_djf.index = pd.to_datetime(ea_valencia_djf.index)
ea_valencia_djf = ea_valencia_djf.loc["1950-01-12":"2021-01-02"]

nao_cdc_djf.index = pd.to_datetime(nao_cdc_djf.index)
pcs_era_djf.index = pd.to_datetime(pcs_era_djf.index)
pcs_echam_djf.index = pd.to_datetime(pcs_echam_djf.index)

# extract the indices (NAO)
df_nao_gib = nao_gilbraltar_djf["NAO"] / nao_gilbraltar_djf["NAO"].std()
df_nao_era = pcs_era_djf["1"].loc["1950-12-01":"2021-02-01"]                                          
df_nao_echam = pcs_echam_djf["1"].loc["1980-12-01":"2014-02-01"]
df_nao_cdc = nao_cdc_djf["NAO"].loc["1950-01-12":"2021-01-02"]

# extract indices EA

#calculate the winter mean and anchor it to the december of the previous year
#df = df_nao_era.groupby(df_nao_era.reset_index().index //3)

df_ea_val = ea_valencia_djf["EA"]
df_ea_era = pcs_era_djf["2"].loc["1950-12-01":"2021-02-01"]                                          
df_ea_echam = pcs_echam_djf["2"].loc["1980-12-01":"2014-02-01"]



winter_year = pd.date_range(start='1950-12-01', end="2020-12-01", freq='12MS')
winter_year_val = pd.date_range(start='1950-12-01', end="2016-12-01", freq='12MS')
winter_year_echam = pd.date_range(start='1981-12-01', end="2015-12-01", freq='12MS')

def calculate_season_avg(df, time_range):
    df = df.groupby(df.reset_index().index // 3).mean()
    df_avg = pd.DataFrame(index=time_range, columns=["season_avg"])
    df_avg["season_avg"] = df.values
    
    return df_avg



df_nao_gib_avg = calculate_season_avg(df=df_nao_gib, time_range=winter_year)
df_nao_era_avg = calculate_season_avg(df=df_nao_era, time_range=winter_year)
df_nao_cdc_avg = calculate_season_avg(df=df_nao_cdc, time_range=winter_year)
df_nao_echam_avg = calculate_season_avg(df=df_nao_echam, time_range=winter_year_echam)

df_ea_val_avg = calculate_season_avg(df=df_ea_val, time_range=winter_year_val)
df_ea_era_avg = calculate_season_avg(df=df_ea_era, time_range=winter_year)
df_ea_echam_avg = calculate_season_avg(df=df_ea_echam, time_range=winter_year_echam)


fig_name ="EA_index_DJF.svg"

apply_style(fontsize=22, style=None, linewidth=2)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize= (18, 10), sharex=False)
ax.plot(df_ea_era_avg, color="black", linewidth=2, label="ERA5")
ax.plot(df_ea_echam_avg*-1, color="blue", linewidth=2, label="ECHAM5-wiso", linestyle="--")
ax.plot(df_ea_val_avg, color="red", linestyle="-", linewidth=2, label="Valencia")
ax.xaxis.set_major_locator(YearLocator(5))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.axhline(y=0, linestyle="--", color=grey, linewidth=2)
ax.legend(bbox_to_anchor=(0.01, 1.02, 1., 0.102), loc=3, ncol=4, borderaxespad=0., frameon = True, 
          fontsize=20)

ax.set_ylabel("EA index", fontweight="bold", fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(path_to_plots, fig_name), bbox_inches="tight", format= "svg")
plt.show()



fig_name ="NAO_index_DJF.svg"

apply_style(fontsize=22, style=None, linewidth=2)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize= (18, 10), sharex=False)
ax.plot(df_nao_era_avg*-1, color="black", linewidth=2, label="ERA5")
ax.plot(df_nao_echam_avg *-1, color="blue", linewidth=2, label="ECHAM5-wiso", linestyle="--")
ax.plot(df_nao_gib_avg, color="red", linestyle="-", linewidth=2, label="Gibraltar")
ax.plot(df_nao_cdc_avg, color="green", linestyle="-", linewidth=2, label="CDC (NOAA)")
ax.xaxis.set_major_locator(YearLocator(5))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.axhline(y=0, linestyle="--", color=grey, linewidth=2)
ax.legend(bbox_to_anchor=(0.01, 1.02, 1., 0.102), loc=3, ncol=4, borderaxespad=0., frameon = True, 
          fontsize=20)

ax.set_ylabel("NAO index", fontweight="bold", fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(path_to_plots, fig_name), bbox_inches="tight", format= "svg")
plt.show()

