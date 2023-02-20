# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:46:22 2023

@author: dboateng

Note that the method for the EOF can be the standard one from Eof.xarray, or the xefos methods (sklearn),
    or the varimax and promax extentions

Perform all the require analysis (eg. EOFs, sorting of indixes, correlation, regional means)
1. Try the routine for extracting pcs and projecting X fields onto the PD eofs
2. Try the correlation scheme for both one time axis with grids and also two 3D arrays for spearman and pearson, 
get the pvalues and also the uncertainties

"""

import os 
import numpy as np 
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt 
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.analysis import extract_var
from pyClimat.stats import EOF_standard
from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance

from pyClimat.data import read_ECHAM_processed, read_from_path, read_ERA_processed

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"
main_path = "D:/Datasets/Model_output_pst/"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021/"

lgm_path = os.path.join(main_path, "LGM")
plio_path = os.path.join(main_path, "PLIO")
mh_path = os.path.join(main_path, "MH")
pi_path = os.path.join(main_path, "PI")
pd_path = os.path.join(main_path, "PD")



# read data
PD_data = read_from_path(pd_path, "PD_1980_2014_monthly.nc", decode=True)
PI_data = read_from_path(pi_path, "PI_1003_1017_monthly.nc", decode=True)
LGM_data = read_from_path(lgm_path,"LGM_1003_1017_monthly.nc", decode=True)
PLIO_data = read_from_path(plio_path, "PLIO_1003_1017_monthly.nc", decode=True)
MH_data = read_from_path(mh_path, "MH_1003_1017_monthly.nc", decode=True)



ERA5_tp_path = os.path.join(ERA5_path, "tp_monthly.nc")
ERA5_t2m_path = os.path.join(ERA5_path, "t2m_monthly.nc")
ERA5_msl_path = os.path.join(ERA5_path, "msl_monthly.nc")
ERA5_z500_path = os.path.join(ERA5_path, "z500_monthly.nc")





def extract_eofs_data(data, figname, units, variable, vmax=15, vmin=-15, plot_covariance=True, is_era=False,
                      path_to_plots=None):
    # analysis for ERA5 Dataset
    
    
    
    # initiate the eof instance 
    EOF = EOF_standard(data=data, weights=True, standardize=True, 
                          extract_region=True, extract_season=True, neofs=4)
    
    # select the region of interest and season
    EOF.select_time_and_region(maxlon=60, minlon=-80, maxlat=80, minlat=20, time="season", 
                                  season="DJF", month="JA") # month="AMJJAS"
    
    # calculate the anomalies and apply norm
    EOF.calculate_anomalies(monthly_anomalies=True)
    
    method = "xeofs"
    
    # fit the eof with the solver
    EOF.eof_solver(method=method, apply_promax=False, apply_varimax=True)
    
    # extract the eofs (as the eigenvectors of the covariance matric of the X field)
    
    eofs = EOF.eofs(eofscaling=2) # 2 - multiply by the singular values or square root of the eigen values
    pcs = EOF.pcs(pscaling=0)
    variance = EOF.explained_variance_ratio()
    
    if plot_covariance:
    # loop through this !!
        plot_eofs(data=eofs, variance=variance, figname=figname + ".svg", units=units, variable=variable, vmax=vmax, 
                  vmin=vmin, is_era=is_era, path_to_plots=path_to_plots)
    
    
    return pcs
    
    
    
def plot_eofs(data, variance, figname, units="m", variable="slp", vmax=15, vmin=-15, is_era=False,
              path_to_plots=None):
    
    apply_style(fontsize=22, style=None, linewidth=2) 
    #projection = ccrs.PlateCarree()
    
    projection = ccrs.AlbersEqualArea(
        central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols=2, 
                                                 figsize=(24, 20), subplot_kw={"projection": projection})
    plot_eofsAsCovariance(variable= variable, data=data.sel(mode=1), mode_var=variance[1], units=units, vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                          ax=ax1, fig=fig, title="[A]", bottom_labels=False)
    
    plot_eofsAsCovariance(variable= variable, data=data.sel(mode=2), mode_var=variance[2], units=units, vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B]", bottom_labels=False, use_AlberEqualArea=True)
    
    if is_era:
        
        plot_eofsAsCovariance(variable= variable, data=data.sel(mode=3) *-1, mode_var=variance[3], units=units, vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                              level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C]", bottom_labels=True, use_AlberEqualArea=True)
        
        plot_eofsAsCovariance(variable= variable, data=data.sel(mode=4) *-1, mode_var=variance[4], units=units, vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                              level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D]", bottom_labels=True, use_AlberEqualArea=True)
        
    else:
        plot_eofsAsCovariance(variable= variable, data=data.sel(mode=3), mode_var=variance[3], units=units, vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                              level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C]", bottom_labels=True, use_AlberEqualArea=True)
        
        plot_eofsAsCovariance(variable= variable, data=data.sel(mode=4), mode_var=variance[4], units=units, vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                              level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D]", bottom_labels=True, use_AlberEqualArea=True)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, figname + ".svg"), format= "svg", bbox_inches="tight", dpi=300)
    
    
    
#ERA5

#def read_ERA(): 
#ERA5_t2m = read_ERA_processed(path=ERA5_t2m_path, varname="t2m")   - 273.15 #°C
#ERA5_tp = read_ERA_processed(path=ERA5_tp_path, varname="tp") * 1000 * 30  #mm/month

from1979to2014 = pd.date_range(start="1979-01-01", end="2014-12-31", freq="MS")
ERA5_msl = read_ERA_processed(path=ERA5_msl_path, varname="msl") / 100 #Pa --> hPa
ERA5_msl = ERA5_msl.sel(time=from1979to2014)


# ERA5_z500 = read_ERA_processed(path=ERA5_z500_path, varname="z500") * 9.81 #--> m
# ERA5_z500 = ERA5_z500.sel(level=500)  # to drop the level dim since it is indexed on it

#extract_eofs_data(data=ERA5_z500, figname="ERA_z500_varimax_wNAO", units="m", variable="Geopotential height")
ERA_pcs = extract_eofs_data(data=ERA5_msl, figname="ERA_msl_varimax_wNAO", units="hPa", variable="Mean Sea Level Pressure", vmax=10, vmin=-10, is_era=True)


# select only the years of echam
#ERA_pcs = ERA_pcs["1980-01-01": "2014-12-01"]


# analysis for echam PD
PD_msl = extract_var(Dataset=PD_data, varname="slp", units="hPa")
PI_msl = extract_var(Dataset=PI_data, varname="slp", units="hPa")
LGM_msl = extract_var(Dataset=LGM_data, varname="slp", units="hPa")
MH_msl = extract_var(Dataset=MH_data, varname="slp", units="hPa")
PLIO_msl = extract_var(Dataset=PLIO_data, varname="slp", units="hPa")
#PD_h500 = extract_var(Dataset=PD_data, varname="geopoth", lev_units="hPa", lev=500)

#extract_eofs_data(data=PD_h500, figname="PD_z500_varimax_wNAO", units="m", variable="Geopotential height")
PD_pcs = extract_eofs_data(data=PD_msl, figname="PD_msl_varimax_wNAO", units="hPa", variable="Mean Sea Level Pressure", vmax=10, vmin=-10)
PI_pcs = extract_eofs_data(data=PI_msl, figname="PI_msl_varimax_wNAO", units="hPa", variable="Mean Sea Level Pressure", vmax=10, vmin=-10)

LGM_pcs = extract_eofs_data(data=LGM_msl, figname="LGM_msl_varimax_wNAO", units="hPa", variable="Mean Sea Level Pressure", vmax=10, vmin=-10)
MH_pcs = extract_eofs_data(data=MH_msl, figname="MH_msl_varimax_wNAO", units="hPa", variable="Mean Sea Level Pressure", vmax=10, vmin=-10)
PLIO_pcs = extract_eofs_data(data=PLIO_msl, figname="PLIO_msl_varimax_wNAO", units="hPa", variable="Mean Sea Level Pressure", vmax=10, vmin=-10)


PD_pcs.index = PD_pcs.index.to_datetimeindex()
# visualize the pcs
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize= (20, 15), sharex=True)

line1 = ax1.plot(ERA_pcs[1], "--", color="black")
line2 = ax2.plot(PD_pcs[1], "--", color="black")



lines = [line1, line2]
axes = [ax1, ax2]

for i,line in enumerate(lines):
    x, y = line[0].get_data()
    
    # Fill above and below the zero line
    axes[i].fill_between(x, y, where=y > 0, color='red', interpolate=True)
    axes[i].fill_between(x, y, where=y < 0, color='blue', interpolate=True)
    
    
    axes[i].xaxis.set_major_locator(YearLocator(5))
    axes[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    axes[i].axhline(y=0, linestyle="--", color=grey, linewidth=2)

# ax1.plot(ERA_pcs[1].rolling(12,min_periods=1, win_type="hann", center=True).mean(), color="black", linewidth=2,)
# ax2.plot(PD_pcs[1].rolling(12,min_periods=1, win_type="hann", center=True).mean(), color="black", linewidth=2)

plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, "era_echam_pcs" + ".svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.show()




