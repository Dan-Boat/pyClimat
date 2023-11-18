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
from pyClimat.data import read_ERA_processed, read_ECHAM_processed, read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff, extract_var, extract_transect

#useful functions
#define variables
t2m = "temp2"
prec = "prec"
v10 = "v10"
u10 = "u10"


def extract_only_prec_temp_winds(data_surface):
    
    data_t2m = extract_var(Dataset=data_surface, varname=t2m, units="°C")
    
    data_prec = extract_var(Dataset=data_surface, varname=prec, units="mm/month")
    
    # data_v10 = extract_var(Dataset=data_surface, varname=v10) # default in m/s
    # data_u10 = extract_var(Dataset=data_surface, varname=u10) # default in m/s
    
    
    return data_prec, data_t2m

def weighted_mean(data):
    weights = np.cos(np.deg2rad(data.lat))
    weights.name = "weights"
    
    data_weighted = data.weighted(weights)
    
    data_mean = data_weighted.mean(dim=("lat", "lon"), skipna=True)
    return data_mean 



def compute_mean(data, save=False, path_to_save=None):
    data_mean = weighted_mean(data)
    
    
    import calendar
    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    df = pd.DataFrame(index = mnames, columns = ["mean"])
    df["mean"] = data_mean
    
    if save == True:
        df.to_csv(path_to_save)
    return df


def extract_sections(data):
    
    max_lon, max_lat, min_lon, min_lat = 30, 20, -20, 10
    
    extract_sahel = extract_transect(data=data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat)
    
    df_sahel = compute_mean(extract_sahel)
    
    max_lon, max_lat, min_lon, min_lat = 30, 10, -20, 5
    
    extract_guinea = extract_transect(data=data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat)
    
    df_guinea = compute_mean(extract_guinea)
    
    max_lon, max_lat, min_lon, min_lat = 30, 30, -20, 20
    
    extract_sahara = extract_transect(data=data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat)
    
    df_sahara = compute_mean(extract_sahara)
    
    return df_sahara, df_sahel, df_guinea


def plot_monthly_period_per_section(PD, ERA, ax, title, ymax, ymin, varname=None):
    
    ax.plot(PD["mean"], "-", color=black, label="ECHAM5", linewidth=4)
    ax.plot(ERA["mean"], "-", color=red, label="ERA5", linewidth=4)
    
    
    if varname is not None:
        ax.set_ylabel(varname, fontweight="bold", fontsize=22)
    
    #ax.grid(False, linestyle="--", color=grey, alpha=0.8)
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc= "left")
    ax.axes.tick_params(which="both", labelsize=25)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
 
    

    
    
# define paths 
main_path = "D:/Datasets/Model_output_pst/"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"
ERA5_path =  "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021/"

ERA5_tp_path = os.path.join(ERA5_path, "tp_monthly.nc")
ERA5_t2m_path = os.path.join(ERA5_path, "t2m_monthly.nc")
# ERA5_v10_path = os.path.join(ERA5_path, "v10_monthly.nc")
# ERA5_u10_path = os.path.join(ERA5_path, "u10_monthly.nc")

#define time of amip 
from1980to2014 = pd.date_range(start="1980-01-01", end="2014-12-31", freq="MS")


#read in postprocessed and analysed data 
ERA5_t2m = read_ERA_processed(path=ERA5_t2m_path, varname="t2m")   - 273.15 #°C
ERA5_tp = read_ERA_processed(path=ERA5_tp_path, varname="tp") * 1000 * 30  #mm/month
# ERA5_v10 = read_ERA_processed(path=ERA5_v10_path, varname="v10") #m/s
# ERA5_u10 = read_ERA_processed(path=ERA5_u10_path, varname="u10") #m/s

pd_path = os.path.join(main_path, "PD")
PD_data = read_from_path(pd_path, "PD_1980_2014_monthly.nc", decode=True)


PD_prec, PD_t2m= extract_only_prec_temp_winds(PD_data)

#ERA
ERA5_t2m_alt = compute_lterm_mean(data=ERA5_t2m, time="month", month="JJAS", time_range=from1980to2014)
ERA5_tp_alt = compute_lterm_mean(data=ERA5_tp, time="month", month="JJAS", time_range=from1980to2014)


ERA5_t2m_mon = compute_lterm_mean(data=ERA5_t2m, time="month", time_range=from1980to2014)
ERA5_tp_mon = compute_lterm_mean(data=ERA5_tp, time="month", time_range=from1980to2014)
# ERA5_v10_alt = compute_lterm_mean(data=ERA5_v10, time="month", month="JJAS", time_range=from1980to2014)
# ERA5_u10_alt = compute_lterm_mean(data=ERA5_u10, time="month", month="JJAS", time_range=from1980to2014)


#PD
PD_t2m_alt = compute_lterm_mean(data=PD_t2m, time="month", month="JJAS")
PD_prec_alt = compute_lterm_mean(data=PD_prec, time="month", month="JJAS")

PD_t2m_mon = compute_lterm_mean(data=PD_t2m, time="month")
PD_prec_mon = compute_lterm_mean(data=PD_prec, time="month")


# PD_v10_alt = compute_lterm_mean(data=PD_v10, time="month", month="JJAS")
# PD_u10_alt = compute_lterm_mean(data=PD_u10, time="month", month="JJAS")

#ERA5-ECHAM5
# interpolate ERA to ECHAM resolution
ERA5_t2m_alt = ERA5_t2m_alt.rename({"longitude": "lon", "latitude":"lat"})
ERA5_tp_alt = ERA5_tp_alt.rename({"longitude": "lon", "latitude":"lat"})


ERA5_t2m_alt_interp = ERA5_t2m_alt.interp(lat=PD_t2m_alt.lat).interp(lon=PD_t2m_alt.lon)
ERA5_tp_alt_interp = ERA5_tp_alt.interp(lat=PD_prec_alt.lat).interp(lon=PD_prec_alt.lon)

t2m_diff = PD_t2m_alt - ERA5_t2m_alt_interp
prec_diff = PD_prec_alt - ERA5_tp_alt_interp

#extract the monthly for the different sections
pd_sahara, pd_sahel, pd_guinea = extract_sections(data=PD_prec_mon)
era_sahara, era_sahel, era_guinea = extract_sections(data=ERA5_tp_mon)

# set up plot function for specific variable

apply_style(fontsize=22, style=None, linewidth=2) 

def plot_monthly_variability():
    apply_style(fontsize=25, style=None, linewidth=3) 
    fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(30, 12), sharey=False)
    
    plot_monthly_period_per_section(PD=pd_sahara , ERA= era_sahara, ax=ax1,
                                    title="Sahara", varname="Precipitation [mm/month]", ymax=5,
                                    ymin=0, )
    
    plot_monthly_period_per_section(PD=pd_sahel , ERA= era_sahel, 
                                    ax=ax2, title="Sahel", ymax=120,
                                    ymin=0, varname="Precipitation [mm/month]")
    
    plot_monthly_period_per_section(PD=pd_guinea , ERA= era_guinea,
                                   ax=ax3, title="Coast of Guinea", ymax=280,
                                    ymin=0, varname="Precipitation [mm/month]")
    
    
    ax2.legend(bbox_to_anchor=(0.01, 1.04, 1., 0.102), loc=3, ncol=4, borderaxespad=0., frameon = True, 
                  fontsize=20)
    
    axes = [ax1,ax2, ax3]
    for ax in axes:
        ax.grid(True, linestyle="--", color="grey")
    plt.tight_layout()
    plt.savefig(os.path.join(path_to_plots, "monthly_sections_sahara_sahel_coast_era_echam.pdf"), format= "pdf", bbox_inches="tight", dpi=300)


def plot_compare():
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
    
    
plot_monthly_variability()
#plot_compare()