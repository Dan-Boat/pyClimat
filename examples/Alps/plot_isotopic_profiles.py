#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 10:48:42 2021

@author: dboateng
This is example script for using Climat in analysing and visualising ECHAM module output (long-term seasonal differece means)
The script contains directories of module outputs and path to save plot
Note: it is structured solely for the personal needs of the author, therefore, it must be adapted advisably.
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xr
import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import scatter_plot_laspe_rate
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, extract_transect, extract_profile



# Path to experiments
module_output_main_path = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

exp_name_aw100e100 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_aw100e0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_aw100e200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_aw100e150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"

# west set up 
exp_name_aw200e100 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
exp_name_aw200e0 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
exp_name_aw200e200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"
exp_name_aw0e0 = "t015_dkrz-mistral_e5w2.3_PI_Alps0_t159l31.6h"


# for supplementary (same but for annual)
years= "1003_1017"
period = "1m"


# reading dataset
aw100e100_data, aw100e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e100, years=years,
                                                  period=period)
aw100e0_data, aw100e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e0, years=years,
                                                  period=period)
aw100e200_data, aw100e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e200, years=years,
                                                  period=period)
aw100e150_data, aw100e150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e150, years=years,
                                                  period=period)

aw200e100_data, aw200e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e100, years=years,
                                                  period=period)
aw200e0_data, aw200e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e0, years=years,
                                                  period=period)
aw200e200_data, aw200e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e200, years=years,
                                                  period=period)

aw0e0_data, aw0e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw0e0, years=years,
                                                  period=period)

#extracting variables and computing long-term means

def extract_vars_and_analysis(data, wiso, time="season", season="JJA", 
                              season_calendar="standard"):

    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    elev = extract_var(Dataset=data, varname="elev", units="m")
    
    d18op_alt = compute_lterm_mean(data=d18op, time=time, season=season,
                                   season_calendar=season_calendar)
    
    elev_alt = compute_lterm_mean(data=elev, time=time, season=season,
                                   season_calendar=season_calendar)
    
    
    
    lon_d18op = extract_profile(d18op_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon",
                                     to_pandas=True)
    lat_d18op = extract_profile(d18op_alt, maxlon=13, minlon=10, maxlat=54, minlat=40, dim="lat", 
                                     to_pandas=True)
    
    lon_elev = extract_profile(elev_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon",
                                     to_pandas=True)
    lat_elev = extract_profile(elev_alt, maxlon=13, minlon=10, maxlat=54, minlat=40, dim="lat", 
                                     to_pandas=True)
    
    return_data = {"d18Op_lon": lon_d18op, "d18Op_lat": lat_d18op,
                   "elev_lon": lon_elev, "elev_lat": lat_elev}
    
    return return_data

# plot function
def plot_profiles_all(varname, units, data_d18Op, data_elev, ax=None, path_to_store=None, filename=None,
                       colors=None, title=None, ax_legend=True,
                       ymin=None, ymax=None, dim="lon", fig=None, labels=None,
                       elev_max=None, elev_min=None, right_label=True, left_label=True,
                       bottom_label=True):
    
    
    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
    
    data_d18Op.plot(ax=ax, linestyle="--", color=colors, linewidth=3) 
    
            
    if left_label:
        ax.set_ylabel(varname + " [" + units + "]", fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_yticklabels([])
    
    if bottom_label:
        if dim == "lon":
             ax.set_xlabel("Longitude [E°]", fontsize=28, fontweight="bold")
        elif dim == "lat":
             ax.set_xlabel("Latitude [N°]", fontsize=28, fontweight="bold")
        else:
            raise ValueError("Define dim as lat or lon")
        
        
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_xticklabels([])

        
    if all(parameter is not None for parameter in [ymax, ymin]):
        ax.set_ylim(ymin, ymax)
       
    if ax_legend:
        ax.legend(frameon=True, fontsize=24, bbox_to_anchor=(0.01, 1.05, 1, 0.102,), loc="best",
                  borderaxespad=0, ncol=2)
    else:
        ax.legend([],[], frameon=False)
        
        
    ax2 = ax.twinx()
    ax2.grid(False)
    data_elev.plot(ax=ax2, linestyle="-", color=colors, linewidth=2, legend=False)
    data_elev.plot(kind="area", color=colors, alpha=0.1, legend=False, stacked=False,
                   ax=ax2) 

    if all(parameter is not None for parameter in [elev_max, elev_min]):
        ax2.set_ylim(elev_min, elev_max)
    
    if right_label == True:
    
        ax2.set_ylabel( "Elevation [m]", fontsize=24)
        
    else:
        ax2.set_yticklabels([])
        
    ax2.tick_params(axis= "y")
    ax2.tick_params(axis= "x")     
       
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="center")
        
    
        
        
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    
    if path_to_store is not None:
        plt.savefig(os.path.join(path_to_store, filename), bbox_inches="tight", format= "svg")


# extract all data 

labels_W1 = ["CTL", "W1E0", "W1E2", "W1E1.5"]
labels_W2 = ["CTL", "W2E1", "W2E0"]

data_W1 = [aw100e100_data, aw100e0_data, aw100e200_data, aw100e150_data]
wiso_W1 = [aw100e100_wiso, aw100e0_wiso, aw100e200_wiso, aw100e150_wiso]

data_W2 = [aw100e100_data, aw200e100_data, aw200e0_data]
wiso_W2 = [aw100e100_wiso, aw200e100_wiso, aw200e0_wiso]

extracts_W1 = {}
extracts_W2 = {}


def extract_all(labels, data_W1, wiso_W1, time="season", season="JJA", season_calendar="standard"):
    extracts_W1 = {}
    
    for i,topo in enumerate(labels):
        extracts_W1[topo] = extract_vars_and_analysis(data=data_W1[i], wiso=wiso_W1[i],
                                                      time=time, season=season, 
                                                      season_calendar=season_calendar)
        
        if i ==0:
            df_d18Op_W1_lon = pd.DataFrame(index=extracts_W1.get(topo)["d18Op_lon"].index.values,
                                           columns=labels_W1)
            df_d18Op_W1_lat = pd.DataFrame(index=extracts_W1.get(topo)["d18Op_lat"].index.values,
                                           columns=labels_W1)
            
            df_elev_W1_lon = pd.DataFrame(index=extracts_W1.get(topo)["elev_lon"].index.values,
                                           columns=labels_W1)
            df_elev_W1_lat = pd.DataFrame(index=extracts_W1.get(topo)["elev_lat"].index.values,
                                           columns=labels_W1)
            
            
        
        df_d18Op_W1_lon[topo] = extracts_W1.get(topo)["d18Op_lon"]
        df_d18Op_W1_lat[topo] = extracts_W1.get(topo)["d18Op_lat"]
        
        df_elev_W1_lon[topo] = extracts_W1.get(topo)["elev_lon"]
        df_elev_W1_lat[topo] = extracts_W1.get(topo)["elev_lat"]
        
    return df_d18Op_W1_lon, df_d18Op_W1_lat, df_elev_W1_lon, df_elev_W1_lat
    
df_d18Op_W1_lon, df_d18Op_W1_lat, df_elev_W1_lon, df_elev_W1_lat = extract_all(labels=labels_W1, 
                                                                               data_W1=data_W1, wiso_W1=wiso_W1) 

df_d18Op_W2_lon, df_d18Op_W2_lat, df_elev_W2_lon, df_elev_W2_lat = extract_all(labels=labels_W2, 
                                                                               data_W1=data_W2, wiso_W1=wiso_W2) 
    
    
colors_W1 = ["black", "blue", "red", "green"]  
colors_W2 = ["black", "darkgoldenrod", "purple"] 

apply_style(fontsize=28, style="seaborn-talk", linewidth=3,)

fig, ((ax1,ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(20, 15),sharey=False, sharex=False)

plot_profiles_all(varname = "$\delta^{18}$Op vs SMOW", units="‰", data_d18Op=df_d18Op_W1_lon,
                  data_elev=df_elev_W1_lon, ax=ax1, ax_legend=True,ymax=-2, ymin=-16, dim="lon",
                  fig=fig, labels=labels_W1,bottom_label=False, right_label=False, elev_max=3500,
                  elev_min=0, colors=colors_W1)

plot_profiles_all(varname = "$\delta^{18}$Op vs SMOW", units="‰", data_d18Op=df_d18Op_W1_lat,
                  data_elev=df_elev_W1_lat, ax=ax2, ax_legend=False,ymax=-2, ymin=-16, dim="lat",
                  fig=fig, labels=labels_W1, bottom_label=False, right_label=True, left_label=False,elev_max=3500,
                  elev_min=0, colors=colors_W1)

plot_profiles_all(varname = "$\delta^{18}$Op vs SMOW", units="‰", data_d18Op=df_d18Op_W2_lon,
                  data_elev=df_elev_W2_lon, ax=ax3, ax_legend=False,ymax=-2, ymin=-16, dim="lon",
                  fig=fig, labels=labels_W2, bottom_label=True, right_label=False, elev_max=3500,
                  elev_min=0, colors=colors_W2)

plot_profiles_all(varname = "$\delta^{18}$Op vs SMOW", units="‰", data_d18Op=df_d18Op_W2_lat,
                  data_elev=df_elev_W2_lat, ax=ax4, ax_legend=False,ymax=-2, ymin=-16, dim="lat",
                  fig=fig, labels=labels_W2, bottom_label=True, left_label=False, elev_max=3500,
                  elev_min=0, colors=colors_W2)
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05, wspace=0.05)
plt.savefig(os.path.join(path_to_plots, "d18Op_profile_JJA.svg"), format= "svg", bbox_inches="tight", dpi=600)


    
    
def supplementary():
    
    fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(18, 13),sharey=False, sharex=False)

    plot_iso_profiles(df_iso=aw100e100_d18op_lon , df_geosp=aw100e100_geosp_lon , dim="lon", iso_color=black, iso_label="CTL",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, title="[A]", 
                      right_labels =False, shade_color=black, shade_alpha=0.2)

    plot_iso_profiles(df_iso=aw100e200_d18op_lon , df_geosp=aw100e200_geosp_lon , dim="lon", iso_color=blue, iso_label="W1E2",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1,
                      right_labels =False,  shade_color=blue, shade_alpha=0.15)

    plot_iso_profiles(df_iso=aw200e200_d18op_lon , df_geosp=aw200e200_geosp_lon , dim="lon", iso_color=red, iso_label="W2E2",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1,
                      right_labels =False, shade_color=red, shade_alpha=0.1)
    
    ax1.grid(visible=True, linestyle="--", linewidth=1.0, color=grey)
    ax1.text(-5, -150, "W", fontsize=20)
    ax1.text(20, -150, "E", fontsize=20)
    
    plot_iso_profiles(df_iso=aw100e100_d18op_lon , df_geosp=aw100e100_geosp_lon , dim="lon", iso_color=black, iso_label=None,
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, title="[B]", 
                      left_labels =False, shade_color=black, shade_alpha=0.15)

    plot_iso_profiles(df_iso=aw100e0_d18op_lon , df_geosp=aw100e0_geosp_lon , dim="lon", iso_color=green, iso_label="W1E0",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2,
                      left_labels =False, shade_color=green, shade_alpha=0.2)

    plot_iso_profiles(df_iso=aw0e0_d18op_lon , df_geosp=aw0e0_geosp_lon , dim="lon", iso_color=purple, iso_label="W0E0",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2,
                      left_labels =False, shade_color=purple, shade_alpha=0.3)
    
    ax2.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    ax2.text(-5, -150, "W", fontsize=20)
    ax2.text(20, -150, "E", fontsize=20)
    
    
    fig.legend(frameon=True, fontsize=22, loc="upper right",)
    plt.tight_layout() 
    plt.subplots_adjust(left=0.04, right=0.86, top=0.94, bottom=0.04)
    plt.savefig(os.path.join(path_to_store, "figS8.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "figS8.png"), format= "png", bbox_inches="tight", dpi=300)
    



    
    
    

