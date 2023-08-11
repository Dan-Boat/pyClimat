#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 15:23:06 2021

@author: dboateng
This is example script for using Climat in analysing and visualising ECHAM module output (isotopic lapse rate)
The script contains directories of module outputs and path to save plot
Note: it is structured solely for the personal needs of the author, therefore, it must be adapted advisably.
"""

#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs
import seaborn as sns


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import scatter_plot_laspe_rate
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, extract_transect, linregression
from pyClimat.variables import extract_var


module_output_main_path = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/plots"

# Path to experiments
exp_name_aw100e100 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_aw100e0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_aw100e200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_aw100e150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"

# west set up 
exp_name_aw200e100 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
exp_name_aw200e0 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
exp_name_aw200e200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"

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

#extracting variables and computing long-term means

def extract_vars_and_analysis(data, wiso, maxlon=5, maxlat=46, minlon=-1, minlat=43,
                              time="season", season="JJA", season_calendar="standard"):

    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    elev = extract_var(Dataset=data , varname="elev", units="m")

    # compute annual means
    d18op_alt = compute_lterm_mean(data=d18op, time=time, season=season, season_calendar=season_calendar)
    elev_alt = compute_lterm_mean(data=elev, time=time, season=season, season_calendar=season_calendar)
    
    elev_transect = extract_transect(data=elev_alt, maxlon=maxlon, minlon=minlon, 
                                  maxlat=maxlat, minlat=minlat, sea_land_mask=True, 
                                  Dataset=data)
    
    
    
    d18op_transect = extract_transect(data=d18op_alt, maxlon=maxlon, minlon=minlon, 
                                  maxlat=maxlat, minlat=minlat, sea_land_mask=True, 
                                  Dataset=data)
   
   
    reg, df, pred = linregression(data_x=elev_transect, data_y=d18op_transect, return_yhat=True, get_ci=True)
    
    
    return_data = {"reg":reg, "df":df, "pred":pred}
    
    
    return return_data



# plot function 
def extract_label(data):
    
    # labels
    slope_pi = data.get("reg").slope*1000
    coef_err_pi = data.get("reg").stderr*1000
    r2_pi = data.get("reg").rvalue*-1
    
    label = "{:.2f} +/-{:.2f} [‰/km], r²={:.2f} (".format(slope_pi, coef_err_pi,r2_pi)
    
    return label




def plot_section_lapse_rate(dataset, makers, colors, labels, ax=None, xlabel=True, ylabel=True, title=None, 
                            ax_legend=True, ymin=None, ymax=None, xmax=None, xmin=None):
    

    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
    
                                                              
    # plot the regression line and points
    
    for i,topo in enumerate(labels):
        data = dataset.get(topo)
        
        label_data = extract_label(data)
        sns.regplot(data=data.get("df"), x="X", y="Y", marker=makers[i], 
                   scatter_kws={"color":colors[i], "s":200, "alpha":0.9}, color=colors[i],ax=ax,
                   label=label_data + topo + ")")
        
        # plot the prediction interval
        
        ax.plot(data.get("df")["X"], data.get("pred")['obs_ci_lower'], linestyle="-",
                linewidth=0.3, color=colors[i],)
        
        ax.plot(data.get("df")["X"], data.get("pred")['obs_ci_upper'], linestyle="-",
                linewidth=0.3, color=colors[i],)
    
    
    # add the labels (with pvalues estimates)
            
            
    if ylabel:
        ax.set_ylabel(u'$\delta^{18}$O ‰ vs SMOW (‰)', fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_yticklabels([])
    
    if xlabel:
        ax.set_xlabel("Elevation (m)", fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_xticklabels([])
        
    if all(parameter is not None for parameter in [ymax, ymin, xmax, xmin]):
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(xmin, xmax)
        
    if ax_legend:
        ax.legend(frameon=True, fontsize=24,
                  loc="upper right", borderaxespad=0,
                  ncol=1) #bbox_to_anchor=(0.01, 1.05, 1, 0.102,)
    else:
        ax.legend([],[], frameon=False)
        
    ax.set_box_aspect(1)
    
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 28, "fontweight":"bold"}, loc="left")
        
    
        
        
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    
# defining coordinates 


# extracting transects 

def extract_all(data, wiso, labels, time="season", season="JJA", season_calendar="standard"):
    
    maxlat_west, minlat_west, maxlon_west, minlon_west = 47, 44, 8, 1
    maxlat_south, minlat_south, maxlon_south, minlon_south = 47, 43, 15, 7.5
    maxlat_north, minlat_north, maxlon_north, minlon_north = 50, 46.5, 16, 5
    
    extracts_north = {}
    extracts_west = {}
    extracts_south = {}
    
    for i,topo in enumerate(labels):
        extracts_north[topo] = extract_vars_and_analysis(data=data[i], wiso=wiso[i],
                                                         maxlon=maxlon_north, maxlat=maxlat_north, minlon=minlon_north, 
                                                         minlat=minlat_north,time=time, season=season, 
                                                         season_calendar=season_calendar)
        
        extracts_west[topo] = extract_vars_and_analysis(data=data[i], wiso=wiso[i],
                                                         maxlon=maxlon_west, maxlat=maxlat_west, minlon=minlon_west, 
                                                         minlat=minlat_west,time=time, season=season, 
                                                         season_calendar=season_calendar)
        
        extracts_south[topo] = extract_vars_and_analysis(data=data[i], wiso=wiso[i],
                                                         maxlon=maxlon_south, maxlat=maxlat_south, minlon=minlon_south, 
                                                         minlat=minlat_south,time=time, season=season, 
                                                         season_calendar=season_calendar)
        
    return extracts_north, extracts_west, extracts_south
        
        
def plot_lapse_rate_JJA():        
        
    labels_W1 = ["CTL", "W1E0", "W1E2",]
    labels_W2 = ["CTL", "W2E1", "W2E0", "W2E2"] 
    
    data_W1 = [aw100e100_data, aw100e0_data, aw100e200_data]
    wiso_W1 = [aw100e100_wiso, aw100e0_wiso, aw100e200_wiso]
    
    data_W2 = [aw100e100_data, aw200e100_data, aw200e0_data, aw200e200_data]
    wiso_W2 = [aw100e100_wiso, aw200e100_wiso, aw200e0_wiso, aw200e200_wiso]  
    
    W1_north, W1_west, W1_south = extract_all(data=data_W1, wiso=wiso_W1, labels=labels_W1)
    W2_north, W2_west, W2_south = extract_all(data=data_W2, wiso=wiso_W2, labels=labels_W2)
    
     
    apply_style(fontsize=28, style="seaborn-paper", linewidth=2,)
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(25, 14), )
    
    plot_section_lapse_rate(dataset=W1_west, makers=["*", "D", "^"], 
                            colors=["black", "red", "green"], 
                            labels=labels_W1, ax=ax1, xlabel=True, ylabel=True, 
                            title="[A] West", ymin=-16, ymax=0, xmax=3500, xmin=0)
    
    plot_section_lapse_rate(dataset=W1_north, makers=["*", "D", "^"], 
                            colors=["black", "red", "green"], 
                            labels=labels_W1, ax=ax2, xlabel=True, ylabel=False, 
                            title="[B] North", ymin=-16, ymax=0, xmax=3500, xmin=0)
    
    
    plot_section_lapse_rate(dataset=W1_south, makers=["*", "D", "^"], 
                            colors=["black", "red", "green"], 
                            labels=labels_W1, ax=ax3, xlabel=True, ylabel=False, 
                            title="[C] South", ymin=-16, ymax=0, xmax=3500, xmin=0)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_plots, "W1_lapse_rate_JJA_ci.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(25, 14), )
    
    plot_section_lapse_rate(dataset=W2_west, makers=["*", "s", "p", "h"], 
                            colors=["black", "darkgoldenrod", "purple", "blue"], 
                            labels=labels_W2, ax=ax1, xlabel=True, ylabel=True, 
                            title="[A] West", ymin=-16, ymax=0, xmax=3500, xmin=0)
    
    plot_section_lapse_rate(dataset=W2_north, makers=["*", "s", "p", "h"], 
                            colors=["black", "darkgoldenrod", "purple", "blue"], 
                            labels=labels_W2, ax=ax2, xlabel=True, ylabel=False, 
                            title="[B] North", ymin=-16, ymax=0, xmax=3500, xmin=0)
    
    
    plot_section_lapse_rate(dataset=W2_south, makers=["*", "s", "p", "h"], 
                            colors=["black", "darkgoldenrod", "purple", "blue"], 
                            labels=labels_W2, ax=ax3, xlabel=True, ylabel=False, 
                            title="[C] South", ymin=-16, ymax=0, xmax=3500, xmin=0)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_plots, "W2_lapse_rate_JJA_ci.svg"), format= "svg", bbox_inches="tight", dpi=600)


# plotting 

def plot_lapse_rate_annual():        
        
    labels_W1 = ["CTL", "W1E0", "W1E2",]
    labels_W2 = ["CTL", "W2E1", "W2E0", "W2E2"] 
    
    data_W1 = [aw100e100_data, aw100e0_data, aw100e200_data]
    wiso_W1 = [aw100e100_wiso, aw100e0_wiso, aw100e200_wiso]
    
    data_W2 = [aw100e100_data, aw200e100_data, aw200e0_data, aw200e200_data]
    wiso_W2 = [aw100e100_wiso, aw200e100_wiso, aw200e0_wiso, aw200e200_wiso]  
    
    W1_north, W1_west, W1_south = extract_all(data=data_W1, wiso=wiso_W1, labels=labels_W1,
                                              time="annual")
    W2_north, W2_west, W2_south = extract_all(data=data_W2, wiso=wiso_W2, labels=labels_W2,
                                              time="annual")
    
     
    apply_style(fontsize=28, style="seaborn-paper", linewidth=2,)
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(25, 14), )
    
    plot_section_lapse_rate(dataset=W1_west, makers=["*", "D", "^"], 
                            colors=["black", "red", "green"], 
                            labels=labels_W1, ax=ax1, xlabel=True, ylabel=True, 
                            title="[A] West", ymin=-20, ymax=0, xmax=3500, xmin=0)
    
    plot_section_lapse_rate(dataset=W1_north, makers=["*", "D", "^"], 
                            colors=["black", "red", "green"], 
                            labels=labels_W1, ax=ax2, xlabel=True, ylabel=False, 
                            title="[B] North", ymin=-20, ymax=0, xmax=3500, xmin=0)
    
    
    plot_section_lapse_rate(dataset=W1_south, makers=["*", "D", "^"], 
                            colors=["black", "red", "green"], 
                            labels=labels_W1, ax=ax3, xlabel=True, ylabel=False, 
                            title="[C] South", ymin=-20, ymax=0, xmax=3500, xmin=0)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_plots, "W1_lapse_rate_annual_ci.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(25, 14), )
    
    plot_section_lapse_rate(dataset=W2_west, makers=["*", "s", "p", "h"], 
                            colors=["black", "darkgoldenrod", "purple", "blue"], 
                            labels=labels_W2, ax=ax1, xlabel=True, ylabel=True, 
                            title="[A] West", ymin=-20, ymax=0, xmax=3500, xmin=0)
    
    plot_section_lapse_rate(dataset=W2_north, makers=["*", "s", "p", "h"], 
                            colors=["black", "darkgoldenrod", "purple", "blue"], 
                            labels=labels_W2, ax=ax2, xlabel=True, ylabel=False, 
                            title="[B] North", ymin=-20, ymax=0, xmax=3500, xmin=0)
    
    
    plot_section_lapse_rate(dataset=W2_south, makers=["*", "s", "p", "h"], 
                            colors=["black", "darkgoldenrod", "purple", "blue"], 
                            labels=labels_W2, ax=ax3, xlabel=True, ylabel=False, 
                            title="[C] South", ymin=-20, ymax=0, xmax=3500, xmin=0)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_plots, "W2_lapse_rate_annual_ci.svg"), format= "svg", bbox_inches="tight", dpi=600)


#plot_lapse_rate_JJA()
plot_lapse_rate_annual()

