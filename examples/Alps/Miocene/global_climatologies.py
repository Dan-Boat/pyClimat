# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:53:03 2023

@author: dboateng
This section analyse the global simulated isotopes, temp, and precipitation for the mio
and the PI simulations
"""

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


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"
W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

# reading data 
# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_PI_data, W1E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)


# extract values and compute long-term means

def extract_vars_and_analysis(data, wiso):
    
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    u10 = extract_var(Dataset=data , varname="u10")
    v10 = extract_var(Dataset=data , varname="v10")
    
    #compute climatologies
    temp2_alt = compute_lterm_mean(data=temp2, time="annual")
    prec_alt = compute_lterm_mean(data=prec, time="annual")
    d18op_alt = compute_lterm_mean(data=d18op, time="annual")
    u10_alt = compute_lterm_mean(data=u10, time="annual")
    v10_alt = compute_lterm_mean(data=v10, time="annual")
    
    
    return_data = {"temperature":temp2_alt, "precipitation":prec_alt, "d18Op":d18op_alt,
                   "Uwinds":u10_alt, "Vwinds":v10_alt}
    
    return return_data
    
# read data 
Mio278_data = extract_vars_and_analysis(data=W1E1_278_data, wiso=W1E1_278_wiso)
Mio450_data = extract_vars_and_analysis(data=W1E1_450_data, wiso=W1E1_450_wiso)
PI_data = extract_vars_and_analysis(data=W1E1_PI_data, wiso=W1E1_PI_wiso)



# plot comments (reduce the resolution of the coastlines)
# plot for isotopes

def plot_d18Op_global(axes=None, fig=None):
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    if axes is None:
        fig,(ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(18,22), subplot_kw={"projection":projection})
        axes = [ax1, ax2, ax3]
    labels = ["(a) PI", "(b) MIO 278ppm", "(c) MIO 450ppm"]
    data = [PI_data ,Mio278_data, Mio450_data]
    
    for i,label in enumerate(labels):
        if i == 0:
            
            plot_annual_mean(variable="$\delta^{18}$Op vs SMOW", data_alt=data[i].get("d18Op"), ax=axes[i],
                             cmap=RdYlBu, units="‰", vmax=2, vmin=-28, 
                            levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.05, 0.05, 0.25, 0.02], 
                            orientation="horizontal", plot_coastlines=True, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, coast_resolution="110m",
                            plot_projection=projection, title=label, center=False)
            
        else:
            plot_annual_mean(variable="$\delta^{18}$Op vs SMOW", data_alt=data[i].get("d18Op"), ax=axes[i],
                             cmap=RdYlBu, units="‰", vmax=2, vmin=-28, 
                            levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, title=label, center=False)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    # plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    # plt.savefig(os.path.join(path_to_plots, "d18Op_global.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
def plot_temperature_global(axes=None, fig=None):
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    if axes is None:
        fig,(ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(18,22), subplot_kw={"projection":projection})
        axes = [ax1, ax2, ax3]
    labels = ["(a) PI", "(b) MIO 278ppm", "(c) MIO 450ppm"]
    data = [PI_data ,Mio278_data, Mio450_data]
    
    for i,label in enumerate(labels):
        if i == 0:
            
            plot_annual_mean(variable="Temperature", data_alt=data[i].get("temperature"), ax=axes[i],
                             cmap=Spectral_r, units="°C", vmax=40, vmin=-10, 
                            levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.35, 0.05, 0.25, 0.02], 
                            orientation="horizontal", plot_coastlines=True, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, coast_resolution="110m",
                            plot_projection=projection, title=label, center=False)
            
        else:
            plot_annual_mean(variable="Temperature", data_alt=data[i].get("temperature"), ax=axes[i],
                             cmap=Spectral_r, units="°C", vmax=40, vmin=-10, 
                            levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, title=label, center=False)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    # plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    # plt.savefig(os.path.join(path_to_plots, "temperature_global.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
def plot_precipitation_global(axes=None, fig=None):
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    if axes is None:
        fig,(ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(18,22), subplot_kw={"projection":projection})
        axes = [ax1, ax2, ax3]
        
    labels = ["(a) PI", "(b) MIO 278ppm", "(c) MIO 450ppm"]
    data = [PI_data ,Mio278_data, Mio450_data]
    
    for i,label in enumerate(labels):
        if i == 0:
            
            plot_annual_mean(variable="Precipitation", data_alt=data[i].get("precipitation"), ax=axes[i],
                             cmap=YlGnBu, units="mm/month", vmax=450, vmin=0, 
                            levels=22, level_ticks=7, add_colorbar=True, cbar_pos= [0.65, 0.05, 0.25, 0.02], 
                            orientation="horizontal", plot_coastlines=True, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, coast_resolution="110m",
                            plot_projection=projection, title=label, center=False)
            
        else:
            plot_annual_mean(variable="Precipitation", data_alt=data[i].get("precipitation"), ax=axes[i],
                             cmap=YlGnBu, units="mm/month", vmax=450, vmin=0, 
                            levels=22, level_ticks=7, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, title=label, center=False)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    # plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    # plt.savefig(os.path.join(path_to_plots, "precipitation_global.svg"), format= "svg", bbox_inches="tight", dpi=600)



def plot_all_on_fig():
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols = 3, figsize=(28, 22),
                                                                           subplot_kw={"projection":  projection})
    axes1=[ax1, ax4, ax7]
    axes2=[ax2, ax5, ax8]
    axes3=[ax3, ax6, ax9]
    
    plot_d18Op_global(axes=axes1, fig=fig)
    plot_precipitation_global(axes=axes2, fig=fig)
    plot_temperature_global(axes=axes3, fig=fig)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, "all_global_climatologies.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
if __name__ == "__main__":
    # plot_precipitation_global()
    # plot_temperature_global()
    # plot_d18Op_global()
    plot_all_on_fig()   
    
    
#plot for temperature 

# plot for precipitatoin (with winds)