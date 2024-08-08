# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 18:50:31 2023

@author: dboateng
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_transect, compute_lterm_diff, compute_lterm_mean
from pyClimat.variables import extract_var
from pyClimat.utils import extract_region


from path_to_data_lm import *

main_path = "D:/Datasets/iGCM_datasets/"

filenames = ["CESM", "ECHAM5", "GISS", "HADCM3", "CCSM"]



def extract_climatologies(filename, time="annual", season=None):
    
    filename_data = filename + "_wiso_vars.nc"
    
    temp = read_from_path(path=main_path, filename=filename_data, decode=True,
                      varname="tsurf")  - 273.15 #°C
    
    d18O = read_from_path(path=main_path, filename=filename_data, decode=True,
                      varname="d18O")  
    
    prec = read_from_path(path=main_path, filename=filename_data, decode=True,
                      varname="prec")
    
    
    date_range_mca = xr.cftime_range(start="0951", end ="1250", freq="MS", calendar="noleap")
    
    date_range_lia = xr.cftime_range(start="1650", end ="1850", freq="MS", calendar="noleap")
    
    if filename == "GISS": #slice only works for HadGEM3 temp.sel(time=slice("1000-01", "1250-01"))
        date_range_mca = xr.cftime_range(start="0951", end ="1250", freq="M", calendar="noleap")
        
        date_range_lia = xr.cftime_range(start="1650", end ="1849", freq="M", calendar="noleap")
        
    elif filename == "HADCM3":
        date_range_mca = temp.sel(time= slice("0951-01", "1250-01")).time
        date_range_lia = temp.sel(time= slice("1650-01", "1850-01")).time
        
        
    
    #compute climatologies
    temp_mca = compute_lterm_mean(data=temp, time=time, time_range=date_range_mca, 
                                  season=season)
    temp_lia = compute_lterm_mean(data=temp, time=time, time_range=date_range_lia, 
                                  season=season)
    temp_diff = compute_lterm_diff(data_control=temp, data_main=temp, time=time,
                                           time_range_control=date_range_mca, 
                                           time_range_main=date_range_lia,
                                           season=season)
    
    
    d18O_mca = compute_lterm_mean(data=d18O, time=time, time_range=date_range_mca, 
                                  season=season)
    d18O_lia = compute_lterm_mean(data=d18O, time=time, time_range=date_range_lia, 
                                  season=season)
    d18O_diff = compute_lterm_diff(data_control=d18O, data_main=d18O, time=time,
                                           time_range_control=date_range_mca, 
                                           time_range_main=date_range_lia,
                                           season=season)
    
    
    prec_mca = compute_lterm_mean(data=prec, time=time, time_range=date_range_mca, 
                                  season=season)
    prec_lia = compute_lterm_mean(data=prec, time=time, time_range=date_range_lia, 
                                  season=season)
    prec_diff = compute_lterm_diff(data_control=prec, data_main=temp, time=time,
                                           time_range_control=date_range_mca, 
                                           time_range_main=date_range_lia,
                                           season=season)
    
    data_mca ={"temp": temp_mca, "d18O": d18O_mca, "prec": prec_mca}
    
    data_lia ={"temp": temp_lia, "d18O": d18O_lia, "prec": prec_lia}
    
    data_diff ={"temp": temp_diff, "d18O": d18O_diff, "prec": prec_diff}
    
    return data_mca, data_lia, data_diff


# test function 
#cesm_mca, cesm_lia, cesm_diff = extract_climatologies(filename="CESM", time="annual")
#giss_mca, giss_lia, giss_diff = extract_climatologies(filename="GISS", time="annual")
hadcm3_mca, hadcm3_lia, hadcm3_diff = extract_climatologies(filename="HADCM3", time="annual")

def plot_d18Op_global(data_mca, data_lia, data_diff, axes=None, fig=None, axes_cbar=True, labels=None):
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    if axes is None:
        fig,(ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(22,18), subplot_kw={"projection":projection})
        axes = [ax1, ax2, ax3]
    
    
    
            
    plot_annual_mean(variable="$\delta^{18}$Op vs SMOW", data_alt=data_mca.get("d18O"), ax=axes[0],
                     cmap=RdYlBu, units="‰", vmax=2, vmin=-30, 
                    levels=22, level_ticks=11, add_colorbar=axes_cbar, cbar_pos= [0.05, 0.05, 0.35, 0.02], 
                    orientation="horizontal", bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[0], center=False, domain="NH Wide")
    
   
    plot_annual_mean(variable="$\delta^{18}$Op vs SMOW", data_alt=data_lia.get("d18O"), ax=axes[1],
                     cmap=RdYlBu, units="‰", vmax=2, vmin=-30, 
                    levels=22, level_ticks=11, add_colorbar=False, 
                    bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[1], center=False, domain="NH Wide")
    
    plot_annual_mean(variable="$\delta^{18}$Op vs SMOW", data_alt=data_diff.get("d18O"), ax=axes[2],
                     cmap="PiYG", units="‰", vmax=0.5, vmin=-0.5, 
                    levels=22, level_ticks=9, add_colorbar=axes_cbar, cbar_pos= [0.45, 0.05, 0.35, 0.02], 
                    orientation="horizontal", bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[2], center=True, domain="NH Wide",
                    label_format="%.1f")
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    
def plot_temp_global(data_mca, data_lia, data_diff, axes=None, fig=None, axes_cbar=True, labels=None):
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    if axes is None:
        fig,(ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(22,18), subplot_kw={"projection":projection})
        axes = [ax1, ax2, ax3]
    
    
    
            
    plot_annual_mean(variable="Temperature", data_alt=data_mca.get("temp"), ax=axes[0],
                     cmap=Spectral_r, units="°C", vmax=30, vmin=-20, 
                    levels=22, level_ticks=11, add_colorbar=axes_cbar, cbar_pos= [0.05, 0.05, 0.35, 0.02], 
                    orientation="horizontal", bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[0], center=False, domain="NH Wide")
    
   
    plot_annual_mean(variable="Temperature", data_alt=data_lia.get("temp"), ax=axes[1],
                     cmap=Spectral_r, units="°C", vmax=30, vmin=-20, 
                    levels=22, level_ticks=11, add_colorbar=False, 
                    bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[1], center=False, domain="NH Wide")
    
    plot_annual_mean(variable="Temperature", data_alt=data_diff.get("temp"), ax=axes[2],
                     cmap=RdBu_r, units="°C", vmax=1.5, vmin=-1.5, 
                    levels=22, level_ticks=7, add_colorbar=axes_cbar, cbar_pos= [0.45, 0.05, 0.35, 0.02], 
                    orientation="horizontal", bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[2], center=True, domain="NH Wide",
                    label_format="%.1f")
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    
def plot_prec_global(data_mca, data_lia, data_diff, axes=None, fig=None, axes_cbar=True, labels=None):
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    if axes is None:
        fig,(ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(22,18), subplot_kw={"projection":projection})
        axes = [ax1, ax2, ax3]
   
    
    
            
    plot_annual_mean(variable="Precipitation", data_alt=data_mca.get("prec"), ax=axes[0],
                     cmap=YlGnBu, units="mm/month", vmax=250, vmin=0, 
                    levels=22, level_ticks=11, add_colorbar=axes_cbar, cbar_pos= [0.05, 0.05, 0.35, 0.02], 
                    orientation="horizontal", bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[0], center=False, domain="NH Wide")
    
   
    plot_annual_mean(variable="Precipitation", data_alt=data_lia.get("prec"), ax=axes[1],
                     cmap=YlGnBu, units="mm/month", vmax=250, vmin=0, 
                    levels=22, level_ticks=11, add_colorbar=False, 
                    bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[1], center=False, domain="NH Wide")
    
    plot_annual_mean(variable="Precipitation", data_alt=data_diff.get("prec"), ax=axes[2],
                     cmap=BrBG, units="mm/month", vmax=100, vmin=-100, 
                    levels=22, level_ticks=9, add_colorbar=axes_cbar, cbar_pos= [0.45, 0.05, 0.35, 0.02], 
                    orientation="horizontal", bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=True, coast_resolution="110m",
                    plot_projection=projection, title=labels[2], center=True, domain="NH Wide",
                    label_format="%.0f")
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 



# plot
apply_style(fontsize=28, style=None, linewidth=2.5)
projection = ccrs.Robinson(central_longitude=0, globe=None)

def plot_d18Op():
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                     figsize=(30, 18), subplot_kw={"projection": projection})
    
    axes_cesm = [ax1, ax2, ax3]
    axes_giss = [ax4, ax5, ax6]
    axes_hadgem = [ax7, ax8, ax9]
    
    plot_d18Op_global(data_mca=cesm_mca, data_lia=cesm_lia, data_diff=cesm_diff, axes=axes_cesm, 
                      fig=fig, axes_cbar=False, labels = ["(a) MCA   [iCESM]", "(b) LIA", "(c) LIA-MCA"])
    plot_d18Op_global(data_mca=giss_mca, data_lia=giss_lia, data_diff=giss_diff, axes_cbar=False, 
                      axes=axes_giss, fig=fig, labels = ["(d) MCA   [GISS-E2-R]", "(e) LIA", "(f) LIA-MCA"])
    plot_d18Op_global(data_mca=hadcm3_mca, data_lia=hadcm3_lia, data_diff=hadcm3_diff, axes_cbar=True,
                      axes=axes_hadgem, fig=fig, labels = ["(g) MCA   [iHadCM3]", "(h) LIA", "(i) LIA-MCA"])
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10)
    plt.savefig(os.path.join(path_to_plots, "d18O_mca_lia.svg"), format= "svg", bbox_inches="tight", dpi=300)


def plot_temp():
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                     figsize=(30, 18), subplot_kw={"projection": projection})
    
    axes_cesm = [ax1, ax2, ax3]
    axes_giss = [ax4, ax5, ax6]
    axes_hadgem = [ax7, ax8, ax9]
    
    plot_temp_global(data_mca=cesm_mca, data_lia=cesm_lia, data_diff=cesm_diff, axes=axes_cesm, 
                      fig=fig, axes_cbar=False, labels = ["(a) MCA   [iCESM]", "(b) LIA", "(c) LIA-MCA"])
    plot_temp_global(data_mca=giss_mca, data_lia=giss_lia, data_diff=giss_diff, axes_cbar=False, 
                      axes=axes_giss, fig=fig, labels = ["(d) MCA   [GISS-E2-R]", "(e) LIA", "(f) LIA-MCA"])
    plot_temp_global(data_mca=hadcm3_mca, data_lia=hadcm3_lia, data_diff=hadcm3_diff, axes_cbar=True,
                      axes=axes_hadgem, fig=fig, labels = ["(g) MCA   [iHadCM3]", "(h) LIA", "(i) LIA-MCA"])
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10)
    plt.savefig(os.path.join(path_to_plots, "temp_mca_lia.svg"), format= "svg", bbox_inches="tight", dpi=300)


def plot_prec():
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                     figsize=(30, 18), subplot_kw={"projection": projection})
    
    axes_cesm = [ax1, ax2, ax3]
    axes_giss = [ax4, ax5, ax6]
    axes_hadgem = [ax7, ax8, ax9]
    
    plot_prec_global(data_mca=cesm_mca, data_lia=cesm_lia, data_diff=cesm_diff, axes=axes_cesm, 
                      fig=fig, axes_cbar=False, labels = ["(a) MCA   [iCESM]", "(b) LIA", "(c) LIA-MCA"])
    plot_prec_global(data_mca=giss_mca, data_lia=giss_lia, data_diff=giss_diff, axes_cbar=False, 
                      axes=axes_giss, fig=fig, labels = ["(d) MCA   [GISS-E2-R]", "(e) LIA", "(f) LIA-MCA"])
    plot_prec_global(data_mca=hadcm3_mca, data_lia=hadcm3_lia, data_diff=hadcm3_diff, axes_cbar=True,
                      axes=axes_hadgem, fig=fig, labels = ["(g) MCA   [iHadCM3]", "(h) LIA", "(i) LIA-MCA"])
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10)
    plt.savefig(os.path.join(path_to_plots, "prec_mca_lia.svg"), format= "svg", bbox_inches="tight", dpi=300)
    
    
if __name__ == "__main__":
    plot_d18Op()
    plot_temp()
    plot_prec()