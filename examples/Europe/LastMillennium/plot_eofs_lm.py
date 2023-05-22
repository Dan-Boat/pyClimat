# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:55:38 2023

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

from path_to_data_lm import *




def read_eof_var(eof_path, vars_path):
    eof = xr.open_dataarray(eof_path)
    vars_data = pd.read_csv(vars_path, index_col=["mode"])
    
    return eof, vars_data 


cesm_eof, cesm_vars = read_eof_var(eof_path=CESM_DJF_EOFs, vars_path=CESM_DJF_VARs)
#ccsm_eof, ccsm_vars = read_eof_var(eof_path=CCSM_DJF_EOFs, vars_path=CCSM_DJF_VARs)
giss_eof, giss_vars = read_eof_var(eof_path=GISS_DJF_EOFs , vars_path=GISS_DJF_VARs)
echam_eof, echam_vars = read_eof_var(eof_path=ECHAM_DJF_EOFs , vars_path=ECHAM_DJF_VARs)
hadesm_eof, hadesm_vars = read_eof_var(eof_path=HADESM_DJF_EOFs , vars_path=HADESM_DJF_VARs)




def plot_NAO():
    apply_style(fontsize=23, style="seaborn-talk", linewidth=3,) 
    #apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.AlbersEqualArea(
        central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))
    
    fig, ((ax1,ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(24, 22), 
                                                                               subplot_kw={"projection": projection})
    
    units="hPa" 
    variable="Mean Sea Level Pressure"
    vmax=50
    vmin=-50
    
    figname ="NAO_DJF_LastMillennium"
    
    plot_eofsAsCovariance(variable= variable, data=cesm_eof.sel(mode=1) * -1, mode_var=cesm_vars.loc[1], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                          ax=ax1, fig=fig, title="iCESM", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=giss_eof.sel(mode=1) * -1, mode_var=giss_vars.loc[1], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax2, fig=fig, title="GISS-E2-R", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=echam_eof.sel(mode=1) * -1, mode_var=echam_vars.loc[1], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax3, fig=fig, title="ECHAM5-wiso", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=hadesm_eof.sel(mode=1) * -1, mode_var=hadesm_vars.loc[1], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax4, fig=fig, title="iHadCM3", bottom_labels=True)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, figname + ".svg"), format= "svg", bbox_inches="tight", dpi=600)


def plot_EA():
    apply_style(fontsize=23, style="seaborn-talk", linewidth=3,) 
    #apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.AlbersEqualArea(
        central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))
    
    fig, ((ax1,ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(24, 22), 
                                                                               subplot_kw={"projection": projection})
    
    units="hPa" 
    variable="Mean Sea Level Pressure"
    vmax=50
    vmin=-50
    
    figname ="EA_DJF_LastMillennium"
    
    plot_eofsAsCovariance(variable= variable, data=cesm_eof.sel(mode=3), mode_var=cesm_vars.loc[3], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                          ax=ax1, fig=fig, title="iCESM", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=giss_eof.sel(mode=3), mode_var=giss_vars.loc[3], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax2, fig=fig, title="GISS-E2-R", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=echam_eof.sel(mode=4), mode_var=echam_vars.loc[4], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax3, fig=fig, title="ECHAM5-wiso", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=hadesm_eof.sel(mode=3), mode_var=hadesm_vars.loc[3], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax4, fig=fig, title="iHadCM3", bottom_labels=True)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, figname + ".svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
    
plot_EA()
plot_NAO()