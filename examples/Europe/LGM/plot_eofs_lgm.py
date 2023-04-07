# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:15:26 2023

@author: dboateng

Plot the performed EOF analysis for the pmip models for MH
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

from path_to_data_lgm import *




# read data

def read_eof_var(eof_path, vars_path):
    eof = xr.open_dataarray(eof_path)
    vars_data = pd.read_csv(vars_path, index_col=["mode"])
    
    return eof, vars_data 
    
awi_eof, awi_vars = read_eof_var(eof_path=AWI_DJF_EOFs , vars_path=AWI_DJF_VARS)
cesm_eof, cesm_vars = read_eof_var(eof_path=CESM_WA_DJF_EOFs , vars_path=CESM_WA_DJF_VARS)
inm_eof, inm_vars = read_eof_var(eof_path=INM_DJF_EOFs , vars_path=INM_DJF_VARS)
miroc_eof, miroc_vars = read_eof_var(eof_path=MIROC_DJF_EOFs , vars_path=MIROC_DJF_VARS)
mpi_esm_eof, mpi_esm_vars = read_eof_var(eof_path=MPI_ESM_DJF_EOFs , vars_path=MPI_ESM_DJF_VARS)




def plot_NAO_MH():
    apply_style(fontsize=23, style="seaborn-talk", linewidth=3,) 
    #apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.AlbersEqualArea(
        central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))
    
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(25, 22), 
                                                                               subplot_kw={"projection": projection})
    
    modes = [1,2]
    units="hPa" 
    variable="Mean Sea Level Pressure"
    vmax=20
    vmin=-20
    
    figname ="NAO_LGM_DJF_PMIP"
    
    plot_eofsAsCovariance(variable= variable, data=awi_eof.sel(mode=1) * -1, mode_var=awi_vars.loc[1], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                          ax=ax1, fig=fig, title="AWI-ESM-1-1-LR", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=cesm_eof.sel(mode=1) * -1, mode_var=cesm_vars.loc[1], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax2, fig=fig, title="CESM2-WACCM-FV2", bottom_labels=True)
    
    
    plot_eofsAsCovariance(variable= variable, data=inm_eof.sel(mode=1) * -1, mode_var=inm_vars.loc[1], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax3, fig=fig, title="INM-CM4-8", bottom_labels=True)
    
    
    plot_eofsAsCovariance(variable= variable, data=miroc_eof.sel(mode=1) * -1, mode_var=miroc_vars.loc[1], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax4, fig=fig, title="MIROC-ES2L", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=mpi_esm_eof.sel(mode=1) * -1, mode_var=mpi_esm_vars.loc[1], units=units, vmax=40,
                          vmin=-40, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax5, fig=fig, title="MPI-ESM1-2-LR", bottom_labels=True)
    
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, figname + ".svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
def plot_EA_MH():
    apply_style(fontsize=23, style="seaborn-talk", linewidth=3,) 
    #apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.AlbersEqualArea(
        central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))
    
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(25, 22), 
                                                                               subplot_kw={"projection": projection})
    
    modes = [1,2]
    units="hPa" 
    variable="Mean Sea Level Pressure"
    vmax=20
    vmin=-20
    
    figname ="EA_LGM_DJF_PMIP"
    
    plot_eofsAsCovariance(variable= variable, data=awi_eof.sel(mode=3), mode_var=awi_vars.loc[3], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                          ax=ax1, fig=fig, title="AWI-ESM-1-1-LR", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=cesm_eof.sel(mode=3), mode_var=cesm_vars.loc[3], units=units,vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax2, fig=fig, title="CESM2-WACCM-FV2", bottom_labels=True)
    
    
    plot_eofsAsCovariance(variable= variable, data=inm_eof.sel(mode=3), mode_var=inm_vars.loc[3], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax3, fig=fig, title="INM-CM4-8", bottom_labels=True)
    
    
    
    plot_eofsAsCovariance(variable= variable, data=miroc_eof.sel(mode=4), mode_var=miroc_vars.loc[4], units=units, vmax=vmax,
                          vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax4, fig=fig, title="MIROC-ES2L", bottom_labels=True)
    
    plot_eofsAsCovariance(variable= variable, data=mpi_esm_eof.sel(mode=3), mode_var=mpi_esm_vars.loc[3], units=units, vmax=40,
                          vmin=-40, cmap=RdBu_r, domain="NH", levels=22,
                          level_ticks=11, cbar=False, use_AlberEqualArea=True,
                          ax=ax5, fig=fig, title="MPI-ESM1-2-LR", bottom_labels=True)
    
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, figname + ".svg"), format= "svg", bbox_inches="tight", dpi=600)
    
plot_EA_MH()
plot_NAO_MH()
