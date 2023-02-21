# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:47:18 2023

@author: dboateng
The control scripts of all the routine functions
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


def extract_eofs_data(data, figname, units, variable, vmax=15, vmin=-15, plot_covariance=True, is_era=False,
                      path_to_plots=None, apply_varimax=True, save_files=False,
                      path_to_files=None, filename=None):
    # analysis for ERA5 Dataset
    
    
    
    # initiate the eof instance 
    EOF = EOF_standard(data=data, weights=True, standardize=True, 
                          extract_region=True, extract_season=True, neofs=4)
    
    # select the region of interest and season
    EOF.select_time_and_region(maxlon=60, minlon=-80, maxlat=80, minlat=10, time="season", 
                                  season="DJF", month="JA") # month="AMJJAS"
    
    # calculate the anomalies and apply norm
    EOF.calculate_anomalies(monthly_anomalies=True)
    
    method = "xeofs"
    
    # fit the eof with the solver
    EOF.eof_solver(method=method, apply_promax=False, apply_varimax=apply_varimax)
    
    # extract the eofs (as the eigenvectors of the covariance matric of the X field)
    
    eofs = EOF.eofs(eofscaling=2) # 2 - multiply by the singular values or square root of the eigen values
    eofs_corr, eofs_pvals = EOF.eofs_as_correlation()
    pcs = EOF.pcs(pscaling=0)
    variance = EOF.explained_variance_ratio()
    
    if save_files:
        
        eofs.to_netcdf(os.path.join(path_to_files, filename + "_eofsAsCovariance.nc"))
        eofs_corr.to_netcdf(os.path.join(path_to_files, filename + "_eofsAsCorrelation.nc"))
        eofs_pvals.to_netcdf(os.path.join(path_to_files, filename + "_eofsAsCorrelation_pvals.nc"))
        pcs.to_csv(os.path.join(path_to_files, filename + "_pcs.csv"))
        variance.to_csv(os.path.join(path_to_files, filename + "_pcs_variance.csv"))
        
    if plot_covariance:
    # loop through this !!
        plot_eofs(data=eofs, variance=variance, figname=figname, units=units, variable=variable, vmax=vmax, 
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