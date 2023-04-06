# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 18:55:06 2023

@author: dboateng
"""

import os 
import numpy as np 
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt 
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.analysis import compute_lterm_diff, compute_lterm_mean
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean

from pyClimat.data import read_from_path

main_path_mh = "D:/Datasets/CMIP6/PMIP/postprocessed/MH"
main_path_lgm = "D:/Datasets/CMIP6/PMIP/postprocessed/LGM"
main_path_plio = "D:/Datasets/CMIP6/PMIP/postprocessed/PLIO"
main_path_pi = "D:/Datasets/CMIP6/PMIP/postprocessed/PI"

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"


#PI

pi_awi_path = os.path.join(main_path_pi, "AWI-ESM-1-1-LR")
pi_cesm_path = os.path.join(main_path_pi, "CESM2")
pi_cesm_waccm_path = os.path.join(main_path_pi, "CESM2-WACCM-FV2")
pi_ec_earth_path = os.path.join(main_path_pi, "EC-Earth3-LR")
pi_giss_path = os.path.join(main_path_pi, "GISS-E2-1-G")
pi_hadGEM_path = os.path.join(main_path_pi, "HadGEM3-GC31-LL")
pi_inm_cm_path = os.path.join(main_path_pi, "INM-CM4-8")
pi_ipsl_path = os.path.join(main_path_pi, "IPSL-CM6A-LR")
pi_miroc_path = os.path.join(main_path_pi, "MIROC-ES2L")
pi_mpi_esm_path = os.path.join(main_path_pi, "MPI-ESM1-2-LR")
pi_norESM_path = os.path.join(main_path_pi, "NorESM1-F")


#MH
mh_awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
mh_cesm_path = os.path.join(main_path_mh, "CESM2")
mh_ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
mh_giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
mh_hadGEM_path = os.path.join(main_path_mh, "HadGEM3-GC31-LL")
mh_ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
mh_miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mh_mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")

#LGM 
lgm_awi_path = os.path.join(main_path_lgm, "AWI-ESM-1-1-LR")
lgm_cesm_waccm_path = os.path.join(main_path_lgm, "CESM2-WACCM-FV2")
lgm_inm_cm_path = os.path.join(main_path_lgm, "INM-CM4-8")
lgm_miroc_path = os.path.join(main_path_lgm, "MIROC-ES2L")
lgm_mpi_esm_path = os.path.join(main_path_lgm, "MPI-ESM1-2-LR")



def compute_anomalies(pi_path, climate_path):
    
    psl_pi = read_from_path(path=pi_path, filename="psl_1m_lterm.nc", 
                          varname="psl", decode=True) / 100 #Pa --> hPa
    
    psl_main = read_from_path(path=climate_path, filename="psl_1m_lterm.nc", 
                          varname="psl", decode=True) / 100 #Pa --> hPa
    
    t2m_pi = read_from_path(path=pi_path, filename="tas_1m_lterm.nc", decode=True,varname="tas",
                              ) - 273.15 #°C
    prec_pi = read_from_path(path=pi_path, filename="pr_1m_lterm.nc", decode=True, 
                               varname="pr",) *60*60*24*30  #mm/month
    
    t2m_main = read_from_path(path=climate_path, filename="tas_1m_lterm.nc", decode=True,varname="tas",
                              ) - 273.15 #°C
    prec_main = read_from_path(path=climate_path, filename="pr_1m_lterm.nc", decode=True, 
                               varname="pr",) *60*60*24*30  #mm/month
    
    psl_diff = compute_lterm_diff(data_control=psl_pi, data_main=psl_main, time= "season", season = "DJF")
    
    t2m_diff = compute_lterm_diff(data_control=t2m_pi, data_main=t2m_main, time= "season", season = "DJF")
    
    prec_diff = compute_lterm_diff(data_control=prec_pi, data_main=prec_main, time= "season", season = "DJF")
    
    return psl_diff, t2m_diff, prec_diff


# calculate anomalies MH
mh_awi_psl_diff, mh_awi_t2m_diff, mh_awi_prec_diff = compute_anomalies(pi_path=pi_awi_path, climate_path=mh_awi_path)
mh_cesm_psl_diff, mh_cesm_t2m_diff, mh_cesm_prec_diff = compute_anomalies(pi_path=pi_cesm_path, climate_path=mh_cesm_path)
mh_ec_earth_psl_diff, mh_ec_earth_t2m_diff, mh_ec_earth_prec_diff = compute_anomalies(pi_path=pi_ec_earth_path, climate_path=mh_ec_earth_path)
mh_giss_psl_diff, mh_giss_t2m_diff, mh_giss_prec_diff = compute_anomalies(pi_path=pi_giss_path, climate_path=mh_giss_path)
mh_ipsl_psl_diff, mh_ipsl_t2m_diff, mh_ipsl_prec_diff = compute_anomalies(pi_path=pi_ipsl_path, climate_path=mh_ipsl_path)
mh_miroc_psl_diff, mh_miroc_t2m_diff, mh_miroc_prec_diff = compute_anomalies(pi_path=pi_miroc_path, climate_path=mh_miroc_path)
mh_mpi_esm_psl_diff, mh_mpi_esm_t2m_diff, mh_mpi_esm_prec_diff = compute_anomalies(pi_path=pi_mpi_esm_path, climate_path=mh_mpi_esm_path)

def plot_prec_anomalies():
    mh_psl_data = [mh_awi_psl_diff, mh_cesm_psl_diff, mh_ec_earth_psl_diff, mh_giss_psl_diff, mh_ipsl_psl_diff, mh_miroc_psl_diff,
                mh_mpi_esm_psl_diff]
    
    mh_prec_data = [mh_awi_prec_diff, mh_cesm_prec_diff, mh_ec_earth_prec_diff, mh_giss_prec_diff, mh_ipsl_prec_diff, mh_miroc_prec_diff,
                mh_mpi_esm_prec_diff]
    
    
    projection = ccrs.PlateCarree()
    apply_style(fontsize=23, style="seaborn-talk", linewidth=2,)
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols = 3, figsize=(28, 22), 
                                                                           subplot_kw={"projection": projection})
    
    
    
    labels_mh = ["AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
                  "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"] 
    
    axes = [ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    
    for i,label in enumerate(labels_mh):
        if i == 0:                                                                         
    
      
            plot_annual_mean(variable="Precipitation anomaly", data_alt=mh_prec_data[i], cmap=BrBG, units="mm/month",
                             ax=axes[i], fig=fig, vmax=50, vmin=-50,levels=22, domain="NH", level_ticks=11, add_colorbar=True,
                             cbar_pos= [0.30, 0.07, 0.30, 0.02], orientation="horizontal", plot_contour=True, c_data= mh_psl_data[i],
                             c_vmax=2, c_vmin=-2, c_levels=9, title=label)  
            
        else:
            plot_annual_mean(variable="Precipitation anomaly", data_alt=mh_prec_data[i], cmap=BrBG, units="mm/month",
                             ax=axes[i], fig=fig, vmax=50, vmin=-50,levels=22, domain="NH", level_ticks=11, add_colorbar=False,
                             plot_contour=True, c_data= mh_psl_data[i],
                             c_vmax=2, c_vmin=-2, c_levels=9, title=label) 
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, "prec_anomalies_mh.svg"), format= "svg", bbox_inches="tight", dpi=300)


def plot_temp_anomalies():
    mh_psl_data = [mh_awi_psl_diff, mh_cesm_psl_diff, mh_ec_earth_psl_diff, mh_giss_psl_diff, mh_ipsl_psl_diff, mh_miroc_psl_diff,
                mh_mpi_esm_psl_diff]
    
    mh_t2m_data = [mh_awi_t2m_diff, mh_cesm_t2m_diff, mh_ec_earth_t2m_diff, mh_giss_t2m_diff, mh_ipsl_t2m_diff, mh_miroc_t2m_diff,
                mh_mpi_esm_t2m_diff]
    
    apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
    projection = ccrs.PlateCarree()
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols = 3, figsize=(28, 22), 
                                                                           subplot_kw={"projection": projection})
    
    
    
    labels_mh = ["AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
                  "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"] 
    
    axes = [ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    
    for i,label in enumerate(labels_mh):
        if i == 0:                                                                         
    
      
            plot_annual_mean(variable="Temperature anomaly", data_alt=mh_t2m_data[i], cmap=RdBu_r, units="°C",
                             ax=axes[i], fig=fig, vmax=5, vmin=-5,levels=22, domain="NH", level_ticks=11, add_colorbar=True,
                             cbar_pos= [0.30, 0.07, 0.30, 0.02], orientation="horizontal", plot_contour=True, c_data= mh_psl_data[i],
                             c_vmax=2, c_vmin=-2, c_levels=9, title=label)  
            
        else:
            plot_annual_mean(variable="Temperature anomaly", data_alt=mh_t2m_data[i], cmap=RdBu_r, units="°C",
                             ax=axes[i], fig=fig, vmax=5, vmin=-5,levels=22, domain="NH", level_ticks=11, add_colorbar=False,
                             plot_contour=True, c_data= mh_psl_data[i],
                             c_vmax=2, c_vmin=-2, c_levels=9, title=label) 
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, "t2m_anomalies_mh.svg"), format= "svg", bbox_inches="tight", dpi=300)




plot_temp_anomalies()
plot_prec_anomalies()
