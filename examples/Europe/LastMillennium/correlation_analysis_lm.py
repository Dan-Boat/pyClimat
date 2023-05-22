# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:56:05 2023

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
from pyClimat.plots import plot_eofsAsCovariance, plot_correlation
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_var, extract_transect
from pyClimat.utils import extract_region


from path_to_data_lm import *

main_path = "D:/Datasets/iGCM_datasets/"




def perform_correlation_nao_ea_with_t2m_prec(pcs_path, path_to_data, nao_mode, ea_mode, 
                                             nao_factor=1, ea_factor=1, filename=None):
    
    
    df = pd.read_csv(pcs_path, parse_dates=["time"])
    
    NAO_indices = xr.DataArray(df[str(nao_mode)] * nao_factor, dims="time", coords={"time": df["time"]})
    EA_indices = xr.DataArray(df[str(ea_mode)] * ea_factor, dims="time", coords={"time": df["time"]})
    
    filename_data = filename + "_wiso_vars.nc"
    
    temp = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="tsurf")  - 273.15 #°C
    
    d18O = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="d18O")  
    prec = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="prec")
    
    
    
    
    t2m_season = extract_region(data=temp, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF") 
    
    prec_season = extract_region(data=prec, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF") 
    
    d18O_season = extract_region(data=d18O, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF")
    
    
    nao_t2m_sval, nao_t2m_pval, nao_t2m_sig = StatCorr(x=t2m_season, y=NAO_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    ea_t2m_sval, ea_t2m_pval, ea_t2m_sig = StatCorr(x=t2m_season, y=EA_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    nao_prec_sval, nao_prec_pval, nao_prec_sig = StatCorr(x=prec_season, y=NAO_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    ea_prec_sval, ea_prec_pval, ea_prec_sig = StatCorr(x=prec_season, y=EA_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    nao_d18O_sval, nao_d18O_pval, nao_d18O_sig = StatCorr(x=d18O_season, y=NAO_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    ea_d18O_sval, ea_d18O_pval, ea_d18O_sig = StatCorr(x=d18O_season, y=EA_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    return_nao_reg = {"nao_reg_temp":nao_t2m_sval, "nao_reg_prec":nao_prec_sval, "nao_reg_d18O":nao_d18O_sval}
    
    return_nao_sig = {"nao_sig_temp":nao_t2m_sig, "nao_sig_prec":nao_prec_sig, "nao_sig_d18O":nao_d18O_sig}
    
    return_ea_reg = {"ea_reg_temp":ea_t2m_sval, "ea_reg_prec":ea_prec_sval, "ea_reg_d18O":ea_d18O_sval}
    
    return_ea_sig = {"ea_sig_temp":ea_t2m_sig, "ea_sig_prec":ea_prec_sig, "ea_sig_d18O":ea_d18O_sig}
    
    return_data = return_nao_reg | return_nao_sig | return_ea_reg | return_ea_sig
    
    return return_data

def plot_climate_index_correlation(labels, data, figname, reg_name="nao_reg_temp", sig_name="nao_sig_temp"):
    
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(24, 18), 
                                                                               subplot_kw={"projection": projection})

    axes = [ax1,ax2, ax3]
    
    for i,label in enumerate(labels):
        if i == 0:

            plot_correlation(variable="Spearman Coefficients", data=data[i].get(reg_name), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                             plot_pvalues=True, pvalue_data=data[i].get(sig_name), bottom_labels=True, title=label)
            
        else:
            

            plot_correlation(variable="Spearman Coefficients", data=data[i].get(reg_name), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, plot_pvalues=True, pvalue_data=data[i].get(sig_name),
                             title=label)

    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, figname), format= "svg", bbox_inches="tight", dpi=300)
    
    
echam_data  = perform_correlation_nao_ea_with_t2m_prec(pcs_path= ECHAM_DJF_PCs, path_to_data=main_path, nao_mode=1, 
                                                 ea_mode=4, nao_factor=-1, ea_factor=1, filename="ECHAM5")

cesm_data  = perform_correlation_nao_ea_with_t2m_prec(pcs_path= CESM_DJF_PCs, path_to_data=main_path, nao_mode=1, 
                                                 ea_mode=3, nao_factor=-1, ea_factor=1, filename="CESM")

giss_data  = perform_correlation_nao_ea_with_t2m_prec(pcs_path= GISS_DJF_PCs, path_to_data=main_path, nao_mode=1, 
                                                 ea_mode=3, nao_factor=-1, ea_factor=1, filename="GISS")

hadcm3_data  = perform_correlation_nao_ea_with_t2m_prec(pcs_path=HADESM_DJF_PCs , path_to_data=main_path, nao_mode=1, 
                                                 ea_mode=3, nao_factor=-1, ea_factor=1, filename="HADCM3")


labels = ["iCESM", "GISS-E2-R", "iHadCM3"]
data = [cesm_data, giss_data, hadcm3_data]

plot_climate_index_correlation(labels=labels, data=data, figname="temp_nao.svg", reg_name="nao_reg_temp", 
                               sig_name="nao_sig_temp")
plot_climate_index_correlation(labels=labels, data=data, figname="temp_ea.svg", reg_name="ea_reg_temp", 
                               sig_name="ea_sig_temp")

plot_climate_index_correlation(labels=labels, data=data, figname="prec_nao.svg", reg_name="nao_reg_prec", 
                               sig_name="nao_sig_prec")
plot_climate_index_correlation(labels=labels, data=data, figname="prec_ea.svg", reg_name="ea_reg_prec", 
                               sig_name="ea_sig_prec")

plot_climate_index_correlation(labels=labels, data=data, figname="d18O_nao.svg", reg_name="nao_reg_d18O", 
                               sig_name="nao_sig_d18O")
plot_climate_index_correlation(labels=labels, data=data, figname="d18O_ea.svg", reg_name="ea_reg_d18O", 
                               sig_name="ea_sig_d18O")

plt.show()



