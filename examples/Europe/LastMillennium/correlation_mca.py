# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 19:24:54 2023

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
from pyClimat.analysis import extract_transect
from pyClimat.variables import extract_var
from pyClimat.utils import extract_region


from path_to_data_lm import *

main_path = "D:/Datasets/iGCM_datasets/"




def perform_correlation_nao_ea_with_t2m_prec(pcs_path, path_to_data, nao_mode, ea_mode, 
                                             nao_factor=1, ea_factor=1, filename=None,
                                             select_time=False, 
                                             time_range=None, use_slice=False, slice_start=None, 
                                             slice_end=None, slice_end_plus=None):
    
    
    df = pd.read_csv(pcs_path, parse_dates=["time"])
    
    NAO_indices = xr.DataArray(df[str(nao_mode)] * nao_factor, dims="time", coords={"time": df["time"]})
    EA_indices = xr.DataArray(df[str(ea_mode)] * ea_factor, dims="time", coords={"time": df["time"]})
    
    if select_time:
        NAO_indices = NAO_indices.sel(time=slice(slice_start, slice_end_plus))
        EA_indices = EA_indices.sel(time=slice(slice_start, slice_end_plus))
        
    filename_data = filename + "_wiso_vars.nc"
    
    temp = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="tsurf")  - 273.15 #°C
    
    d18O = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="d18O")  
    prec = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="prec")
    
    
    
    maxlon = 150
    minlon = -140
    maxlat = 86
    minlat = 24
    
    t2m_season = extract_region(data=temp, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF", select_time=select_time, time_range=time_range, use_slice=use_slice, 
                                 slice_start=slice_start, slice_end=slice_end) 
    
    prec_season = extract_region(data=prec, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF", select_time=select_time, time_range=time_range, use_slice=use_slice, 
                                 slice_start=slice_start, slice_end=slice_end) 
    
    d18O_season = extract_region(data=d18O, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF", select_time=select_time, time_range=time_range, use_slice=use_slice, 
                                 slice_start=slice_start, slice_end=slice_end)
    
    
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
    
    apply_style(fontsize=28, style=None, linewidth=2.5)
    #projection = ccrs.PlateCarree()
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    fig, (ax1,ax2, ax3) = plt.subplots(nrows = 3, ncols = 1, figsize=(22, 28), 
                                                                               subplot_kw={"projection": projection})

    axes = [ax1,ax2, ax3]
    
    for i,label in enumerate(labels):
        if i == 0:

            plot_correlation(variable="Spearman Coefficients", data=data[i].get(reg_name), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="NH Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                             plot_pvalues=True, pvalue_data=data[i].get(sig_name), bottom_labels=True, title=label,
                             plot_projection=projection, extend="both")
            
        else:
            

            plot_correlation(variable="Spearman Coefficients", data=data[i].get(reg_name), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="NH Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, plot_pvalues=True, pvalue_data=data[i].get(sig_name),
                             title=label, plot_projection=projection, extend="both")

    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, figname), format= "svg", bbox_inches="tight", dpi=300)



# date_range_mca = xr.cftime_range(start="0851", end ="1160", freq="MS", calendar="noleap")

# date_range_lia = xr.cftime_range(start="1489", end ="1850", freq="MS", calendar="noleap")  
  
    
cesm_data  = perform_correlation_nao_ea_with_t2m_prec(pcs_path= CESM_DJF_PCs, path_to_data=main_path, nao_mode=1, 
                                                  ea_mode=3, nao_factor=-1, ea_factor=1, filename="CESM",
                                                  select_time=True, 
                                                  time_range=None, use_slice=True, slice_start="0850", 
                                                  slice_end="1160", slice_end_plus="1161")


giss_data  = perform_correlation_nao_ea_with_t2m_prec(pcs_path= GISS_DJF_PCs, path_to_data=main_path, nao_mode=1, 
                                                 ea_mode=3, nao_factor=-1, ea_factor=1, filename="GISS",
                                                 select_time=True, 
                                                 time_range=None, use_slice=True, slice_start="0850", 
                                                 slice_end="1160", slice_end_plus="1161")

hadcm3_data  = perform_correlation_nao_ea_with_t2m_prec(pcs_path=HADESM_DJF_PCs , path_to_data=main_path, nao_mode=1, 
                                                 ea_mode=3, nao_factor=-1, ea_factor=1, filename="HADCM3",
                                                 select_time=True, 
                                                 time_range=None, use_slice=True, slice_start="0850", 
                                                 slice_end="1160", slice_end_plus="1161")


labels = ["iCESM [MCA]", "GISS-E2-R [MCA]", "iHadCM3 [MCA]"]
data = [cesm_data, giss_data, hadcm3_data]


apply_style(fontsize=28, style=None, linewidth=2.5)
 
plot_climate_index_correlation(labels=labels, data=data, figname="temp_nao_nh_mca.svg", reg_name="nao_reg_temp", 
                               sig_name="nao_sig_temp")
plot_climate_index_correlation(labels=labels, data=data, figname="temp_ea_nh_mca.svg", reg_name="ea_reg_temp", 
                               sig_name="ea_sig_temp")

plot_climate_index_correlation(labels=labels, data=data, figname="prec_nao_nh_mca.svg", reg_name="nao_reg_prec", 
                               sig_name="nao_sig_prec")
plot_climate_index_correlation(labels=labels, data=data, figname="prec_ea_nh_mca.svg", reg_name="ea_reg_prec", 
                               sig_name="ea_sig_prec")

plot_climate_index_correlation(labels=labels, data=data, figname="d18O_nao_nh_mca.svg", reg_name="nao_reg_d18O", 
                               sig_name="nao_sig_d18O")
plot_climate_index_correlation(labels=labels, data=data, figname="d18O_ea_nh_mca.svg", reg_name="ea_reg_d18O", 
                               sig_name="ea_sig_d18O")

plt.show()



