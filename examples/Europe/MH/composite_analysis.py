# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 16:26:35 2023

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
from pyClimat.plots import plot_correlation, plot_annual_mean
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_var, extract_transect
from pyClimat.utils import extract_region


from path_to_data_mh import *

main_path_mh = "D:/Datasets/CMIP6/PMIP/postprocessed/MH"
awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
cesm_path = os.path.join(main_path_mh, "CESM2")
ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")


def get_indices(pcs_path, nao_mode, ea_mode, 
                                             nao_factor=1, ea_factor=1):
    
    
    df = pd.read_csv(pcs_path, parse_dates=["time"])
    
    # extract the indices for OP and EG
    
    df_eq_index = (df[str(nao_mode)] * nao_factor > 0) & (df[str(ea_mode)] * ea_factor > 0) | (df[str(nao_mode)] * nao_factor < 0) & (df[str(ea_mode)] * ea_factor < 0)
    df_op_index = (df[str(nao_mode)] * nao_factor > 0) & (df[str(ea_mode)] * ea_factor < 0) | (df[str(nao_mode)] * nao_factor < 0) & (df[str(ea_mode)] * ea_factor > 0)
    
    
    df_eq_positive_index = (df[str(nao_mode)] * nao_factor > 0) & (df[str(ea_mode)] * ea_factor > 0)
    
    # extract them from the pcs 
    
    df_eq = df[df_eq_index]
    df_op = df[df_op_index]
    
    df_eq_positive = df[df_eq_positive_index]
    
    #create xarray for the pcs
    
    
    EQ_indices = xr.DataArray(df_eq[str(nao_mode)] * nao_factor, dims="time", coords={"time": df_eq["time"]})
    EQ_positive_indices = xr.DataArray(df_eq_positive[str(nao_mode)] * nao_factor, dims="time", coords={"time": df_eq_positive["time"]})
    OP_indices = xr.DataArray(df_op[str(nao_mode)] * nao_factor, dims="time", coords={"time": df_op["time"]})
    
    return EQ_indices, EQ_positive_indices, OP_indices


    
def perform_correlation_composite(atmos_index, path_to_data, pcs_path=None):   
    
    if pcs_path is not None:
        df = pd.read_csv(pcs_path, parse_dates=["time"])
        

    t2m_data = read_from_path(path=path_to_data, filename="tas_monthly.nc", decode=True,varname="tas",
                              ) - 273.15 #°C
    prec_data = read_from_path(path=path_to_data, filename="pr_monthly.nc", decode=True, 
                               varname="pr",) *60*60*24*30  #mm/month
    
    t2m_season = extract_region(data=t2m_data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF") 
    
    prec_season = extract_region(data=prec_data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF")
    
    
    t2m_season["time"] = df.time.values
    
    t2m_index = t2m_season.sel(time = t2m_season.time.isin(atmos_index.time.values))
    
    prec_season["time"] = df.time.values
    
    prec_index = prec_season.sel(time = prec_season.time.isin(atmos_index.time.values))
    
    
    
    nao_t2m_sval, nao_t2m_pval, nao_t2m_sig = StatCorr(x=t2m_index, y=atmos_index, dim="time",
                                                       return_sig=True, sig=0.05)
    
    
    
    nao_prec_sval, nao_prec_pval, nao_prec_sig = StatCorr(x=prec_index, y=atmos_index, dim="time",
                                                       return_sig=True, sig=0.05)
    
    
    
    return nao_t2m_sval, nao_t2m_sig, nao_prec_sval, nao_prec_sig


def cal_diff_OP_EQ(atmos_index_op, atmos_index_eq, path_to_data, pcs_path=None):   
    if pcs_path is not None:
        df = pd.read_csv(pcs_path, parse_dates=["time"])
        

    t2m_data = read_from_path(path=path_to_data, filename="tas_monthly.nc", decode=True,varname="tas",
                              ) - 273.15 #°C
    prec_data = read_from_path(path=path_to_data, filename="pr_monthly.nc", decode=True, 
                               varname="pr",) *60*60*24*30  #mm/month
    
    t2m_season = extract_region(data=t2m_data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF") 
    
    prec_season = extract_region(data=prec_data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF")
    
    
    t2m_season["time"] = df.time.values
    t2m_op = t2m_season.sel(time = t2m_season.time.isin(atmos_index_op.time.values)).mean(dim="time")
    t2m_eq = t2m_season.sel(time = t2m_season.time.isin(atmos_index_eq.time.values)).mean(dim="time")
    
    
    prec_season["time"] = df.time.values
    prec_op = prec_season.sel(time = prec_season.time.isin(atmos_index_op.time.values)).mean(dim="time")
    prec_eq = prec_season.sel(time = prec_season.time.isin(atmos_index_eq.time.values)).mean(dim="time")
    
    t2m_diff = t2m_op - t2m_eq
    prec_diff = prec_op - prec_eq
    
    return t2m_diff.sortby("lon"), prec_diff.sortby("lon")
    

    
def extract_all_for_model(pcs_path, path_to_data, nao_mode, ea_mode, 
                                             nao_factor=1, ea_factor=1):
    
    eq, eq_pos, op = get_indices(pcs_path, nao_mode, ea_mode, 
                                                 nao_factor=1, ea_factor=1)
    t2m_eq_sval, t2m_eq_sig, prec_eq_sval, prec_eq_sig = perform_correlation_composite(atmos_index=eq, path_to_data=path_to_data,
                                                                                       pcs_path=pcs_path)
    t2m_op_sval, t2m_op_sig, prec_op_sval, prec_op_sig = perform_correlation_composite(atmos_index=op, path_to_data=path_to_data,
                                                                                       pcs_path=pcs_path)
    t2m_diff, prec_diff = cal_diff_OP_EQ(atmos_index_op=op, atmos_index_eq=eq_pos, path_to_data=path_to_data, pcs_path=pcs_path)
    
    data = {"t2m_eq_sval":t2m_eq_sval, "t2m_eq_sig":t2m_eq_sig, "t2m_op_sval":t2m_op_sval, "t2m_op_sig":t2m_op_sig, "t2m_diff":t2m_diff,
            "prec_eq_sval":prec_eq_sval, "prec_eq_sig":prec_eq_sig, "prec_op_sval":prec_op_sval, "prec_op_sig":prec_op_sig, "prec_diff":prec_diff}
    
    return data


awi_data = extract_all_for_model(pcs_path=AWI_DJF_PCS, path_to_data=awi_path, nao_mode=1, ea_mode=2, nao_factor=-1, ea_factor=-1)  
cesm_data = extract_all_for_model(pcs_path=CESM2_DJF_PCS, path_to_data=cesm_path, nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1) 
ec_earth_data = extract_all_for_model(pcs_path=EC_Earth3_DJF_PCS, path_to_data=ec_earth_path, nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)  
#giss_data = extract_all_for_model(pcs_path=CESM2_DJF_PCS, path_to_data=cesm_path, nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)  
#ipsl_data = extract_all_for_model(pcs_path=CESM2_DJF_PCS, path_to_data=cesm_path, nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)   



labels_mh = ["AWI-ESM-1-1-LR", "CESM2",]# "EC-Earth3-LR", "GISS-E2-1-G", 
              #"IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"]

mh_data = [awi_data, cesm_data]


apply_style(fontsize=22, style=None, linewidth=2)
projection = ccrs.PlateCarree()
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                 figsize=(28, 25), subplot_kw={"projection": projection})

axes_op = [ax2, ax5, ax8]
axes_eq = [ax1, ax4, ax7]
axes_diff = [ax3, ax6, ax9]

for i,label in enumerate(labels_mh):
    if i == 0:
        
        plot_correlation(variable="Spearman Coefficients", data=mh_data[i].get("t2m_eq_sval"), units="-", vmax=1,
                         vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                         level_ticks=7, ax=axes_eq[i], fig=fig, cbar_pos= [0.05, 0.01, 0.30, 0.02],
                         plot_pvalues=True, pvalue_data=mh_data[i].get("t2m_eq_sig"), bottom_labels=True, 
                         title="NAO-t2m (EQ)-"+ label)
        
        plot_correlation(variable="Spearman Coefficients", data=mh_data[i].get("t2m_op_sval"), units="-", vmax=1,
                         vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                         level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=mh_data[i].get("t2m_op_sig"), bottom_labels=False,
                         left_labels=False, title="NAO-t2m (OP)-"+ label)
        
        plot_annual_mean(variable='Temperature', data_alt=mh_data[i].get("t2m_diff"), cmap=RdBu, units="°C",
             ax=axes_diff[i], fig=fig, vmax=4, vmin=-4, levels=22, domain="Europe Wide", level_ticks=11, 
             cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] t2m (OP-EQ)",
             left_labels=False, bottom_labels=False, label_format="%.1f", add_colorbar=True)
    else:
        plot_correlation(variable="Spearman Coefficients", data=mh_data[i].get("t2m_eq_sval"), units="-", vmax=1,
                         vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False,
                         level_ticks=7, ax=axes_eq[i], fig=fig,
                         plot_pvalues=True, pvalue_data=mh_data[i].get("t2m_eq_sig"), bottom_labels=True, 
                         title="NAO-t2m (EQ)-"+ label)
        
        plot_correlation(variable="Spearman Coefficients", data=mh_data[i].get("t2m_op_sval"), units="-", vmax=1,
                         vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                         level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=mh_data[i].get("t2m_op_sig"), bottom_labels=False,
                         left_labels=False, title="NAO-t2m (OP)-"+ label)
        
        plot_annual_mean(variable='Temperature', data_alt=mh_data[i].get("t2m_diff"), cmap=RdBu, units="°C",
             ax=axes_diff[i], fig=fig, vmax=4, vmin=-4, levels=22, domain="Europe Wide", level_ticks=11, 
             cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] t2m (OP-EQ)",
             left_labels=False, bottom_labels=False, label_format="%.1f", add_colorbar=False)
        
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
plt.savefig(os.path.join("composite_nao_ea_mh_t2m.svg", figname), format= "svg", bbox_inches="tight", dpi=300)




apply_style(fontsize=22, style=None, linewidth=2)
projection = ccrs.PlateCarree()
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                 figsize=(28, 25), subplot_kw={"projection": projection})

axes_op = [ax2, ax5, ax8]
axes_eq = [ax1, ax4, ax7]
axes_diff = [ax3, ax6, ax9]

for i,label in enumerate(labels_mh):
    if i == 0:
        
        plot_correlation(variable="Spearman Coefficients", data=mh_data[i].get("prec_eq_sval"), units="-", vmax=1,
                         vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                         level_ticks=7, ax=axes_eq[i], fig=fig, cbar_pos= [0.05, 0.01, 0.30, 0.02],
                         plot_pvalues=True, pvalue_data=mh_data[i].get("prec_eq_sig"), bottom_labels=True, 
                         title="NAO-Prec (EQ)-"+ label)
        
        plot_correlation(variable="Spearman Coefficients", data=mh_data[i].get("prec_op_sval"), units="-", vmax=1,
                         vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                         level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=mh_data[i].get("prec_op_sig"), bottom_labels=False,
                         left_labels=False, title="NAO-Prec (OP)-"+ label)
        
        plot_annual_mean(variable='Precipitation', data_alt=mh_data[i].get("prec_diff"), cmap=BrBG, units="mm/month",
             ax=axes_diff[i], fig=fig, vmax=50, vmin=-50, levels=22, domain="Europe Wide", level_ticks=9, 
             cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] Prec (OP-EQ)",
             left_labels=False, bottom_labels=False, add_colorbar=True)
    else:
        plot_correlation(variable="Spearman Coefficients", data=mh_data[i].get("prec_eq_sval"), units="-", vmax=1,
                         vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False,
                         level_ticks=7, ax=axes_eq[i], fig=fig,
                         plot_pvalues=True, pvalue_data=mh_data[i].get("prec_eq_sig"), bottom_labels=True, 
                         title="NAO-Prec (EQ)-"+ label)
        
        plot_correlation(variable="Spearman Coefficients", data=mh_data[i].get("prec_op_sval"), units="-", vmax=1,
                         vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                         level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=mh_data[i].get("prec_op_sig"), bottom_labels=False,
                         left_labels=False, title="NAO-Prec (OP)-"+ label)
        
        plot_annual_mean(variable='Precipitation', data_alt=mh_data[i].get("prec_diff"), cmap=BrBG, units="mm/month",
             ax=axes_diff[i], fig=fig, vmax=50, vmin=-50, levels=22, domain="Europe Wide", level_ticks=9, 
             cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] Prec (OP-EQ)",
             left_labels=False, bottom_labels=False, add_colorbar=False)
        
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
plt.savefig(os.path.join("composite_nao_ea_mh_prec.svg", figname), format= "svg", bbox_inches="tight", dpi=300)
    