# -*- coding: utf-8 -*-
"""
Created on Mon May 22 19:14:50 2023

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


from path_to_data_lm import *

main_path = "D:/Datasets/iGCM_datasets/"




def get_indices(pcs_path, nao_mode, ea_mode, nao_factor, ea_factor):
    
    
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


    
def perform_correlation_composite(atmos_index, path_to_data, pcs_path=None, filename=None):   
    
    if pcs_path is not None:
        df = pd.read_csv(pcs_path, parse_dates=["time"])
        
        
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
        
    
    
    t2m_season["time"] = df.time.values
    
    t2m_index = t2m_season.sel(time = t2m_season.time.isin(atmos_index.time.values))
    
    prec_season["time"] = df.time.values
    
    prec_index = prec_season.sel(time = prec_season.time.isin(atmos_index.time.values))
    
    d18O_season["time"] = df.time.values
    
    d18O_index = d18O_season.sel(time = d18O_season.time.isin(atmos_index.time.values))
    
    
    
    nao_t2m_sval, nao_t2m_pval, nao_t2m_sig = StatCorr(x=t2m_index, y=atmos_index, dim="time",
                                                       return_sig=True, sig=0.05)
    
    
    
    nao_prec_sval, nao_prec_pval, nao_prec_sig = StatCorr(x=prec_index, y=atmos_index, dim="time",
                                                       return_sig=True, sig=0.05)
    
    nao_d18O_sval, nao_d18O_pval, nao_d18O_sig = StatCorr(x=d18O_index, y=atmos_index, dim="time",
                                                       return_sig=True, sig=0.05)
    
    
    
    return_nao_reg = {"nao_reg_temp":nao_t2m_sval, "nao_reg_prec":nao_prec_sval, "nao_reg_d18O":nao_d18O_sval}
    
    return_nao_sig = {"nao_sig_temp":nao_t2m_sig, "nao_sig_prec":nao_prec_sig, "nao_sig_d18O":nao_d18O_sig}
    
    
    return_data = return_nao_reg | return_nao_sig
    
    return return_data


def cal_diff_OP_EQ(atmos_index_op, atmos_index_eq, path_to_data, pcs_path=None, filename=None):   
    if pcs_path is not None:
        df = pd.read_csv(pcs_path, parse_dates=["time"])
        

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
    
    
    t2m_season["time"] = df.time.values
    t2m_op = t2m_season.sel(time = t2m_season.time.isin(atmos_index_op.time.values)).mean(dim="time")
    t2m_eq = t2m_season.sel(time = t2m_season.time.isin(atmos_index_eq.time.values)).mean(dim="time")
    
    
    prec_season["time"] = df.time.values
    prec_op = prec_season.sel(time = prec_season.time.isin(atmos_index_op.time.values)).mean(dim="time")
    prec_eq = prec_season.sel(time = prec_season.time.isin(atmos_index_eq.time.values)).mean(dim="time")
    
    
    d18O_season["time"] = df.time.values
    d18O_op = d18O_season.sel(time = d18O_season.time.isin(atmos_index_op.time.values)).mean(dim="time")
    d18O_eq = d18O_season.sel(time = d18O_season.time.isin(atmos_index_eq.time.values)).mean(dim="time")
    
    t2m_diff = t2m_op - t2m_eq
    prec_diff = prec_op - prec_eq
    d18O_diff = d18O_op - d18O_eq
    
    return {"t2m_diff":t2m_diff.sortby("lon"), "prec_diff":prec_diff.sortby("lon"), "d18O_diff":d18O_diff.sortby("lon")}
    

    
def extract_all_for_model(pcs_path, path_to_data, nao_mode, ea_mode, 
                                             nao_factor, ea_factor, filename):
    
    eq, eq_pos, op = get_indices(pcs_path, nao_mode, ea_mode, 
                                                 nao_factor, ea_factor)
    data_eq = perform_correlation_composite(atmos_index=eq, path_to_data=path_to_data, pcs_path=pcs_path,
                                            filename=filename)
    data_op = perform_correlation_composite(atmos_index=op, path_to_data=path_to_data, pcs_path=pcs_path,
                                            filename=filename)
    data_diff = cal_diff_OP_EQ(atmos_index_op=op, atmos_index_eq=eq_pos, path_to_data=path_to_data, pcs_path=pcs_path,
                               filename=filename)
    
    
    return data_eq, data_op, data_diff



cesm_data_eq, cesm_data_op, cesm_data_diff = extract_all_for_model(pcs_path=CESM_DJF_PCs , path_to_data=main_path,
                                                    nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1, filename="CESM")


giss_data_eq, giss_data_op, giss_data_diff  = extract_all_for_model(pcs_path= GISS_DJF_PCs, path_to_data=main_path, nao_mode=1, 
                                                 ea_mode=3, nao_factor=-1, ea_factor=1, filename="GISS")

hadcm3_data_eq, hadcm3_data_op, hadcm3_data_diff  = extract_all_for_model(pcs_path=HADESM_DJF_PCs , path_to_data=main_path, nao_mode=1, 
                                                 ea_mode=3, nao_factor=-1, ea_factor=1, filename="HADCM3")




labels = ["iCESM", "GISS-E2-R", "iHadCM3"]

data_eq = [cesm_data_eq, giss_data_eq, hadcm3_data_eq]
data_op = [cesm_data_op, giss_data_op, hadcm3_data_op]
data_diff = [cesm_data_diff, giss_data_diff, hadcm3_data_diff]




def plot_nao_t2m():
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                     figsize=(22, 28), subplot_kw={"projection": projection})
    
    axes_op = [ax2, ax5, ax8,]
    axes_eq = [ax1, ax4, ax7,]
    axes_diff = [ax3, ax6, ax9]
    
    for i,label in enumerate(labels):
        if i == 0:
            
            plot_correlation(variable="Spearman Coefficients", data=data_eq[i].get("nao_reg_temp"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_eq[i], fig=fig, cbar_pos= [0.05, 0.01, 0.30, 0.02],
                             plot_pvalues=True, pvalue_data=data_eq[i].get("nao_sig_temp"), bottom_labels=True, 
                             title="NAO-t2m (EQ)-"+ label)
            
            plot_correlation(variable="Spearman Coefficients", data=data_op[i].get("nao_reg_temp"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=data_op[i].get("nao_sig_temp"), bottom_labels=True,
                             left_labels=True, title="NAO-t2m (OP)-"+ label)
            
            plot_annual_mean(variable='Temperature', data_alt=data_diff[i].get("t2m_diff"), cmap=RdBu, units="°C",
                 ax=axes_diff[i], fig=fig, vmax=4, vmin=-4, levels=22, domain="Europe Wide", level_ticks=11, 
                 cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] t2m (OP-EQ)",
                 left_labels=True, bottom_labels=True, label_format="%.1f", add_colorbar=True)
        else:
            plot_correlation(variable="Spearman Coefficients", data=data_eq[i].get("nao_reg_temp"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False,
                             level_ticks=7, ax=axes_eq[i], fig=fig,
                             plot_pvalues=True, pvalue_data=data_eq[i].get("nao_sig_temp"), bottom_labels=True, 
                             title="NAO-t2m (EQ)-"+ label)
            
            plot_correlation(variable="Spearman Coefficients", data=data_op[i].get("nao_reg_temp"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=data_op[i].get("nao_sig_temp"), bottom_labels=True,
                             left_labels=True, title="NAO-t2m (OP)-"+ label)
            
            plot_annual_mean(variable='Temperature', data_alt=data_diff[i].get("t2m_diff"), cmap=RdBu, units="°C",
                 ax=axes_diff[i], fig=fig, vmax=4, vmin=-4, levels=22, domain="Europe Wide", level_ticks=11, 
                 cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] t2m (OP-EQ)",
                 left_labels=True, bottom_labels=True, label_format="%.1f", add_colorbar=False)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, "composite_nao_ea_lm_t2m.svg"), format= "svg", bbox_inches="tight", dpi=300)


def plot_nao_prec():
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                     figsize=(22, 28), subplot_kw={"projection": projection})
    
    axes_op = [ax2, ax5, ax8]
    axes_eq = [ax1, ax4, ax7]
    axes_diff = [ax3, ax6, ax9]
    
    for i,label in enumerate(labels):
        if i == 0:
            
            plot_correlation(variable="Spearman Coefficients", data=data_eq[i].get("nao_reg_prec"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_eq[i], fig=fig, cbar_pos= [0.05, 0.01, 0.30, 0.02],
                             plot_pvalues=True, pvalue_data=data_eq[i].get("nao_sig_prec"), bottom_labels=True, 
                             title="NAO-Prec (EQ)-"+ label)
            
            plot_correlation(variable="Spearman Coefficients", data=data_op[i].get("nao_reg_prec"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=data_op[i].get("nao_sig_prec"), bottom_labels=True,
                             left_labels=True, title="NAO-Prec (OP)-"+ label)
            
            plot_annual_mean(variable='Precipitation', data_alt=data_diff[i].get("prec_diff"), cmap=BrBG, units="mm/month",
                 ax=axes_diff[i], fig=fig, vmax=50, vmin=-50, levels=22, domain="Europe Wide", level_ticks=9, 
                 cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] Prec (OP-EQ)",
                 left_labels=True, bottom_labels=True, add_colorbar=True)
        else:
            plot_correlation(variable="Spearman Coefficients", data=data_eq[i].get("nao_reg_prec"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False,
                             level_ticks=7, ax=axes_eq[i], fig=fig,
                             plot_pvalues=True, pvalue_data=data_eq[i].get("nao_sig_prec"), bottom_labels=True, 
                             title="NAO-Prec (EQ)-"+ label)
            
            plot_correlation(variable="Spearman Coefficients", data=data_op[i].get("nao_reg_prec"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=data_op[i].get("nao_sig_prec"), bottom_labels=True,
                             left_labels=True, title="NAO-Prec (OP)-"+ label)
            
            plot_annual_mean(variable='Precipitation', data_alt=data_diff[i].get("prec_diff"), cmap=BrBG, units="mm/month",
                 ax=axes_diff[i], fig=fig, vmax=50, vmin=-50, levels=22, domain="Europe Wide", level_ticks=9, 
                 cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] Prec (OP-EQ)",
                 left_labels=True, bottom_labels=True, add_colorbar=False)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10)
    plt.savefig(os.path.join(path_to_plots, "composite_nao_ea_lm_prec.svg"), format= "svg", bbox_inches="tight", dpi=300)
    
def plot_nao_d18O():
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols=3, 
                                                     figsize=(22, 28), subplot_kw={"projection": projection})
    
    axes_op = [ax2, ax5, ax8]
    axes_eq = [ax1, ax4, ax7]
    axes_diff = [ax3, ax6, ax9]
    
    for i,label in enumerate(labels):
        if i == 0:
            
            plot_correlation(variable="Spearman Coefficients", data=data_eq[i].get("nao_reg_d18O"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_eq[i], fig=fig, cbar_pos= [0.05, 0.01, 0.30, 0.02],
                             plot_pvalues=True, pvalue_data=data_eq[i].get("nao_sig_d18O"), bottom_labels=True, 
                             title="NAO-$\delta^{18}$Op (EQ)-"+ label)
            
            plot_correlation(variable="Spearman Coefficients", data=data_op[i].get("nao_reg_d18O"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=data_op[i].get("nao_sig_d18O"), bottom_labels=True,
                             left_labels=True, title="NAO-$\delta^{18}$Op (OP)-"+ label)
            
            plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=data_diff[i].get("d18O_diff"), cmap="PiYG", units="mm/month",
                 ax=axes_diff[i], fig=fig, vmax=2, vmin=-2, levels=22, domain="Europe Wide", level_ticks=11, 
                 cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] $\delta^{18}$Op (OP-EQ)",
                 left_labels=True, bottom_labels=True, add_colorbar=True, label_format="%.1f")
        else:
            plot_correlation(variable="Spearman Coefficients", data=data_eq[i].get("nao_reg_d18O"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False,
                             level_ticks=7, ax=axes_eq[i], fig=fig,
                             plot_pvalues=True, pvalue_data=data_eq[i].get("nao_sig_d18O"), bottom_labels=True, 
                             title="NAO-$\delta^{18}$Op (EQ)-"+ label)
            
            plot_correlation(variable="Spearman Coefficients", data=data_op[i].get("nao_reg_d18O"), units="-", vmax=0.8,
                             vmin=-0.8, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes_op[i], fig=fig, plot_pvalues=True, pvalue_data=data_op[i].get("nao_sig_d18O"), bottom_labels=True,
                             left_labels=True, title="NAO-$\delta^{18}$Op (OP)-"+ label)
            
            plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=data_diff[i].get("d18O_diff"), cmap="PiYG", units="mm/month",
                 ax=axes_diff[i], fig=fig, vmax=2, vmin=-2, levels=22, domain="Europe Wide", level_ticks=11, 
                 cbar_pos = [0.60, 0.01, 0.30, 0.02], title="[F] $\delta^{18}$Op (OP-EQ)",
                 left_labels=True, bottom_labels=True, add_colorbar=False, label_format="%.1f")
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10)
    plt.savefig(os.path.join(path_to_plots, "composite_nao_ea_lm_d18O.svg"), format= "svg", bbox_inches="tight", dpi=300)
  
if __name__ == "__main__":
    plot_nao_prec()      
    plot_nao_t2m()
    plot_nao_d18O()