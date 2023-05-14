# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 20:10:57 2023

@author: dboateng
1. compare correlation with echam d18Op and dD and NAO and EA index, 
extract regions to check distribution of correlation, and perform correlation with
regional means. Perferfom the analysis with ERA and GNIP datasets
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
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_var, extract_transect
from pyClimat.utils import extract_region


from paths_to_data import *
echam_pd_data_path = "D:/Datasets/Model_output_pst/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021"
        
        



def perform_correlation(atmos_index, varname, units=None, sig=0.05, lev=None, lev_units=None, save_sval=False,
                        path_to_save=main_path_to_data, filename=None, season="DJF"):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data = extract_var(Dataset=echam_data, varname=varname, units=units, Dataset_wiso=echam_wiso,
                       lev=lev, lev_units=lev_units)
    
    data_season = extract_region(data=data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season=season)
    
    sval, pval, sig = StatCorr(x=data_season, y=atmos_index, dim="time",
                                                 return_sig=True, sig=0.05)
    
    if save_sval:
        results = xr.merge([sval, sig])
        results.to_netcdf(os.path.join(path_to_save, filename))
        
    return sval, pval, sig
    


from pyClimat.plots import plot_correlation

# read the save data

def plot_correlation_DJF():
    # read all the required datsets (for winter)
    df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

    #select the node ands convert to xarray

    nao_index_echam = xr.DataArray(df_echam_pcs[str(1)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})

    ea_index_echam = xr.DataArray(df_echam_pcs[str(2)], dims="time", coords={"time": df_echam_pcs["time"]})
    
    
    d18op_nao_sval, d18op_nao_pval, d18op_nao_sig = perform_correlation(atmos_index=nao_index_echam, 
                                                                        varname="d18op", units="per mil",
                                                                        filename="d18op_nao_scc.nc", 
                                                                        season="DJF")
    d18op_ea_sval, d18op_ea_pval, d18op_ea_sig = perform_correlation(atmos_index=ea_index_echam, 
                                                                        varname="d18op", units="per mil",
                                                                        filename="d18op_ea_scc.nc", 
                                                                        season="DJF")


    t2m_nao_sval, t2m_nao_pval, t2m_nao_sig = perform_correlation(atmos_index=nao_index_echam, 
                                                                        varname="temp2", units="°C",
                                                                        filename="t2m_nao_scc.nc", 
                                                                        season="DJF")
    t2m_ea_sval, t2m_ea_pval, t2m_ea_sig = perform_correlation(atmos_index=ea_index_echam, 
                                                                        varname="temp2", units="°C",
                                                                        filename="t2m_ea_scc.nc", 
                                                                        season="DJF")

    prec_nao_sval, prec_nao_pval, prec_nao_sig = perform_correlation(atmos_index=nao_index_echam, 
                                                                        varname="prec", units="mm/month",
                                                                        filename="prec_nao_scc.nc", 
                                                                        season="DJF")
    prec_ea_sval, prec_ea_pval, prec_ea_sig = perform_correlation(atmos_index=ea_index_echam, 
                                                                        varname="prec", units="mm/month",
                                                                        filename="prec_ea_scc.nc", 
                                                                        season="DJF")
    
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols=3, 
                                                     figsize=(28, 22), subplot_kw={"projection": projection})
    
    plot_correlation(variable="Spearman Coefficients", data=d18op_nao_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax1, fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                     plot_pvalues=True, pvalue_data=d18op_nao_sig, bottom_labels=False, title="[A] NAO-$\delta^{18}$Op")
    
    plot_correlation(variable="Spearman Coefficients", data=d18op_ea_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax4, fig=fig, plot_pvalues=True, pvalue_data=d18op_ea_sig,
                     title="[D] EA-$\delta^{18}$Op")
    
    plot_correlation(variable="Spearman Coefficients", data=t2m_nao_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax2, fig=fig, plot_pvalues=True, pvalue_data=t2m_nao_sig, bottom_labels=False,
                     left_labels=False, title="[B] NAO-t2m")
    
    plot_correlation(variable="Spearman Coefficients", data=t2m_ea_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax5, fig=fig, plot_pvalues=True, pvalue_data=t2m_ea_sig, left_labels=False,
                     title="[E] EA-t2m")
    
    plot_correlation(variable="Spearman Coefficients", data=prec_nao_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax3, fig=fig, plot_pvalues=True, pvalue_data=prec_nao_sig, bottom_labels=False,
                     left_labels=False, title="[C] NAO-prec")
    
    plot_correlation(variable="Spearman Coefficients", data=prec_ea_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax6, fig=fig, plot_pvalues=True, pvalue_data=prec_ea_sig,
                     left_labels=False, title="[F] EA-prec")
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.001)
    plt.savefig(os.path.join(path_to_plots, "corr_nao_ea_d18op_djf.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.show()    
    
    
def plot_correlation_JJA():
    
    #summer 
    # read all the required datsets (for winter)
    df_echam_pcs_jja = pd.read_csv(PCs_ECHAM_JJA_path, parse_dates=["time"])

    #select the node ands convert to xarray

    nao_index_echam_jja = xr.DataArray(df_echam_pcs_jja[str(1)], dims="time", coords={"time": df_echam_pcs_jja["time"]})

    ea_index_echam_jja = xr.DataArray(df_echam_pcs_jja[str(2)] *-1, dims="time", coords={"time": df_echam_pcs_jja["time"]})
    
    d18op_nao_sval, d18op_nao_pval, d18op_nao_sig = perform_correlation(atmos_index=nao_index_echam_jja, 
                                                                        varname="d18op", units="per mil",
                                                                        filename="d18op_nao_scc.nc", season="JJA")
    d18op_ea_sval, d18op_ea_pval, d18op_ea_sig = perform_correlation(atmos_index=ea_index_echam_jja, 
                                                                        varname="d18op", units="per mil",
                                                                        filename="d18op_ea_scc.nc", season="JJA")


    t2m_nao_sval, t2m_nao_pval, t2m_nao_sig = perform_correlation(atmos_index=nao_index_echam_jja, 
                                                                        varname="temp2", units="°C",
                                                                        filename="t2m_nao_scc.nc", season="JJA")
    t2m_ea_sval, t2m_ea_pval, t2m_ea_sig = perform_correlation(atmos_index=ea_index_echam_jja, 
                                                                        varname="temp2", units="°C",
                                                                        filename="t2m_ea_scc.nc", season="JJA")

    prec_nao_sval, prec_nao_pval, prec_nao_sig = perform_correlation(atmos_index=nao_index_echam_jja, 
                                                                        varname="prec", units="mm/month",
                                                                        filename="prec_nao_scc.nc", season="JJA")
    prec_ea_sval, prec_ea_pval, prec_ea_sig = perform_correlation(atmos_index=ea_index_echam_jja, 
                                                                        varname="prec", units="mm/month",
                                                                        filename="prec_ea_scc.nc", season="JJA")
    
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols=3, 
                                                     figsize=(28, 22), subplot_kw={"projection": projection})
    
    plot_correlation(variable="Spearman Coefficients", data=d18op_nao_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax1, fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                     plot_pvalues=True, pvalue_data=d18op_nao_sig, bottom_labels=False, title="[A] NAO-$\delta^{18}$Op")
    
    plot_correlation(variable="Spearman Coefficients", data=d18op_ea_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax4, fig=fig, plot_pvalues=True, pvalue_data=d18op_ea_sig,
                     title="[D] EA-$\delta^{18}$Op")
    
    plot_correlation(variable="Spearman Coefficients", data=t2m_nao_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax2, fig=fig, plot_pvalues=True, pvalue_data=t2m_nao_sig, bottom_labels=False,
                     left_labels=False, title="[B] NAO-t2m")
    
    plot_correlation(variable="Spearman Coefficients", data=t2m_ea_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax5, fig=fig, plot_pvalues=True, pvalue_data=t2m_ea_sig, left_labels=False,
                     title="[E] EA-t2m")
    
    plot_correlation(variable="Spearman Coefficients", data=prec_nao_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax3, fig=fig, plot_pvalues=True, pvalue_data=prec_nao_sig, bottom_labels=False,
                     left_labels=False, title="[C] NAO-prec")
    
    plot_correlation(variable="Spearman Coefficients", data=prec_ea_sval, units="-", vmax=1,
                     vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax6, fig=fig, plot_pvalues=True, pvalue_data=prec_ea_sig,
                     left_labels=False, title="[F] EA-prec")
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.001)
    plt.savefig(os.path.join(path_to_plots, "corr_nao_ea_d18op_jja.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.show()

if __name__ == "__main__":
    plot_correlation_JJA()
    plot_correlation_DJF()