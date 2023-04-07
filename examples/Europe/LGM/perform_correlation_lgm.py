# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:19:30 2023

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



main_path_lgm = "D:/Datasets/CMIP6/PMIP/postprocessed/LGM"
from path_to_data_lgm import *

awi_path = os.path.join(main_path_lgm, "AWI-ESM-1-1-LR")
cesm_waccm_path = os.path.join(main_path_lgm, "CESM2-WACCM-FV2")
inm_cm_path = os.path.join(main_path_lgm, "INM-CM4-8")
miroc_path = os.path.join(main_path_lgm, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_lgm, "MPI-ESM1-2-LR")


def perform_correlation_nao_ea_with_t2m_prec(pcs_path, path_to_data, nao_mode, ea_mode, 
                                             nao_factor=1, ea_factor=1):
    
    
    df = pd.read_csv(pcs_path, parse_dates=["time"])
    
    NAO_indices = xr.DataArray(df[str(nao_mode)] * nao_factor, dims="time", coords={"time": df["time"]})
    EA_indices = xr.DataArray(df[str(ea_mode)] * ea_factor, dims="time", coords={"time": df["time"]})
    

    t2m_data = read_from_path(path=path_to_data, filename="tas_monthly.nc", decode=True,varname="tas",
                              ) - 273.15 #°C
    prec_data = read_from_path(path=path_to_data, filename="pr_monthly.nc", decode=True, 
                               varname="pr",) *60*60*24*30  #mm/month
    
    t2m_season = extract_region(data=t2m_data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF") 
    
    prec_season = extract_region(data=prec_data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF") 
    
    nao_t2m_sval, nao_t2m_pval, nao_t2m_sig = StatCorr(x=t2m_season, y=NAO_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    ea_t2m_sval, ea_t2m_pval, ea_t2m_sig = StatCorr(x=t2m_season, y=EA_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    nao_prec_sval, nao_prec_pval, nao_prec_sig = StatCorr(x=prec_season, y=NAO_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    ea_prec_sval, ea_prec_pval, ea_prec_sig = StatCorr(x=prec_season, y=EA_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    
    results = {"nao_t2m_sval":nao_t2m_sval, "nao_t2m_sig":nao_t2m_sig, "ea_t2m_sval":ea_t2m_sval, 
               "ea_t2m_sig":ea_t2m_sig, "nao_prec_sval": nao_prec_sval, "nao_prec_sig": nao_prec_sig,
               "ea_prec_sval":ea_prec_sval, "ea_prec_sig":ea_prec_sig}
    
    return results
    
    



awi_data = perform_correlation_nao_ea_with_t2m_prec(pcs_path=AWI_DJF_PCS, path_to_data=awi_path,
                                                    nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)

cesm_data = perform_correlation_nao_ea_with_t2m_prec(pcs_path=CESM_WA_DJF_PCS, path_to_data=cesm_waccm_path, 
                                                     nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)


inm_data = perform_correlation_nao_ea_with_t2m_prec(pcs_path=INM_DJF_PCS, path_to_data=inm_cm_path, 
                                                     nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)



miroc_data = perform_correlation_nao_ea_with_t2m_prec(pcs_path=MIROC_DJF_PCS, path_to_data=miroc_path,
                                                    nao_mode=1, ea_mode=4, nao_factor=-1, ea_factor=1)

mpi_data = perform_correlation_nao_ea_with_t2m_prec(pcs_path=MPI_ESM_DJF_PCS, path_to_data=mpi_esm_path, 
                                                    nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)



# generate plots 

def plot_climate_index_correlation(corr_name, corr_sig_name, figname):
    
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(25, 22), 
                                                                               subplot_kw={"projection": projection})
    labels_lgm = ["AWI-ESM-1-1-LR", "CESM2-WACCM-FV2", "INM-CM4-8",
                  "MIROC-ES2L", "MPI-ESM1-2-LR"]
    
    lgm_data = [awi_data, cesm_data, inm_data, miroc_data, mpi_data]
    
    axes = [ax1,ax2, ax3, ax4, ax5, ax6]
    
    for i,label in enumerate(labels_lgm):
        if i == 0:

            plot_correlation(variable="Spearman Coefficients", data=lgm_data[i].get(corr_name), units="-", vmax=1,
                             vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                             plot_pvalues=True, pvalue_data=lgm_data[i].get(corr_sig_name), bottom_labels=True, title=label)
            
        else:
            

            plot_correlation(variable="Spearman Coefficients", data=lgm_data[i].get(corr_name), units="-", vmax=1,
                             vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, plot_pvalues=True,
                             pvalue_data=lgm_data[i].get(corr_sig_name),
                             title=label)

    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, figname), format= "svg", bbox_inches="tight", dpi=300)

 

plot_climate_index_correlation(corr_name="nao_t2m_sval", corr_sig_name="nao_t2m_sig", figname="corr_nao_t2m_lgm_pmip.svg",)

plot_climate_index_correlation(corr_name="ea_t2m_sval", corr_sig_name="ea_t2m_sig", figname="corr_ea_t2m_lgm_pmip.svg",)

plot_climate_index_correlation(corr_name="nao_prec_sval", corr_sig_name="nao_prec_sig", figname="corr_nao_prec_lgm_pmip.svg")

plot_climate_index_correlation(corr_name="ea_prec_sval", corr_sig_name="ea_prec_sig", figname="corr_ea_prec_lgm_pmip.svg")

