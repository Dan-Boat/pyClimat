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



main_path_mh = "D:/Datasets/CMIP6/PMIP/postprocessed/MH"
from path_to_data_mh import *

awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
cesm_path = os.path.join(main_path_mh, "CESM2")
ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")


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
    
    
    return nao_t2m_sval, nao_t2m_sig, ea_t2m_sval, ea_t2m_sig, nao_prec_sval, nao_prec_sig, ea_prec_sval, ea_prec_sig
    
    



awi_nao_t2m_sval, awi_nao_t2m_sig, awi_ea_t2m_sval, awi_ea_t2m_sig, awi_nao_prec_sval, awi_nao_prec_sig, awi_ea_prec_sval, awi_ea_prec_sig = perform_correlation_nao_ea_with_t2m_prec(pcs_path=AWI_DJF_PCS, 
                                                                                                path_to_data=awi_path, nao_mode=1, ea_mode=2, nao_factor=-1, ea_factor=-1)

cesm_nao_t2m_sval, cesm_nao_t2m_sig, cesm_ea_t2m_sval, cesm_ea_t2m_sig, cesm_nao_prec_sval, cesm_nao_prec_sig, cesm_ea_prec_sval, cesm_ea_prec_sig = perform_correlation_nao_ea_with_t2m_prec(pcs_path=CESM2_DJF_PCS, 
                                                                                                path_to_data=cesm_path, nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)

ec_earth_nao_t2m_sval, ec_earth_nao_t2m_sig, ec_earth_ea_t2m_sval, ec_earth_ea_t2m_sig, ec_earth_nao_prec_sval, ec_earth_nao_prec_sig, ec_earth_ea_prec_sval, ec_earth_ea_prec_sig = perform_correlation_nao_ea_with_t2m_prec(pcs_path=EC_Earth3_DJF_PCS, 
                                                                                                path_to_data=ec_earth_path, nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)


giss_nao_t2m_sval, giss_nao_t2m_sig, giss_ea_t2m_sval, giss_ea_t2m_sig, giss_nao_prec_sval, giss_nao_prec_sig, giss_ea_prec_sval, giss_ea_prec_sig = perform_correlation_nao_ea_with_t2m_prec(pcs_path=GISS_DJF_PCS, 
                                                                                                path_to_data=giss_path, nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)


ipsl_nao_t2m_sval, ipsl_nao_t2m_sig, ipsl_ea_t2m_sval, ipsl_ea_t2m_sig, ipsl_nao_prec_sval, ipsl_nao_prec_sig, ipsl_ea_prec_sval, ipsl_ea_prec_sig = perform_correlation_nao_ea_with_t2m_prec(pcs_path=IPSL_DJF_PCS, 
                                                                                                path_to_data=ipsl_path, nao_mode=1, ea_mode=2, nao_factor=-1, ea_factor=-1)


miroc_nao_t2m_sval, miroc_nao_t2m_sig, miroc_ea_t2m_sval, miroc_ea_t2m_sig, miroc_nao_prec_sval, miroc_nao_prec_sig, miroc_ea_prec_sval, miroc_ea_prec_sig = perform_correlation_nao_ea_with_t2m_prec(pcs_path=MIROC_DJF_PCS, 
                                                                                                path_to_data=miroc_path, nao_mode=1, ea_mode=2, nao_factor=-1, ea_factor=-1)

mpi_esm_nao_t2m_sval, mpi_esm_nao_t2m_sig, mpi_esm_ea_t2m_sval, mpi_esm_ea_t2m_sig, mpi_esm_nao_prec_sval, mpi_esm_nao_prec_sig, mpi_esm_ea_prec_sval, mpi_esm_ea_prec_sig = perform_correlation_nao_ea_with_t2m_prec(pcs_path=MPI_ESM_DJF_PCS, 
                                                                                                path_to_data=mpi_esm_path, nao_mode=1, ea_mode=2, nao_factor=-1, ea_factor=-1)


nao_t2m_data = [awi_nao_t2m_sval, cesm_nao_t2m_sval, ec_earth_nao_t2m_sval, giss_nao_t2m_sval, ipsl_nao_t2m_sval, miroc_nao_t2m_sval,
                mpi_esm_nao_t2m_sval]


nao_t2m_data_sig = [awi_nao_t2m_sig, cesm_nao_t2m_sig, ec_earth_nao_t2m_sig, giss_nao_t2m_sig, ipsl_nao_t2m_sig, miroc_nao_t2m_sig,
                mpi_esm_nao_t2m_sig]

nao_prec_data_sig = [awi_nao_prec_sig, cesm_nao_prec_sig, ec_earth_nao_prec_sig, giss_nao_prec_sig, ipsl_nao_prec_sig, miroc_nao_prec_sig,
                mpi_esm_nao_prec_sig]

nao_prec_data = [awi_nao_prec_sval, cesm_nao_prec_sval, ec_earth_nao_prec_sval, giss_nao_prec_sval, ipsl_nao_prec_sval, miroc_nao_prec_sval,
                mpi_esm_nao_prec_sval]



ea_t2m_data = [awi_ea_t2m_sval, cesm_ea_t2m_sval, ec_earth_ea_t2m_sval, giss_ea_t2m_sval, ipsl_ea_t2m_sval, miroc_ea_t2m_sval,
                mpi_esm_ea_t2m_sval]


ea_t2m_data_sig = [awi_ea_t2m_sig, cesm_ea_t2m_sig, ec_earth_ea_t2m_sig, giss_ea_t2m_sig, ipsl_ea_t2m_sig, miroc_ea_t2m_sig,
                mpi_esm_ea_t2m_sig]

ea_prec_data_sig = [awi_ea_prec_sig, cesm_ea_prec_sig, ec_earth_ea_prec_sig, giss_ea_prec_sig, ipsl_ea_prec_sig, miroc_ea_prec_sig,
                mpi_esm_ea_prec_sig]

ea_prec_data = [awi_ea_prec_sval, cesm_ea_prec_sval, ec_earth_ea_prec_sval, giss_ea_prec_sval, ipsl_ea_prec_sval, miroc_ea_prec_sval,
                mpi_esm_ea_prec_sval]





# generate plots 

def plot_climate_index_correlation(labels, data, figname, pval_data):
    
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols = 3, figsize=(28, 25), 
                                                                               subplot_kw={"projection": projection})

    axes = [ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    
    for i,label in enumerate(labels):
        if i == 0:

            plot_correlation(variable="Spearman Coefficients", data=data[i], units="-", vmax=1,
                             vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                             plot_pvalues=True, pvalue_data=pval_data[i], bottom_labels=True, title=label)
            
        else:
            

            plot_correlation(variable="Spearman Coefficients", data=data[i], units="-", vmax=1,
                             vmin=-1, cmap="PRGn", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, plot_pvalues=True, pvalue_data=pval_data[i],
                             title=label)

    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, figname), format= "svg", bbox_inches="tight", dpi=300)




labels_mh = ["AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
              "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"]
 

plot_climate_index_correlation(labels=labels_mh, data=nao_t2m_data, figname="corr_nao_t2m_mh_pmip.svg",
                               pval_data=nao_t2m_data_sig)

plot_climate_index_correlation(labels=labels_mh, data=ea_t2m_data, figname="corr_ea_t2m_mh_pmip.svg",
                               pval_data=ea_t2m_data_sig)

plot_climate_index_correlation(labels=labels_mh, data=nao_prec_data, figname="corr_nao_prec_mh_pmip.svg",
                               pval_data=nao_prec_data_sig)

plot_climate_index_correlation(labels=labels_mh, data=ea_prec_data, figname="corr_ea_prec_mh_pmip.svg",
                               pval_data=ea_prec_data_sig)

