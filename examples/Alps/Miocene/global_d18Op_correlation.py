# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 11:10:04 2023

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_correlation
from pyClimat.data import read_ECHAM_processed, read_from_path
from pyClimat.analysis import extract_var, compute_lterm_mean, compute_lterm_diff
from pyClimat.stats import  StatCorr
from pyClimat.utils import extract_region


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_path = "D:/Datasets/Model_output_pst/a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h/output_processed/"
W1E1_450_path = "D:/Datasets/Model_output_pst/a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h/output_processed/"
W1E1_PI_path = "D:/Datasets/Model_output_pst/a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h/output_processed/"

# reading data 
# read data (monthly)

# main data
filename_lterm = "1003_1017_monthly.nc"

W1E1_278_data = read_from_path(os.path.join(W1E1_278_path, "MONTHLY_MEANS"), filename_lterm, decode=True)
W1E1_450_data = read_from_path(os.path.join(W1E1_450_path, "MONTHLY_MEANS"), filename_lterm, decode=True)
W1E1_PI_data = read_from_path(os.path.join(W1E1_PI_path, "MONTHLY_MEANS"), filename_lterm, decode=True)


#wiso data 
filename_lterm_wiso = "1003_1017_monthly_wiso.nc"
W1E1_278_wiso = read_from_path(os.path.join(W1E1_278_path, "MONTHLY_MEANS_WISO"), filename_lterm_wiso, decode=True)
W1E1_450_wiso = read_from_path(os.path.join(W1E1_450_path, "MONTHLY_MEANS_WISO"), filename_lterm_wiso, decode=True)
W1E1_PI_wiso = read_from_path(os.path.join(W1E1_PI_path, "MONTHLY_MEANS_WISO"), filename_lterm_wiso, decode=True)




def perform_correlation(data, wiso, varname, units=None, lev=None, lev_units=None):
    
    d18Op = extract_var(Dataset=data, varname="d18op", units="per mil", Dataset_wiso=wiso,
                           )
    
    var = extract_var(Dataset=data, varname=varname, units=units, lev=lev,
                      lev_units=lev_units)
    
    d18Op_season = extract_region(data=d18Op, maxlon=48, minlon=-48, maxlat=80, minlat=30, time="annual", 
                                 )
    
    var_season = extract_region(data=var, maxlon=48, minlon=-48, maxlat=80, minlat=30, time="annual", 
                                )
        
    
    sval, pval, sig = StatCorr(x=d18Op_season, y=var_season, dim="time",
                                                 return_sig=True, sig=0.05)
    
    return sval, sig

def perform_for_all(varname, units=None, lev=None, lev_units=None):
    
    pi_sval, pi_sig = perform_correlation(W1E1_PI_data, W1E1_PI_wiso, varname=varname, units=units, 
                                          lev=lev, lev_units=lev_units)
    
    mio278_sval, mio278_sig = perform_correlation(W1E1_278_data, W1E1_278_wiso, varname=varname, units=units, 
                                          lev=lev, lev_units=lev_units)
    
    mio450_sval, mio450_sig = perform_correlation(W1E1_450_data, W1E1_450_wiso, varname=varname, units=units, 
                                          lev=lev, lev_units=lev_units)
    
    return_data = {"pi_sval":pi_sval, "mio278_sval":mio278_sval, "mio450_sval":mio450_sval, 
                   "pi_sig":pi_sig, "mio278_sig":mio278_sig,"mio450_sig":mio450_sig,}
    
    return return_data

prec_cor_d18Op = perform_for_all(varname="prec", units="mm/month")
temp_cor_d18Op = perform_for_all(varname="temp2", units="°C")


apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)
#projection = ccrs.PlateCarree()

fig,((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(24,15), 
                                                      subplot_kw={"projection":projection})


axes = [ax1, ax2, ax3, ax4, ax5, ax6]

sval = ["pi_sval", "mio278_sval", "mio450_sval", "pi_sval", "mio278_sval", "mio450_sval"]

sig = ["pi_sig", "mio278_sig", "mio450_sig", "pi_sig", "mio278_sig", "mio450_sig"]

labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

for i,label in enumerate(labels):
    
    if i > 2:
        data = temp_cor_d18Op
        
    else:
        data = prec_cor_d18Op
    
    if i == 0:

        plot_correlation(variable="Spearman Coefficients", data=data.get(sval[i]), units="-", vmax=1,
                         vmin=-1, cmap="coolwarm", domain="Europe Wide", levels=22,cbar=True, cbar_orientation="horizontal",
                         level_ticks=7, ax=axes[i], fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                         plot_pvalues=True, pvalue_data=data.get(sig[i]), bottom_labels=True, title=label,
                         plot_projection=projection, plot_coastlines=False, sea_land_mask=mio_slm)
        
    else:
        

        plot_correlation(variable="Spearman Coefficients", data=data.get(sval[i]), units="-", vmax=1,
                         vmin=-1, cmap="coolwarm", domain="Europe Wide", levels=22,cbar=False, cbar_orientation="horizontal",
                         level_ticks=7, ax=axes[i], fig=fig, plot_pvalues=True, pvalue_data=data.get(sig[i]),
                         title=label, plot_projection=projection, plot_coastlines=False, sea_land_mask=mio_slm)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, "d18Op_correlation_with_prec_temp.svg"), format= "svg", bbox_inches="tight", dpi=600)