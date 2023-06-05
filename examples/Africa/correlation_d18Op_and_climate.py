# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 20:00:30 2023

@author: dboateng
"""

import os
import xarray as xr 
import numpy as np
import pandas as pd 


from pyClimat.analysis import extract_var
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pyClimat.data import read_from_path
from pyClimat.plot_utils import *
from pyClimat.plots import plot_correlation

from pyClimat.stats import  StatCorr
from pyClimat.utils import extract_region

main_path = "D:/Datasets/Model_output_pst"

path_to_store = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

lgm_path = os.path.join(main_path, "LGM", )
plio_path = os.path.join(main_path, "PLIO")
mh_path = os.path.join(main_path, "MH", )
pi_path = os.path.join(main_path, "PI",)


filename_lterm = "1003_1017_monthly.nc"

PI_data = read_from_path(pi_path, "PI_"+ filename_lterm, decode=True)
LGM_data = read_from_path(lgm_path, "LGM_"+ filename_lterm, decode=True)
PLIO_data = read_from_path(plio_path, "PLIO_"+ filename_lterm, decode=True)
MH_data = read_from_path(mh_path, "MH_"+ filename_lterm, decode=True)

filename_lterm = "1003_1017_monthly_wiso.nc"

PI_wiso = read_from_path(pi_path, "PI_"+ filename_lterm, decode=True)
LGM_wiso = read_from_path(lgm_path, "LGM_"+ filename_lterm, decode=True)
PLIO_wiso = read_from_path(plio_path, "PLIO_"+ filename_lterm, decode=True)
MH_wiso = read_from_path(mh_path, "MH_"+ filename_lterm, decode=True)

def perform_correlation(data, wiso, varname, units=None, lev=None, lev_units=None):
    d18Op = extract_var(Dataset=data, varname="d18op", units="per mil", Dataset_wiso=wiso,
                           )
    
    var = extract_var(Dataset=data, varname=varname, units=units, lev=lev,
                      lev_units=lev_units)
        
    d18Op_season = extract_region(data=d18Op, maxlon=40, minlon=-25, maxlat=35, minlat=-5, time="month", 
                                 month="JJAS")
    
    var_season = extract_region(data=var, maxlon=40, minlon=-25, maxlat=35, minlat=-5, time="month", 
                                 month="JJAS")
    
    sval, pval, sig = StatCorr(x=d18Op_season, y=var_season, dim="time",
                                                 return_sig=True, sig=0.05)
    
    return sval, sig

def perform_for_all(varname, units=None, lev=None, lev_units=None):
    
    pi_sval, pi_sig = perform_correlation(PI_data, PI_wiso, varname=varname, units=units, 
                                          lev=lev, lev_units=lev_units)
    
    mh_sval, mh_sig = perform_correlation(MH_data, MH_wiso, varname=varname, units=units, 
                                          lev=lev, lev_units=lev_units)
    
    lgm_sval, lgm_sig = perform_correlation(LGM_data, LGM_wiso, varname=varname, units=units, 
                                          lev=lev, lev_units=lev_units)
    
    plio_sval, plio_sig = perform_correlation(PLIO_data, PLIO_wiso, varname=varname, units=units, 
                                          lev=lev, lev_units=lev_units)
    
    return_data = {"pi_sval":pi_sval, "mh_sval":mh_sval, "lgm_sval":lgm_sval, "plio_sval":plio_sval, 
                   "pi_sig":pi_sig, "mh_sig":mh_sig,"lgm_sig":lgm_sig, "plio_sig":plio_sig}
    
    return return_data


prec_cor_d18Op = perform_for_all(varname="prec", units="mm/month")
temp_cor_d18Op = perform_for_all(varname="temp2", units="°C")



def plot_climate_correlation(data, figname):
    
    apply_style(fontsize=22, style=None, linewidth=2)
    projection = ccrs.PlateCarree()
    fig, (ax1,ax2, ax3, ax4) = plt.subplots(nrows = 1, ncols = 4, figsize=(28, 10), 
                                                                               subplot_kw={"projection": projection})

    axes = [ax1,ax2, ax3, ax4]
    sval = ["pi_sval", "mh_sval", "lgm_sval", "plio_sval"]
    sig = ["pi_sig", "mh_sig", "lgm_sig", "plio_sig"]
    labels = ["PI", "MH", "LGM", "PLIO"]
    
    for i,label in enumerate(labels):
        if i == 0:

            plot_correlation(variable="Spearman Coefficients", data=data.get(sval[i]), units="-", vmax=1,
                             vmin=-1, cmap="PRGn", domain="West Africa", levels=22,cbar=True, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                             plot_pvalues=True, pvalue_data=data.get(sig[i]), bottom_labels=True, title=label)
            
        else:
            

            plot_correlation(variable="Spearman Coefficients", data=data.get(sval[i]), units="-", vmax=1,
                             vmin=-1, cmap="PRGn", domain="West Africa", levels=22,cbar=False, cbar_orientation="horizontal",
                             level_ticks=7, ax=axes[i], fig=fig, plot_pvalues=True, pvalue_data=data.get(sig[i]),
                             title=label)

    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10, hspace=0.01)
    plt.savefig(os.path.join(path_to_store, figname), format= "svg", bbox_inches="tight", dpi=600)
    
    
plot_climate_correlation(data=prec_cor_d18Op, figname="corr_with_d18Op_prec_JJAS.svg")
plot_climate_correlation(data=temp_cor_d18Op, figname="corr_with_d18Op_temp_JJAS.svg")
    
