# -*- coding: utf-8 -*-
"""
Created on Sun May 14 18:17:09 2023

@author: dboateng
"""
import os
import xarray as xr 
import numpy as np
import pandas as pd 


from pyClimat.analysis import compute_lterm_diff, compute_lterm_mean
from pyClimat.variables import extract_var

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pyClimat.data import read_from_path
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean

path_to_input_lgm = "D:/Datasets/ECHAM5/Inputs/T159_LGM"
path_to_input_mh = "D:/Datasets/ECHAM5/Inputs/T159_MH"
path_to_input_plio = "D:/Datasets/ECHAM5/Inputs/T159_PLIO"
path_to_input_pi = "D:/Datasets/ECHAM5/Inputs/T159_PI_Alps_100_east"

main_path = "D:/Datasets/Model_output_pst"

path_to_store = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

data_lgm = xr.open_dataset(os.path.join(path_to_input_lgm, "T159_VGRATCLIM_lgm_veg.nc"),
                           decode_cf=True, use_cftime=True)
data_mh = xr.open_dataset(os.path.join(path_to_input_mh, "T159_VGRATCLIM_MH_2.nc"),
                          decode_cf=True, use_cftime=True)
data_plio = xr.open_dataset(os.path.join(path_to_input_plio, "T159_VGRATCLIM_plio.nc"),
                            decode_cf=True, use_cftime=True)
data_pi = xr.open_dataset(os.path.join(path_to_input_pi, "T159_VGRATCLIM.nc"),
                          decode_cf=True, use_cftime=True)



lgm_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
plio_path = os.path.join(main_path, "PLIO", "MONTHLY_MEANS")
mh_path = os.path.join(main_path, "MH", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")
pd_path = os.path.join(main_path, )

filename_lterm = "1003_1017_1m_mlterm.nc"

PI_data = read_from_path(pi_path, filename_lterm)
LGM_data = read_from_path(lgm_path, filename_lterm)
PLIO_data = read_from_path(plio_path, filename_lterm)
MH_data = read_from_path(mh_path, filename_lterm)


def extract_veg(data):
    veg_data = extract_var(Dataset=data, varname="vegetation_frac")
    
    veg = compute_lterm_mean(data=veg_data, time="month", month="JJAS")
    
    return veg


# extract variables
def extract_fluxes_vars(data, data_pi):
    toa = extract_var(Dataset=data, varname="TOA")
    sw_land = extract_var(Dataset=data, varname="SW_land")
    
    toa_pi = extract_var(Dataset=data_pi, varname="TOA")
    sw_land_pi = extract_var(Dataset=data_pi, varname="SW_land")
    
    toa_diff = compute_lterm_diff(data_control=toa_pi, data_main=toa, time="month", month="JJAS")
    sw_land_diff = compute_lterm_diff(data_control=sw_land_pi, data_main=sw_land, time="month", month="JJAS")
    
    # toa_alt = compute_lterm_mean(data=sf, time="month", month="JJAS")
    # sw_land_alt = compute_lterm_mean(data=lf, time="month", month="JJAS")

    
        
    
    return toa_diff, sw_land_diff


toa_lgm, sw_land_lgm = extract_fluxes_vars(data_pi = PI_data, data=LGM_data)
toa_mh, sw_land_mh = extract_fluxes_vars(data_pi = PI_data, data=MH_data)
toa_plio, sw_land_plio = extract_fluxes_vars(data_pi = PI_data, data=PLIO_data)





def plot_vegetation():
    veg_pi = extract_veg(data_pi)
    veg_mh = extract_veg(data_mh)
    veg_lgm = extract_veg(data_lgm)
    veg_plio = extract_veg(data_plio)
    
    
    apply_style(fontsize=22, style=None, linewidth=2) 
        
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    fig, ((ax1,ax2),(ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(24, 18), 
                                        subplot_kw={"projection": projection})
    
    plot_annual_mean(variable="Vegetaion Ratio (fractional)", data_alt=veg_pi, cmap=vegetation, units="-", 
                     fig=fig, ax=ax1, vmax=0.9, vmin=0.1, levels=23, level_ticks=11, add_colorbar=False, time="JJAS",
                     bottom_labels=False, left_labels=False, title="[A] PI")
    
    plot_annual_mean(variable="Vegetaion Ratio (fractional)", data_alt=veg_mh, cmap=vegetation, units="-", 
                     fig=fig, ax=ax2, vmax=0.9, vmin=0.1, levels=23, level_ticks=11, add_colorbar=False, time="JJAS",
                     bottom_labels=False, left_labels=False, title="[B] MH")
    
    plot_annual_mean(variable="Vegetaion Ratio (fractional)", data_alt=veg_lgm, cmap=vegetation, units="-", 
                     fig=fig, ax=ax3, vmax=0.9, vmin=0.1, levels=23, level_ticks=11, add_colorbar=False, time="JJAS",
                     bottom_labels=False, left_labels=False, title="[C] LGM")
    
    plot_annual_mean(variable="Vegetaion Ratio (fractional)", data_alt=veg_plio, cmap=vegetation, units="-", 
                     fig=fig, ax=ax4, vmax=0.9, vmin=0.1, levels=23, level_ticks=11, time="JJAS",orientation= "horizontal",
                     bottom_labels=False, left_labels=False, add_colorbar=True, cbar_pos = [0.35, 0.01, 0.30, 0.03],
                     label_format="%.1f", title="[D] mPlio")
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.15)
    plt.savefig(os.path.join(path_to_store, "vegetation.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
plot_vegetation() 
