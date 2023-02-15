# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 11:27:11 2023

@author: dboateng

1. Plot the annual means of the control experiments (and summer means if required)
2. Add the simulated anomalies for the different scenarios
3. 
"""
# import models
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
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, compute_lterm_diff


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"


# read data

CTL_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

# experiments
W2E1_PI_filename = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
W2E1_Mio278_filename = "a017_hpc-bw_e5w2.3_t159_MIO_W2E1_278ppm_t159l31.6h"
W2E1_Mio450_filename = "a016_hpc-bw_e5w2.3_t159_MIO_W2E1_450ppm_t159l31.6h"


years = "1003_1017"
period = "1m"

CTL_data, CTL_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=CTL_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_PI_data, W2E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_278_data, W2E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_450_data, W2E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)



def extract_relevant_vars(data, wiso, CTL_data=None, CTL_wiso=None, calc_diff=False):
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    u10 = extract_var(Dataset=data , varname="u10")
    v10 = extract_var(Dataset=data , varname="v10")
    
    # compute annual means
    
    if calc_diff:
        if CTL_data is not None:
            temp2_ctl = extract_var(Dataset=CTL_data, varname="temp2", units="°C")
            prec_ctl = extract_var(Dataset= CTL_data , varname="prec", units="mm/month")
            d18op_ctl = extract_var(Dataset=CTL_data , varname="d18op", units="per mil", Dataset_wiso= CTL_wiso)
            u10_ctl = extract_var(Dataset=CTL_data , varname="u10")
            v10_ctl = extract_var(Dataset=CTL_data , varname="v10")
            
        else:
            raise ValueError("The CTL data must be defined")
            
            
        temp2_alt = compute_lterm_diff(data_control=temp2_ctl, data_main=temp2, time="annual")
        prec_alt = compute_lterm_diff(data_control=prec_ctl, data_main=prec, time="annual")
        d18op_alt = compute_lterm_diff(data_control=d18op_ctl, data_main=d18op, time="annual")
        u10_alt = compute_lterm_diff(data_control=u10_ctl, data_main=u10, time="annual")
        v10_alt = compute_lterm_diff(data_control=v10_ctl, data_main=v10, time="annual")
        
    else:
        
        temp2_alt = compute_lterm_mean(data=temp2, time="annual")
        prec_alt = compute_lterm_mean(data=prec, time="annual")
        d18op_alt = compute_lterm_mean(data=d18op, time="annual")
        u10_alt = compute_lterm_mean(data=u10, time="annual")
        v10_alt = compute_lterm_mean(data=v10, time="annual")
   
    
    
    return temp2_alt, prec_alt, d18op_alt, u10_alt, v10_alt


# extract slm for miocene

slm_mio = W1E1_278_data.slm[0]
#xtract variables

# annual means
CTL_temp2, CTL_prec, CTL_d18op, CTL_u10, CTL_v10 = extract_relevant_vars(data=CTL_data, wiso=CTL_wiso)
W1E1_278_temp2, W1E1_278_prec, W1E1_278_d18op, W1E1_278_u10, W1E1_278_v10 = extract_relevant_vars(data=W1E1_278_data, wiso=W1E1_278_wiso)
W1E1_450_temp2, W1E1_450_prec, W1E1_450_d18op, W1E1_450_u10, W1E1_450_v10 = extract_relevant_vars(data=W1E1_450_data, wiso=W1E1_450_wiso)

# annual difference
W2E1_PI_temp2, W2E1_PI_prec, W2E1_PI_d18op, W2E1_PI_u10, W2E1_PI_v10 = extract_relevant_vars(data=W2E1_PI_data, wiso=W2E1_PI_wiso, 
                                                                                             CTL_data=CTL_data, CTL_wiso=CTL_wiso, calc_diff=True)
W2E1_278_temp2, W2E1_278_prec, W2E1_278_d18op, W2E1_278_u10, W2E1_278_v10 = extract_relevant_vars(data=W2E1_278_data, wiso=W2E1_278_wiso,
                                                                                                  CTL_data=W1E1_278_data, CTL_wiso=W1E1_278_wiso, calc_diff=True)
W2E1_450_temp2, W2E1_450_prec, W2E1_450_d18op, W2E1_450_u10, W2E1_450_v10 = extract_relevant_vars(data=W2E1_450_data, wiso=W2E1_450_wiso,
                                                                                                  CTL_data=W1E1_450_data, CTL_wiso=W1E1_450_wiso, calc_diff=True)


# plot
apply_style(fontsize=22, style=None, linewidth=2) 
    
projection = ccrs.PlateCarree()

def plot_prec():
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(28, 24), subplot_kw={"projection":  projection})
    
    plot_annual_mean(ax=ax1, variable="Precipitation", data_alt=CTL_prec, cmap=YlGnBu, units="mm/month", vmax=180, vmin=20, domain="Europe", 
                      levels=22, level_ticks=8, title="W1E1 (PI)", left_labels=True, bottom_labels=True, add_colorbar=True, cbar_pos = [0.15, 0.25, 0.25, 0.02],
                      plot_projection=projection, plot_winds=True, data_u=CTL_u10, data_v=CTL_v10,  orientation= "horizontal", fig=fig, plot_coastlines=True,
                      show_arrow_scale=True)
    
    plot_annual_mean(ax=ax2, variable="Precipitation", data_alt=W1E1_278_prec, cmap=YlGnBu, units="mm/month", vmax=180, vmin=20, domain="Europe", 
                      levels=22, level_ticks=8, title="W1E1 (MIO 278ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      plot_winds=True, data_u=W1E1_278_u10, data_v=W1E1_278_v10, fig=fig, plot_coastlines=False, 
                      sea_land_mask=slm_mio)
    plot_annual_mean(ax=ax3, variable="Precipitation", data_alt=W1E1_450_prec, cmap=YlGnBu, units="mm/month", vmax=180, vmin=20, domain="Europe", 
                      levels=22, level_ticks=8, title="W1E1 (MIO 450ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      plot_winds=True, data_u=W1E1_450_u10, data_v=W1E1_450_v10, fig=fig, plot_coastlines=False, 
                      sea_land_mask=slm_mio)
    
    
    #anomalies
    
    plot_annual_mean(ax=ax4, variable="Precipitation difference", data_alt=W2E1_PI_prec, cmap=BrBG, units="mm/month", vmax=80, vmin=-80, domain="Europe", 
                      levels=22, level_ticks=9, title="W2E1 - W1E1 (PI)", left_labels=True, bottom_labels=True, add_colorbar=True, cbar_pos = [0.45, 0.25, 0.25, 0.02],
                      plot_projection=projection, plot_winds=True, data_u=W2E1_PI_u10, data_v=W2E1_PI_v10,  orientation= "horizontal", fig=fig, plot_coastlines=True)
    
    plot_annual_mean(ax=ax5, variable="Precipitation difference", data_alt=W2E1_278_prec, cmap=BrBG, units="mm/month", vmax=80, vmin=-80, domain="Europe", 
                      levels=22, level_ticks=9, title="W2E1 - W1E1 (MIO 278ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      plot_winds=True, data_u=W2E1_278_u10, data_v=W2E1_278_v10, fig=fig, plot_coastlines=False, 
                      sea_land_mask=slm_mio)
    
    plot_annual_mean(ax=ax6, variable="Precipitation difference", data_alt=W2E1_450_prec, cmap=BrBG, units="mm/month", vmax=80, vmin=-80, domain="Europe", 
                      levels=22, level_ticks=9, title="W2E1 - W1E1 (MIO 450ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      plot_winds=True, data_u=W2E1_450_u10, data_v=W2E1_450_v10, fig=fig, plot_coastlines=False, 
                      sea_land_mask=slm_mio)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.20)
    plt.savefig(os.path.join(path_to_plots, "precipitation_CTL_mio_pi.svg"), format= "svg", bbox_inches="tight", dpi=300)


def plot_temp2():
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(28, 24), subplot_kw={"projection":  projection})
    
    plot_annual_mean(ax=ax1, variable="Temperature", data_alt=CTL_temp2, cmap=Spectral_r, units="°C", vmax=22, vmin=-5, domain="Europe", 
                      levels=22, level_ticks=11, title="CTL (PI)", left_labels=True, bottom_labels=True, add_colorbar=True, cbar_pos = [0.15, 0.25, 0.25, 0.02],
                      plot_projection=projection, orientation= "horizontal", fig=fig, plot_coastlines=True, center=False)
    
    plot_annual_mean(ax=ax2, variable="Temperature", data_alt=W1E1_278_temp2, cmap=Spectral_r, units="°C", vmax=22, vmin=-5, domain="Europe", 
                      levels=22, level_ticks=11, title="W1E1 (MIO 278ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                       fig=fig, plot_coastlines=False, sea_land_mask=slm_mio, center=False)
    plot_annual_mean(ax=ax3, variable="Temperature", data_alt=W1E1_450_temp2, cmap=Spectral_r, units="°C", vmax=22, vmin=-5, domain="Europe", 
                      levels=22, level_ticks=11, title="W1E1 (MIO 450ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      fig=fig, plot_coastlines=False, sea_land_mask=slm_mio, center=False)
    
    
    #anomalies
    
    plot_annual_mean(ax=ax4, variable="Temperature difference", data_alt=W2E1_PI_temp2, cmap=RdBu_r, units="°C", vmax=10, vmin=-10, domain="Europe", 
                      levels=22, level_ticks=11, title="W2E1 - W1E1 (PI)", left_labels=True, bottom_labels=True, add_colorbar=True, cbar_pos = [0.45, 0.25, 0.25, 0.02],
                      plot_projection=projection, orientation= "horizontal", fig=fig, plot_coastlines=True)
    
    plot_annual_mean(ax=ax5, variable="Temperature difference", data_alt=W2E1_278_temp2, cmap=RdBu_r, units="°C", vmax=10, vmin=-10, domain="Europe", 
                      levels=22, level_ticks=11, title="W2E1 - W1E1 (MIO 278ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      fig=fig, plot_coastlines=False, sea_land_mask=slm_mio)
    
    plot_annual_mean(ax=ax6, variable="Temperature difference", data_alt=W2E1_450_temp2, cmap=RdBu_r, units="°C", vmax=10, vmin=-10, domain="Europe", 
                      levels=22, level_ticks=11, title="W2E1 - W1E1 (MIO 450ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      fig=fig, plot_coastlines=False, sea_land_mask=slm_mio)
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.20)
    plt.savefig(os.path.join(path_to_plots, "temperature_CTL_mio_pi.svg"), format= "svg", bbox_inches="tight", dpi=300)


def plot_d18op():
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(28, 24), subplot_kw={"projection":  projection})
    
    plot_annual_mean(ax=ax1, variable="$\delta^{18}$Op vs SMOW", data_alt=CTL_d18op, cmap=RdYlBu, units="‰", vmax=2, vmin=-16, domain="Europe", 
                      levels=22, level_ticks=11, title="CTL (PI)", left_labels=True, bottom_labels=True, add_colorbar=True, cbar_pos = [0.15, 0.25, 0.25, 0.02],
                      plot_projection=projection, orientation= "horizontal", fig=fig, plot_coastlines=True, center=False)
    
    plot_annual_mean(ax=ax2, variable="$\delta^{18}$Op vs SMOW", data_alt=W1E1_278_d18op, cmap=RdYlBu, units="‰", vmax=2, vmin=-16, domain="Europe", 
                      levels=22, level_ticks=11, title="W1E1 (MIO 278ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                       fig=fig, plot_coastlines=False, sea_land_mask=slm_mio, center=False)
    plot_annual_mean(ax=ax3, variable="$\delta^{18}$Op vs SMOW", data_alt=W1E1_450_d18op, cmap=RdYlBu, units="‰", vmax=2, vmin=-16, domain="Europe", 
                      levels=22, level_ticks=11, title="W1E1 (MIO 450ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      fig=fig, plot_coastlines=False, sea_land_mask=slm_mio, center=False)
    
    
    
    #anomalies 
    
    
    plot_annual_mean(ax=ax4, variable="$\delta^{18}$Op vs SMOW difference", data_alt=W2E1_PI_d18op, cmap="PRGn", units="‰", vmax=8, vmin=-8, domain="Europe", 
                      levels=22, level_ticks=11, title="W2E1 - W1E1 (PI)", left_labels=True, bottom_labels=True, add_colorbar=True, cbar_pos = [0.45, 0.25, 0.25, 0.02],
                      plot_projection=projection, orientation= "horizontal", fig=fig, plot_coastlines=True)
    
    plot_annual_mean(ax=ax5, variable="$\delta^{18}$Op vs SMOW difference", data_alt=W2E1_278_d18op, cmap="PRGn", units="‰", vmax=8, vmin=-8, domain="Europe", 
                      levels=22, level_ticks=11, title="W2E1 - W1E1 (MIO 278ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      fig=fig, plot_coastlines=False, sea_land_mask=slm_mio)
    
    plot_annual_mean(ax=ax6, variable="$\delta^{18}$Op vs SMOW difference", data_alt=W2E1_450_d18op, cmap="PRGn", units="‰", vmax=8, vmin=-8, domain="Europe", 
                      levels=22, level_ticks=11, title="W2E1 - W1E1 (MIO 450ppm)", left_labels=True, bottom_labels=True, add_colorbar=False, plot_projection=projection, 
                      fig=fig, plot_coastlines=False, sea_land_mask=slm_mio)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.20)
    plt.savefig(os.path.join(path_to_plots, "d18op_CTL_mio_pi.svg"), format= "svg", bbox_inches="tight", dpi=300)


if __name__== "__main__":
    plot_prec()
    plot_temp2()
    plot_d18op()
    
