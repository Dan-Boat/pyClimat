# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:55:20 2023

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
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var



path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

# experiments
W2E1_Mio278_filename = "a017_hpc-bw_e5w2.3_t159_MIO_W2E1_278ppm_t159l31.6h"
W2E1_Mio450_filename = "a016_hpc-bw_e5w2.3_t159_MIO_W2E1_450ppm_t159l31.6h"

W2E15_Mio278_filename = "a023_dkrz-levante_e5w2.3_t159_MIO_W2E1.5_278ppm_t159l31.6h"
W2E15_Mio450_filename = "a022_hpc-bw_e5w2.3_t159_MIO_W2E1.5_450ppm_t159l31.6h"


W2E0_Mio278_filename="a019_hpc-bw_e5w2.3_t159_MIO_W2E0_278ppm_t159l31.6h"
W2E0_Mio450_filename="a018_hpc-bw_e5w2.3_t159_MIO_W2E0_450ppm_t159l31.6h"

W2E2_Mio278_filename="a020_dkrz-levante_e5w2.3_t159_MIO_W2E2_278ppm_t159l31.6h"
W2E2_Mio450_filename="a021_dkrz-levante_e5w2.3_t159_MIO_W2E2_450ppm_t159l31.6h"


# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E1_278_data, W2E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_450_data, W2E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)

W2E15_278_data, W2E15_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E15_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E15_450_data, W2E15_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E15_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E0_278_data, W2E0_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_Mio278_filename, 
                                                    years=years, period=period, read_wiso=True)

W2E0_450_data, W2E0_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E2_278_data, W2E2_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio278_filename, 
                                                    years=years, period=period, read_wiso=True)

W2E2_450_data, W2E2_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


def extract_prec_and_analysis(exp_data, exp_wiso, W1E1_data, W1E1_wiso):
    
    # extract prec
    exp_prec_data = extract_var(Dataset=exp_data, varname="prec", units="mm/month")
    w1e1_prec_data = extract_var(Dataset=W1E1_data, varname="prec", units="mm/month")
    
    simulated_prec_change_alt = compute_lterm_diff(data_control=w1e1_prec_data, data_main=exp_prec_data,
                                          time="annual")
    simulated_climatologies = compute_lterm_mean(data=exp_prec_data, time="annual")
    
    
    return_data_monthly = {"control_mon": w1e1_prec_data, "simulated_mon": exp_prec_data,}
    return_data_ltm = {"simulated_change_ltm": simulated_prec_change_alt, "climatology_means": simulated_climatologies}
    
    return_data = return_data_monthly | return_data_ltm
    
    
    return return_data


W2E0_278_diff = extract_prec_and_analysis(exp_data=W2E0_278_data, exp_wiso=W2E0_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)


W2E0_450_diff = extract_prec_and_analysis(exp_data=W2E0_450_data, exp_wiso=W2E0_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso) 
      
W2E1_278_diff = extract_prec_and_analysis(exp_data=W2E1_278_data, exp_wiso=W2E1_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)

W2E1_450_diff = extract_prec_and_analysis(exp_data=W2E1_450_data, exp_wiso=W2E1_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso)

W2E15_278_diff = extract_prec_and_analysis(exp_data=W2E15_278_data, exp_wiso=W2E15_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)

W2E15_450_diff = extract_prec_and_analysis(exp_data=W2E15_450_data, exp_wiso=W2E15_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso)

W2E2_278_diff = extract_prec_and_analysis(exp_data=W2E1_278_data, exp_wiso=W2E1_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)

W2E2_450_diff = extract_prec_and_analysis(exp_data=W2E1_450_data, exp_wiso=W2E1_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso)

#plot the anomalies after correcting the lapse rate

def plot_prec_simulated():
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    fig,((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, figsize=(22,28), 
                                                                    subplot_kw={"projection":projection})
    
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    labels = ["(a) W2E0 - W1E1 (MIO 278ppm)",  "(b) W2E0 - W1E1 (MIO 450ppm)",  "(c) W2E1 - W1E1 (MIO 278ppm)", 
              "(d) W2E1 - W1E1 (MIO 450ppm)", "(e) W2E1.5 - W1E1 (MIO 278ppm)",  "(f) W2E1.5 - W1E1 (MIO 450ppm)",
              "(g) W2E2 - W1E1 (MIO 278ppm)",  "(h) W2E2 - W1E1 (MIO 450ppm)"]
    
    data = [W2E0_278_diff, W2E0_450_diff, W2E1_278_diff, W2E1_450_diff, W2E15_278_diff, W2E15_450_diff,
            W2E2_278_diff, W2E2_450_diff,]
    
    
    for i,label in enumerate(labels):
        if i == 0:
            
            plot_annual_mean(variable="Precipitation difference", data_alt=data[i].get("simulated_change_ltm"), ax=axes[i],
                             cmap=PrecAno, units="mm/month", vmax=60, vmin=-60, 
                            levels=22, level_ticks=9, add_colorbar=True, cbar_pos= [0.30, 0.05, 0.45, 0.02], 
                            orientation="horizontal", plot_coastlines=False, bottom_labels=True,
                            left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", compare_data1=data[i].get("simulated_mon"), 
                            compare_data2=data[i].get("control_mon"), max_pvalue=0.1, plot_stats=True,
                            hatches=".", title=label)
            
        else:
            plot_annual_mean(variable="Precipitation difference", data_alt=data[i].get("simulated_change_ltm"), ax=axes[i],
                             cmap=PrecAno, units="mm/month", vmax=60, vmin=-60, 
                            levels=22, level_ticks=9, add_colorbar=False, plot_coastlines=False, bottom_labels=True,
                            left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", compare_data1=data[i].get("simulated_mon"), 
                            compare_data2=data[i].get("control_mon"), max_pvalue=0.1, plot_stats=True,
                            hatches=".", title=label)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    plt.savefig(os.path.join(path_to_plots, "prec_change.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
def plot_prec_climatologies():
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    fig,((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, figsize=(22,28), 
                                                                    subplot_kw={"projection":projection})
    
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    labels = ["(a) W2E0 (MIO 278ppm)",  "(b) W2E0 (MIO 450ppm)",  "(c) W2E1 (MIO 278ppm)", 
              "(d) W2E1 (MIO 450ppm)", "(e) W2E1.5 (MIO 278ppm)",  "(f) W2E1.5 (MIO 450ppm)",
              "(g) W2E2 (MIO 278ppm)",  "(h) W2E2 (MIO 450ppm)"]
    
    data = [W2E0_278_diff, W2E0_450_diff, W2E1_278_diff, W2E1_450_diff, W2E15_278_diff, W2E15_450_diff,
            W2E2_278_diff, W2E2_450_diff,]
    
    
    for i,label in enumerate(labels):
        if i == 0:
            
            plot_annual_mean(variable="Precipitation", data_alt=data[i].get("climatology_means"), ax=axes[i],
                             cmap=YlGnBu, units="mm/month", vmax=300, vmin=0, 
                            levels=22, level_ticks=7, add_colorbar=True, cbar_pos= [0.30, 0.05, 0.45, 0.02], 
                            orientation="horizontal", plot_coastlines=False, bottom_labels=True,
                            left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", title=label)
            
        else:
            plot_annual_mean(variable="Precipitation", data_alt=data[i].get("climatology_means"), ax=axes[i],
                             cmap=YlGnBu, units="mm/month", vmax=300, vmin=0, 
                            levels=22, level_ticks=7, add_colorbar=False, plot_coastlines=False, bottom_labels=True,
                            left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", title=label)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    plt.savefig(os.path.join(path_to_plots, "prec_climatologies.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
#plot_prec_simulated()
plot_prec_climatologies()