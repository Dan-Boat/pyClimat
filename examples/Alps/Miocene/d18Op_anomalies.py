# -*- coding: utf-8 -*-
"""
Created on Tue May 16 17:31:00 2023

@author: dboateng
This script plots the difference between the expected isotopic difference and the actual simulation difference.
The aim is to highlight what the proxies might miss using present-day lapse rate for reconstructing topography
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
from pyClimat.analysis import extract_var, compute_lterm_mean, compute_lterm_diff


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


def extract_d18Op_elev_and_analysis(exp_data, exp_wiso, W1E1_data, W1E1_wiso, diff="missing",
                                    lapse_rate=-0.002):
    
    # extract d18op
    exp_d18op_data = extract_var(Dataset=exp_data, varname="d18op", units="per mil", Dataset_wiso=exp_wiso)
    w1e1_d18op_data = extract_var(Dataset=W1E1_data, varname="d18op", units="per mil", Dataset_wiso=W1E1_wiso)
    
    #extract elevation 
    exp_elev_data = extract_var(Dataset=exp_data, varname="elev", units="m",)
    w1e1_elev_data = extract_var(Dataset=W1E1_data, varname="elev", units="m")
    
    # changes in elevation 
    elev_change = exp_elev_data - w1e1_elev_data
    
    expected_d18op_change = elev_change * lapse_rate # per mil 
    
    simulated_d18op_change = exp_d18op_data - w1e1_d18op_data
    
    missing = simulated_d18op_change - expected_d18op_change
    
    # calculate annual means of missing
    missing_alt = compute_lterm_mean(data=missing, time="annual")
    expected_d18op_change_alt = compute_lterm_mean(data=expected_d18op_change, time="annual")
    simulated_d18op_change_alt = compute_lterm_mean(data=simulated_d18op_change, time="annual")
    
    
    return_data_monthly = {"control_mon": w1e1_d18op_data, "simulated_change_mon": simulated_d18op_change,
                           "missing_mon": missing, "simulated_mon": exp_d18op_data, "expected_mon": expected_d18op_change}
    
    return_data_ltm = {"simulated_change_ltm": simulated_d18op_change_alt, "expected_change_ltm": expected_d18op_change_alt,
                       "missing_ltm": missing_alt}
    
    return_data = return_data_monthly | return_data_ltm
    
    
    return return_data



# extract data for the topo_exps
W2E0_278_diff = extract_d18Op_elev_and_analysis(exp_data=W2E0_278_data, exp_wiso=W2E0_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)  

W2E0_450_diff = extract_d18Op_elev_and_analysis(exp_data=W2E0_450_data, exp_wiso=W2E0_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso) 
      
W2E1_278_diff = extract_d18Op_elev_and_analysis(exp_data=W2E1_278_data, exp_wiso=W2E1_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)

W2E1_450_diff = extract_d18Op_elev_and_analysis(exp_data=W2E1_450_data, exp_wiso=W2E1_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso)

W2E15_278_diff = extract_d18Op_elev_and_analysis(exp_data=W2E15_278_data, exp_wiso=W2E15_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)

W2E15_450_diff = extract_d18Op_elev_and_analysis(exp_data=W2E15_450_data, exp_wiso=W2E15_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso)

W2E2_278_diff = extract_d18Op_elev_and_analysis(exp_data=W2E1_278_data, exp_wiso=W2E1_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)

W2E2_450_diff = extract_d18Op_elev_and_analysis(exp_data=W2E1_450_data, exp_wiso=W2E1_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso)

#plot the anomalies after correcting the lapse rate

def plot_d18Op_missing():
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
            
            plot_annual_mean(variable="$\delta^{18}$Op vs SMOW difference", data_alt=data[i].get("missing_ltm"), ax=axes[i],
                             cmap="PRGn", units="‰", vmax=5, vmin=-5, 
                            levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.30, 0.05, 0.45, 0.02], 
                            orientation="horizontal", plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", compare_data1=data[i].get("simulated_change_mon"), 
                            compare_data2=data[i].get("expected_mon"), max_pvalue=0.1, plot_stats=True,
                            hatches=".", title=label)
            
        else:
            plot_annual_mean(variable="$\delta^{18}$Op vs SMOW difference", data_alt=data[i].get("missing_ltm"), ax=axes[i],
                             cmap="PRGn", units="‰", vmax=5, vmin=-5, 
                            levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", compare_data1=data[i].get("simulated_change_mon"), 
                            compare_data2=data[i].get("expected_mon"), max_pvalue=0.1, plot_stats=True,
                            hatches=".", title=label)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    plt.savefig(os.path.join(path_to_plots, "d18Op_missing.svg"), format= "svg", bbox_inches="tight", dpi=600)

def plot_d18Op_simulated():
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
            
            plot_annual_mean(variable="$\delta^{18}$Op vs SMOW difference", data_alt=data[i].get("simulated_change_ltm"), 
                             ax=axes[i], cmap="PRGn", units="‰", vmax=8, vmin=-8, 
                            levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.30, 0.05, 0.45, 0.02], 
                            orientation="horizontal", plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", compare_data1=data[i].get("simulated_mon"), 
                            compare_data2=data[i].get("control_mon"), max_pvalue=0.1, plot_stats=True,
                            hatches=".", title=label)
            
        else:
            plot_annual_mean(variable="$\delta^{18}$Op vs SMOW difference", data_alt=data[i].get("simulated_change_ltm"), ax=axes[i],
                             cmap="PRGn", units="‰", vmax=8, vmin=-8, 
                            levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", compare_data1=data[i].get("simulated_mon"), 
                            compare_data2=data[i].get("control_mon"), max_pvalue=0.1, plot_stats=True,
                            hatches=".", title=label)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    plt.savefig(os.path.join(path_to_plots, "d18Op_simulated_change.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
def plot_d18Op_expected():
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
            
            plot_annual_mean(variable="$\delta^{18}$Op vs SMOW difference", data_alt=data[i].get("expected_change_ltm"), 
                             ax=axes[i], cmap="PRGn", units="‰", vmax=8, vmin=-8, 
                            levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.30, 0.05, 0.45, 0.02], 
                            orientation="horizontal", plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", title=label)
            
        else:
            plot_annual_mean(variable="$\delta^{18}$Op vs SMOW difference", data_alt=data[i].get("expected_change_ltm"), ax=axes[i],
                             cmap="PRGn", units="‰", vmax=8, vmin=-8, 
                            levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", title=label)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    plt.savefig(os.path.join(path_to_plots, "d18Op_expected_change.svg"), format= "svg", bbox_inches="tight", dpi=600)



if __name__ == "__main__":
    plot_d18Op_missing()
    plot_d18Op_expected()
    plot_d18Op_simulated()

