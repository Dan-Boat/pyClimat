# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 18:07:14 2023

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
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/atmospheric_dynamics/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
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





def extract_prec_and_analysis(exp_name, control_exp_name, return_climatologies=False):
    
   
    # read data (long-term means) - change to monthly for stats 
    years = "1003_1017"
    period = "1m"
    
    exp_data = read_ECHAM_processed(main_path=path_to_data, exp_name=exp_name, years=years, 
                                          period=period, read_wiso=False)
    
    control_data = read_ECHAM_processed(main_path=path_to_data, exp_name=control_exp_name, years=years, 
                                          period=period, read_wiso=False)
    

    
    exp_prec_data = extract_var(Dataset=exp_data, varname="prec", units="mm/month")
    w1e1_prec_data = extract_var(Dataset=control_data, varname="prec", units="mm/month")
    
    # extract aprc
    exp_aprc_data = extract_var(Dataset=exp_data, varname="aprc", units="mm/month")
    w1e1_aprc_data = extract_var(Dataset=control_data, varname="aprc", units="mm/month")
    
    # extract aprl
    exp_aprl_data = extract_var(Dataset=exp_data, varname="aprl", units="mm/month")
    w1e1_aprl_data = extract_var(Dataset=control_data, varname="aprl", units="mm/month")
    
    if return_climatologies:
        simulated_prec_alt_main = compute_lterm_mean(data=exp_prec_data, time="season", 
                                                     season_calendar="standard")
        simulated_prec_alt_ctl = compute_lterm_mean(data=w1e1_prec_data, time="season", 
                                                     season_calendar="standard")
        simulated_aprc_alt_main = compute_lterm_mean(data=exp_aprc_data, time="season", 
                                                     season_calendar="standard")
        simulated_aprc_alt_ctl = compute_lterm_mean(data=w1e1_aprc_data, time="season", 
                                                     season_calendar="standard")
        simulated_aprl_alt_main = compute_lterm_mean(data=exp_aprl_data, time="season", 
                                                     season_calendar="standard")
        simulated_aprl_alt_ctl = compute_lterm_mean(data=w1e1_aprl_data, time="season", 
                                                     season_calendar="standard")
        
        return_data_prec = {"control": simulated_prec_alt_main, "main": simulated_prec_alt_ctl}
        return_data_aprc = {"control": simulated_aprc_alt_main, "main": simulated_aprc_alt_ctl}
        return_data_aprl = {"control": simulated_aprl_alt_main, "main": simulated_aprl_alt_ctl}
        
        return return_data_prec, return_data_aprc, return_data_aprl
    
    else:
        
        simulated_prec_change_alt = compute_lterm_diff(data_control=w1e1_prec_data, data_main=exp_prec_data,
                                              time="season", season_calendar="standard")
        
        simulated_aprc_change_alt = compute_lterm_diff(data_control=w1e1_aprc_data, data_main=exp_aprc_data,
                                              time="season", season_calendar="standard")
        
        simulated_aprl_change_alt = compute_lterm_diff(data_control=w1e1_aprl_data, data_main=exp_aprl_data,
                                              time="season", season_calendar="standard")
        
        
        return_data_prec = {"control_mon": w1e1_prec_data, "simulated_mon": exp_prec_data,
                            "simulated_change_ltm": simulated_prec_change_alt}
        
        return_data_aprc = {"control_mon": w1e1_aprc_data, "simulated_mon": exp_aprc_data,
                            "simulated_change_ltm": simulated_aprc_change_alt}
        
        return_data_aprl = {"control_mon": w1e1_aprl_data, "simulated_mon": exp_aprl_data,
                            "simulated_change_ltm": simulated_aprl_change_alt}
        
        
        
        return return_data_prec, return_data_aprc, return_data_aprl
    
    
def plot_prec_simulated(data, labels, season, filename, vmax=None, vmin=None, levels=None, level_tick=None):
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    fig,((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(22,28), 
                                                                    subplot_kw={"projection":projection})
    
    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    
    
    
    for i,label in enumerate(labels):
        
        if i == 0:
            
            plot_annual_mean(variable="Precipitation", data_alt=data[i].get("control").sel(season=season), ax=axes[i],
                             cmap=YlGnBu, units="mm/month", vmax=200, vmin=20, 
                            levels=22, level_ticks=10, add_colorbar=True, cbar_pos= [0.09, 0.05, 0.35, 0.02], 
                            orientation="horizontal", plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", title=label)
            
        elif i == 1:
            
            plot_annual_mean(variable="Precipitation difference", data_alt=data[i].get("simulated_change_ltm").sel(season=season), ax=axes[i],
                             cmap=BrBG, units="mm/month", vmax=80, vmin=-80, 
                            levels=22, level_ticks=9, add_colorbar=True, cbar_pos= [0.45, 0.05, 0.35, 0.02], 
                            orientation="horizontal", plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", title=label)
            
        else:
            plot_annual_mean(variable="Precipitation difference", data_alt=data[i].get("simulated_change_ltm").sel(season=season), ax=axes[i],
                             cmap=BrBG, units="mm/month", vmax=80, vmin=-80, 
                            levels=22, level_ticks=9, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                            left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                            plot_projection=projection, domain="Europe", title=label)
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    plt.savefig(os.path.join(path_to_plots, filename), format= "svg", bbox_inches="tight", dpi=600)
    


Mio278_PI_prec, Mio278_PI_aprc, Mio278_PI_aprl = extract_prec_and_analysis(exp_name=W1E1_278_filename, 
                                                                           control_exp_name=W1E1_PI_filename,
                                                                           return_climatologies=True)

Mio278_PI_diff_prec, Mio278_PI_diff_aprc, Mio278_PI_diff_aprl = extract_prec_and_analysis(exp_name=W1E1_278_filename, 
                                                                           control_exp_name=W1E1_PI_filename)

Mio278_W2E0_diff_prec, Mio278_W2E0_diff_aprc, Mio278_W2E0_diff_aprl = extract_prec_and_analysis(exp_name=W2E0_Mio278_filename, 
                                                                           control_exp_name=W1E1_278_filename)
Mio278_W2E1_diff_prec, Mio278_W2E1_diff_aprc, Mio278_W2E1_diff_aprl = extract_prec_and_analysis(exp_name=W2E1_Mio278_filename, 
                                                                           control_exp_name=W1E1_278_filename)
Mio278_W2E15_diff_prec, Mio278_W2E15_diff_aprc, Mio278_W2E15_diff_aprl = extract_prec_and_analysis(exp_name=W2E15_Mio278_filename, 
                                                                           control_exp_name=W1E1_278_filename)
Mio278_W2E2_diff_prec, Mio278_W2E2_diff_aprc, Mio278_W2E2_diff_aprl = extract_prec_and_analysis(exp_name=W2E2_Mio278_filename, 
                                                                           control_exp_name=W1E1_278_filename)



labels = ["(a) W1E1 (MIO 278ppm)", "(b) W1E1 (MIO 278ppm) - W1E1 (PI)", "(c) W2E0 - W1E1 (MIO 278ppm)",  
          "(d) W2E1 - W1E1 (MIO 278ppm)", "(e) W2E1.5 - W1E1 (MIO 278ppm)", "(f) W2E2 - W1E1 (MIO 278ppm)"]

data_prec = [Mio278_PI_prec, Mio278_PI_diff_prec, Mio278_W2E0_diff_prec, Mio278_W2E1_diff_prec,
        Mio278_W2E15_diff_prec, Mio278_W2E2_diff_prec]

data_aprc = [Mio278_PI_aprc, Mio278_PI_diff_aprc, Mio278_W2E0_diff_aprc, Mio278_W2E1_diff_aprc,
        Mio278_W2E15_diff_aprc, Mio278_W2E2_diff_aprc]

data_aprl = [Mio278_PI_aprl, Mio278_PI_diff_aprl, Mio278_W2E0_diff_aprl, Mio278_W2E1_diff_aprl,
        Mio278_W2E15_diff_aprl, Mio278_W2E2_diff_aprl]

plot_prec_simulated(data_prec, labels, "JJA", "prec_diff_JJA.svg")
plot_prec_simulated(data_prec, labels, "DJF", "prec_diff_DJF.svg")

plot_prec_simulated(data_aprc, labels, "JJA", "aprc_diff_JJA.svg")
plot_prec_simulated(data_aprc, labels, "DJF", "aprc_diff_DJF.svg")

plot_prec_simulated(data_aprl, labels, "JJA", "aprl_diff_JJA.svg")
plot_prec_simulated(data_aprl, labels, "DJF", "aprl_diff_DJF.svg")



