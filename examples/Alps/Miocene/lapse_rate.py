# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 15:30:34 2023

@author: dboateng
1. Calculate the lapse rate across the Alps for the control experiments
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
from pyClimat.plots import scatter_plot_laspe_rate
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, extract_transect, linregression


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"


# read data

CTL_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

years = "1003_1017"
period = "1m"

CTL_data, CTL_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=CTL_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

def extract_relevant_vars_sections(data, wiso):

    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    elev = extract_var(Dataset=data , varname="elev", units="m")

    # compute annual means
    d18op_alt = compute_lterm_mean(data=d18op, time="annual")
    elev_alt = compute_lterm_mean(data=elev, time="annual")
    
    maxlat_west, minlat_west, maxlon_west, minlon_west = 47, 44, 8, 1
    maxlat_south, minlat_south, maxlon_south, minlon_south = 47, 43, 15, 7.5
    maxlat_north, minlat_north, maxlon_north, minlon_north = 50, 46.5, 16, 5
    
    elev_north = extract_transect(data=elev_alt, maxlon=maxlon_north, minlon=minlon_north , 
                                  maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, 
                                  Dataset=data)
    
    elev_south = extract_transect(data=elev_alt, maxlon=maxlon_south, minlon=minlon_south, 
                                  maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, 
                                  Dataset=data)
    
    
    elev_west = extract_transect(data=elev_alt, maxlon=maxlon_west, minlon=minlon_west, 
                                  maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, 
                                  Dataset=data)
    
    
    d18op_north = extract_transect(data=d18op_alt, maxlon=maxlon_north, minlon=minlon_north , 
                                  maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, 
                                  Dataset=data)
    
    d18op_south = extract_transect(data=d18op_alt, maxlon=maxlon_south, minlon=minlon_south, 
                                  maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, 
                                  Dataset=data)
    
    
    d18op_west = extract_transect(data=d18op_alt, maxlon=maxlon_west, minlon=minlon_west, 
                                  maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, 
                                  Dataset=data)
   
   
    west_reg, west_df = linregression(data_x=elev_west, data_y=d18op_west, return_yhat=True)
    north_reg, north_df = linregression(data_x=elev_north, data_y=d18op_north, return_yhat=True)
    south_reg, south_df = linregression(data_x=elev_south, data_y=d18op_south, return_yhat=True)
    
    
    
    return north_reg, north_df, south_reg, south_df, west_reg, west_df





#plot 

apply_style(fontsize=22, style=None, linewidth=2)

def plot_lape_rate_per_section():

    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 8))
    
    #ax1 (west)
    
    scatter_plot_laspe_rate(ax=ax1, reg_params= CTL_reg_west , df_x_y_yhat=CTL_df_west , color=black, marker= "*", label= "CTL (PI)",
                           title="[A] West", xmax=1500, xmin=0,
                            ymax=-2, ymin= -12, bottom_labels=True)
    scatter_plot_laspe_rate(ax=ax1, reg_params= W1E1_278_reg_west , df_x_y_yhat=W1E1_278_df_west , color=red, marker= "D", label= "W1E1 (MIO 278ppm)",
                           bottom_labels=True)
    scatter_plot_laspe_rate(ax=ax1, reg_params= W1E1_450_reg_west , df_x_y_yhat=W1E1_450_df_west , color=green, marker= "^", label= "W1E1 (MIO 450ppm)",
                           bottom_labels=True)
    
    
    
    ax1.legend(frameon=True, fontsize=20, loc="upper left", framealpha=0.5, ncol=1)
    ax1.grid(visible=False)
    
    
    #ax2 (north)
    scatter_plot_laspe_rate(ax=ax2, reg_params= CTL_reg_north , df_x_y_yhat=CTL_df_north , color=black, marker= "*", label= "CTL (PI)",
                            left_labels=False, xmax=1500, xmin=0, title= "[B] North",
                             ymax=-2, ymin= -12,)
    scatter_plot_laspe_rate(ax=ax2, reg_params= W1E1_278_reg_north , df_x_y_yhat=W1E1_278_df_north , color=red, marker= "D", label= "W1E1 (MIO 278ppm)",
                            left_labels=False)
    
    scatter_plot_laspe_rate(ax=ax2, reg_params= W1E1_450_reg_north , df_x_y_yhat=W1E1_450_df_north , color=green, marker= "^", label= "W1E1 (MIO 450ppm)",
                            left_labels=False)
    
    
    ax2.legend(frameon=True, fontsize=20, loc="upper left", framealpha=0.5, ncol=1)
    ax2.grid(visible=False)
    
    #ax3 (south)
    scatter_plot_laspe_rate(ax=ax3, reg_params= CTL_reg_south , df_x_y_yhat=CTL_df_south , color=black, marker= "*", label= "CTL (PI)",
                            left_labels=False, xmax=1500, xmin=0, title= "[C] South",
                             ymax=-2, ymin= -12, add_label=True)
    scatter_plot_laspe_rate(ax=ax3, reg_params= W1E1_278_reg_south , df_x_y_yhat=W1E1_278_df_south , color=red, marker= "^", label= "W1E1 (MIO 278ppm)",
                           left_labels=False, add_label=True)
    scatter_plot_laspe_rate(ax=ax3, reg_params= W1E1_450_reg_south, df_x_y_yhat=W1E1_450_df_south, color=green, marker= "D", label= "W1E1 (MIO 450ppm)",
                           left_labels=False, add_label=True)
    
    
    ax3.legend(frameon=True, fontsize=20, loc="upper left", framealpha=0.5, ncol=1)
    ax3.grid(visible=False)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.88, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_plots, "lapse_rate_CTL_pi_mio.svg"), format= "svg", bbox_inches="tight", dpi=600)
 
    
if __name__ == '__main__':

    CTL_reg_north, CTL_df_north, CTL_reg_south, CTL_df_south, CTL_reg_west, CTL_df_west = extract_relevant_vars_sections(data=CTL_data, wiso=CTL_wiso)
    W1E1_278_reg_north, W1E1_278_df_north, W1E1_278_reg_south, W1E1_278_df_south, W1E1_278_reg_west, W1E1_278_df_west = extract_relevant_vars_sections(data=W1E1_278_data, wiso=W1E1_278_wiso)
    W1E1_450_reg_north, W1E1_450_df_north, W1E1_450_reg_south, W1E1_450_df_south, W1E1_450_reg_west, W1E1_450_df_west = extract_relevant_vars_sections(data=W1E1_450_data, wiso=W1E1_450_wiso)

    plot_lape_rate_per_section()                                                                                                                          