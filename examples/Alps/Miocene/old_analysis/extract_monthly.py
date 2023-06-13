# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:58:41 2023

@author: dboateng

This script extracts the montly means of surface temperature, deep soil temperature and
soil wetness across the Alps 
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

import calendar


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, compute_lterm_diff, extract_transect


#define paths

path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/"


# read data

W1E1_Mio278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_Mio450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

# experiments
W2E1_Mio278_filename = "a017_hpc-bw_e5w2.3_t159_MIO_W2E1_278ppm_t159l31.6h"
W2E1_Mio450_filename = "a016_hpc-bw_e5w2.3_t159_MIO_W2E1_450ppm_t159l31.6h"

W2E0_Mio278_filename="a019_hpc-bw_e5w2.3_t159_MIO_W2E0_278ppm_t159l31.6h"
W2E0_Mio450_filename="a018_hpc-bw_e5w2.3_t159_MIO_W2E0_450ppm_t159l31.6h"


W2E2_Mio278_filename="a020_dkrz-levante_e5w2.3_t159_MIO_W2E2_278ppm_t159l31.6h"
W2E2_Mio450_filename="a021_dkrz-levante_e5w2.3_t159_MIO_W2E2_450ppm_t159l31.6h"

years = "1003_1017"
period = "1m"

#read data
W1E1_278_data = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_Mio278_filename, years=years,
                                     period=period, read_wiso=False)

W1E1_450_data = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_Mio450_filename, years=years,
                                     period=period, read_wiso=False)

W2E2_278_data = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio278_filename, years=years,
                                     period=period, read_wiso=False)

W2E2_450_data = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio450_filename, years=years,
                                     period=period, read_wiso=False)

# function to calculate the mean, std, max, min
def weighted_mean_std_max_min(data):
    weights = np.cos(np.deg2rad(data.lat))
    weights.name = "weights"
    
    data_weighted = data.weighted(weights)
    
    mean = data_weighted.mean(dim=("lat", "lon"), skipna=True)
    std = data_weighted.std(dim=("lat", "lon"), skipna=True)
    max_value = data.max(dim=("lat", "lon"), skipna=True)
    min_value = data.min(dim=("lat", "lon"), skipna=True)
    
    stats = [mean, std, max_value, min_value]
    return stats

#function that extracts the variables, then their transects and monthly means
def extract_monthly(dataset, varname, units=None, path_to_save=None, filename=None):
    var_data = extract_var(Dataset=dataset, varname=varname, units=units)
    
    max_lon, max_lat, min_lon, min_lat = 7, 45, 4, 40
    
    region_data = extract_transect(data=var_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat,
                                   minlat=min_lat, sea_land_mask=True, Dataset=dataset)
    
    stat_data = weighted_mean_std_max_min(region_data)
    
    month_names = [calendar.month_abbr[im+1] for im in np.arange(12)]
    
    column_names = ["mean", "std", "max", "min"]
    df = pd.DataFrame(index=month_names, columns= column_names)
    
    for i,name in enumerate(column_names):
        df[name] = stat_data[i]
    
    
    # save data in csv
    if filename is not None:
        filename = filename + ".csv"
        
        df.to_csv(filename)
        
    
    
# mio_450_datasets = [W1E1_450_data, W2E2_450_data]
# exp_names = ["W1E1_Mio450", "W2E2_Mio450"]

mio_278_datasets = [W1E1_278_data, W2E2_278_data]
exp_names = ["W1E1_Mio278", "W2E2_Mio278"]
varnames = ["tsurf", "ws", "tsoil"]   

for i,exp_name in enumerate(exp_names):
    
    for varname in varnames:
        units="°C"
        
        if varname =="ws":
            units=None
            
        extract_monthly(dataset=mio_278_datasets[i], varname=varname, units=units, path_to_save=path_to_plots,
                        filename=exp_name + "_" + varname)    
    
    
# plots the monthly means and variablity (using max, min or thier std)


#save the files in csv 