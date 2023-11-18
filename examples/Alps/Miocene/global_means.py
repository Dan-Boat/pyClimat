# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:50:30 2023

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
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff, compute_spatial_means
from pyClimat.variables import extract_var

path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"
W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

# reading data 
# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_PI_data, W1E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)


# temp means (trials)
def print_global_means(dataset, wiso):
    glob_temp = compute_spatial_means(dataset=dataset, varname="temp2", units="°C")
    glob_temp_land = compute_spatial_means(dataset=dataset, varname="temp2", units="°C", land=True)
    glob_temp_ocean = compute_spatial_means(dataset=dataset, varname="temp2", units="°C", ocean=True)
    
    glob_prec = compute_spatial_means(dataset=dataset, varname="prec", units="mm/month")
    glob_prec_land = compute_spatial_means(dataset=dataset, varname="prec", units="mm/month", land=True)
    glob_prec_ocean = compute_spatial_means(dataset=dataset, varname="prec", units="mm/month", ocean=True)
    
    glob_d18Op = compute_spatial_means(dataset=dataset, varname="d18op", units="per mil", 
                                       Dataset_wiso= wiso)
    glob_d18Op_land = compute_spatial_means(dataset=dataset, varname="d18op", units="per mil",
                                            Dataset_wiso= wiso, land=True)
    glob_d18Op_ocean = compute_spatial_means(dataset=dataset, varname="d18op", units="per mil", 
                                             Dataset_wiso= wiso, ocean=True)
    
    print(" Global temperature: {:.2f} , for land: {:.2f} , for ocean: {:.2f} (°C)".format(glob_temp, glob_temp_land,
                                                                                      glob_temp_ocean))
    
    print(" Global precipitation: {:.2f} , for land: {:.2f} , for ocean: {:.2f} (mm/month)".format(glob_prec, glob_prec_land,
                                                                                      glob_prec_ocean))
    
    print(" Global d18Op: {:.2f} , for land: {:.2f} , for ocean: {:.2f} (per mil)".format(glob_d18Op, glob_d18Op_land,
                                                                                      glob_d18Op_ocean))
    
    



print_global_means(dataset=W1E1_PI_data, wiso=W1E1_PI_wiso)
print_global_means(dataset=W1E1_278_data, wiso=W1E1_278_wiso)
print_global_means(dataset=W1E1_450_data, wiso=W1E1_450_wiso)
