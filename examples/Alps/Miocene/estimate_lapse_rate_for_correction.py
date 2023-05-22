# -*- coding: utf-8 -*-
"""
Created on Sat May 20 15:06:02 2023

@author: dboateng
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
import seaborn as sns


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import scatter_plot_laspe_rate
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, extract_transect, linregression


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


def extract_vars_and_analysis(data, wiso, maxlon=20, maxlat=50, minlon=-1, minlat=41):

    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    temp = extract_var(Dataset=data , varname="temp2", units="°C")
    elev = extract_var(Dataset=data , varname="elev", units="m")

    # compute annual means
    d18op_alt = compute_lterm_mean(data=d18op, time="annual")
    temp_alt = compute_lterm_mean(data=temp, time="annual")
    elev_alt = compute_lterm_mean(data=elev, time="annual")
    
    elev_transect = extract_transect(data=elev_alt, maxlon=maxlon, minlon=minlon, 
                                  maxlat=maxlat, minlat=minlat, sea_land_mask=True, 
                                  Dataset=data)
    
    
    
    d18op_transect = extract_transect(data=d18op_alt, maxlon=maxlon, minlon=minlon, 
                                  maxlat=maxlat, minlat=minlat, sea_land_mask=True, 
                                  Dataset=data)
    
    temp_transect = extract_transect(data=temp_alt, maxlon=maxlon, minlon=minlon, 
                                  maxlat=maxlat, minlat=minlat, sea_land_mask=True, 
                                  Dataset=data)
   
   
    reg_d18op = linregression(data_x=elev_transect, data_y=d18op_transect, return_yhat=False)
    
    reg_temp = linregression(data_x=elev_transect, data_y=temp_transect, return_yhat=False)
    
    
    return_data = {"reg":reg_temp, "reg_d18op": reg_d18op}
    
    return return_data


P1_reg = extract_vars_and_analysis(data=W1E1_PI_data, wiso=W1E1_PI_wiso)
Mio278_reg = extract_vars_and_analysis(data=W1E1_278_data, wiso=W1E1_278_wiso)
Mio450_reg = extract_vars_and_analysis(data=W1E1_450_data, wiso=W1E1_450_wiso)


print(None)



