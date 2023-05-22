# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:42:02 2023

@author: dboateng
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_var, extract_transect
from pyClimat.utils import extract_region
from main import extract_eofs_data

path_to_plots="C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium"
path_to_files="C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium"

#paths
main_path = "D:/Datasets/iGCM_datasets/"
giss_data_path = "D:/Datasets/iGCM_datasets/GISS_wiso_vars.nc"
hadesm_data_path = "D:/Datasets/iGCM_datasets/HADCM3_wiso_vars.nc"
mpi_esm_data_path = "D:/Datasets/iGCM_datasets/ECHAM5_wiso_vars.nc"
ccsm_data_path = "D:/Datasets/iGCM_datasets/CCSM_wiso_vars.nc"
cesm_data_path = "D:/Datasets/iGCM_datasets/CESM_wiso_vars.nc"


# read data 

mpi_esm_data = read_from_path(path=main_path, filename= "ECHAM5_wiso_vars.nc", decode=True,)
giss_data = read_from_path(path=main_path, filename="GISS_wiso_vars.nc", decode=True)
hadesm_data = read_from_path(path=main_path, filename="HADCM3_wiso_vars.nc", decode=True)
ccsm_data = read_from_path(path=main_path, filename="CCSM_wiso_vars.nc", decode=True)
cesm_data = read_from_path(path=main_path, filename="CESM_wiso_vars.nc", decode=True)


def perform_eof_and_store(apply_varimax=False, filename=None, path_to_data=None, vmax=20, vmin=-20, data=None,
                          method="xeofs", season="DJF"):
    
    filename_data = filename + "_wiso_vars.nc"
    data = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="slp") / 100 
    
    
    
    if hasattr(data, "lat_2"):
        data = data.rename({"lat_2":"lat"})
        
    if hasattr(data, "lon_2"):
        data = data.rename({"lon_2":"lon"})
        
    data_pcs = extract_eofs_data(data=data, figname=filename, 
                                 units="hPa", variable="Mean Sea Level Pressure", vmax=vmax, vmin=vmin,
                                 path_to_plots=path_to_plots, apply_varimax=apply_varimax, save_files=True,
                                 filename=filename, path_to_files=path_to_files, standardize=False, 
                                 monthly_anomalies=True, method=method, season=season, is_era=False)



filenames = ["ECHAM5"]

for filename in filenames:
    perform_eof_and_store(filename=filename, path_to_data=main_path, season="DJF")
    #perform_eof_and_store(filename=filename, path_to_data=main_path, season="JJA")


print(None)


