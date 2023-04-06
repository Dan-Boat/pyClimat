# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 11:21:30 2023

@author: dboateng
"""

import os 
import numpy as np 
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt 
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.analysis import extract_var
from pyClimat.stats import EOF_standard
from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance

from pyClimat.data import read_ECHAM_processed, read_from_path, read_ERA_processed
from main import extract_eofs_data

main_path_plio = "D:/Datasets/CMIP6/PMIP/postprocessed/PLIO"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/PLIO"
path_to_files = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/PLIO"
main_path_echam = "D:/Datasets/Model_output_pst/PLIO"

cesm_path = os.path.join(main_path_plio, "CESM2")
ec_earth_path = os.path.join(main_path_plio, "EC-Earth3-LR")
giss_path = os.path.join(main_path_plio, "GISS-E2-1-G")
ipsl_path = os.path.join(main_path_plio, "IPSL-CM6A-LR")
norESM_path = os.path.join(main_path_plio, "NorESM1-F")

def perform_eof_and_store(apply_varimax=False, filesname=None, path_to_data=None, vmax=20, vmin=-20, data=None,
                          method="xeofs", season="DJF"):
    
    if data is None:
        data = read_from_path(path=path_to_data, filename="psl_monthly.nc", 
                              varname="psl", decode=True) / 100 #Pa --> hPa
    
    
    
    data_pcs = extract_eofs_data(data=data, figname=filesname, 
                                 units="hPa", variable="Mean Sea Level Pressure", vmax=vmax, vmin=vmin,
                                 path_to_plots=path_to_plots, apply_varimax=apply_varimax, save_files=True,
                                 filename=filesname, path_to_files=path_to_files, standardize=False, 
                                 monthly_anomalies=True, method=method, season=season, is_era=False)

season = "DJF"

echam_data = read_from_path(main_path_echam, "PLIO_1003_1017_monthly.nc", decode=True) 
PLIO_msl = extract_var(Dataset=echam_data, varname="slp", units="hPa")

labels_plio = ["CESM2", "EC-Earth3-LR", "GISS-E2-1-G", "IPSL-CM6A-LR", "NorESM1-F"]

data_paths = [cesm_path, ec_earth_path, giss_path, ipsl_path, norESM_path]

for i,model in enumerate(labels_plio):
    perform_eof_and_store(apply_varimax=False, filesname= model + "_standard_eof_plio" + season, path_to_data=data_paths[i],
                          vmax=20, vmin=-20)
    perform_eof_and_store(apply_varimax=True, filesname= model + "_varimax_eof_plio" + season, path_to_data=data_paths[i])
    


    
# perform_eof_and_store(apply_varimax=False, filesname="ECHAM5-wiso" + "_standard_eof_plio", data=PLIO_msl,
#                       vmax=20, vmin=-20)
# perform_eof_and_store(apply_varimax=True, filesname="ECHAM5-wiso" + "_varimax_eof_plio", data=PLIO_msl,
#                       vmax=15, vmin=-15)