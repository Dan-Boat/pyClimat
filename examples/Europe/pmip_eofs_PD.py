# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 16:07:36 2023

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

def perform_eof_and_store(apply_varimax, filesname, path_to_data=None, vmax=15, vmin=-15, data=None, 
                          method="xeofs", season="JJA", is_era=False):
    
    if data is None:
        data = read_from_path(path=path_to_data, filename="psl_monthly.nc", 
                              varname="psl", decode=True) / 100 #Pa --> hPa
    
    
    
    data_pcs = extract_eofs_data(data=data, figname=filesname, 
                                 units="hPa", variable="Mean Sea Level Pressure", vmax=vmax, vmin=vmin,
                                 path_to_plots=path_to_plots, apply_varimax=apply_varimax, save_files=True,
                                 filename=filesname, path_to_files=path_to_files, standardize=False,
                                 monthly_anomalies=True, method=method, season=season, is_era=is_era)
    
    

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"
main_path_echam = "D:/Datasets/Model_output_pst/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021/"
path_to_files = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/PD"


echam_data = read_from_path(main_path_echam, "PD_1980_2014_monthly.nc", decode=True)    
pd_msl = extract_var(Dataset=echam_data, varname="slp", units="hPa")


ERA5_msl_path = os.path.join(ERA5_path, "msl_monthly.nc")
ERA5_msl = read_ERA_processed(path=ERA5_msl_path, varname="msl") / 100 #Pa --> hPa

season="JJA"

perform_eof_and_store(apply_varimax=False, filesname="ECHAM5-wiso" + "_standard_eof_pd_" + season, data=pd_msl,
                      vmax=10, vmin=-10, method="xeofs")

perform_eof_and_store(apply_varimax=True, filesname="ECHAM5-wiso" + "_varimax_eof_pd_" + season, data=pd_msl,
                      vmax=10, vmin=-10, method="xeofs")


perform_eof_and_store(apply_varimax=False, filesname="ERA5" + "_standard_eof_" + season, data=ERA5_msl,
                      vmax=10, vmin=-10, method="xeofs")
perform_eof_and_store(apply_varimax=True, filesname="ERA5" + "_varimax_eof_" + season, data=ERA5_msl,
                      vmax=10, vmin=-10, method="xeofs")