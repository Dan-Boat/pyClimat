# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 10:54:03 2023

@author: dboateng
Extract the relevant teleconnections for the midHolocene using the PMIP4
model simulations
Store the eofs (as corr or corvar), pcs, variance and the plot for the first four modes
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

main_path_mh = "D:/Datasets/CMIP6/PMIP/postprocessed/MH"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/MH"
path_to_files = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/MH"
main_path_echam = "D:/Datasets/Model_output_pst/MH"


awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
cesm_path = os.path.join(main_path_mh, "CESM2")
ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")



def perform_eof_and_store(apply_varimax, filesname, path_to_data=None, vmax=15, vmin=-15, data=None):
    
    if data is None:
        data = read_from_path(path=path_to_data, filename="psl_monthly.nc", 
                              varname="psl", decode=True) / 100 #Pa --> hPa
    
    
    
    data_pcs = extract_eofs_data(data=data, figname=filesname, 
                                 units="hPa", variable="Mean Sea Level Pressure", vmax=vmax, vmin=vmin,
                                 path_to_plots=path_to_plots, apply_varimax=apply_varimax, save_files=True,
                                 filename=filesname, path_to_files=path_to_files)


echam_data = read_from_path(main_path_echam, "MH_1003_1017_monthly.nc", decode=True)    
mh_msl = extract_var(Dataset=echam_data, varname="slp", units="hPa")
    
labels_mh = ["AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
              "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"]

data_paths = [awi_path, cesm_path, ec_earth_path, giss_path, ipsl_path, miroc_path, mpi_esm_path] 


# for i,model in enumerate(labels_mh):
#     perform_eof_and_store(apply_varimax=False, filesname= model + "_standard_eof_mh", path_to_data=data_paths[i],
#                           vmax=20, vmin=-20)
#     perform_eof_and_store(apply_varimax=True, filesname= model + "_varimax_eof_mh", path_to_data=data_paths[i])
    


    
perform_eof_and_store(apply_varimax=False, filesname="ECHAM5-wiso" + "_standard_eof_mh", data=mh_msl,
                      vmax=20, vmin=-20)
perform_eof_and_store(apply_varimax=True, filesname="ECHAM5-wiso" + "_varimax_eof_mh", data=mh_msl,
                      vmax=15, vmin=-15)