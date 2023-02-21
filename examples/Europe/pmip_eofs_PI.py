# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 16:07:18 2023

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

main_path_pi = "D:/Datasets/CMIP6/PMIP/postprocessed/PI"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/PI"
path_to_files = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/PI"
main_path_echam = "D:/Datasets/Model_output_pst/PI"


awi_path = os.path.join(main_path_pi, "AWI-ESM-1-1-LR")
cesm_path = os.path.join(main_path_pi, "CESM2")
cesm_waccm_path = os.path.join(main_path_pi, "CESM2-WACCM-FV2")
ec_earth_path = os.path.join(main_path_pi, "EC-Earth3-LR")
giss_path = os.path.join(main_path_pi, "GISS-E2-1-G")
hadGEM_path = os.path.join(main_path_pi, "HadGEM3-GC31-LL")
inm_cm_path = os.path.join(main_path_pi, "INM-CM4-8")
ipsl_path = os.path.join(main_path_pi, "IPSL-CM6A-LR")
miroc_path = os.path.join(main_path_pi, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_pi, "MPI-ESM1-2-LR")
norESM_path = os.path.join(main_path_pi, "NorESM1-F")




def perform_eof_and_store(apply_varimax, filesname, path_to_data=None, vmax=15, vmin=-15, data=None):
    
    if data is None:
        data = read_from_path(path=path_to_data, filename="psl_monthly.nc", 
                              varname="psl", decode=True) / 100 #Pa --> hPa
    
    
    
    data_pcs = extract_eofs_data(data=data, figname=filesname, 
                                 units="hPa", variable="Mean Sea Level Pressure", vmax=vmax, vmin=vmin,
                                 path_to_plots=path_to_plots, apply_varimax=apply_varimax, save_files=True,
                                 filename=filesname, path_to_files=path_to_files)


echam_data = read_from_path(main_path_echam, "PI_1003_1017_monthly.nc", decode=True)    
pi_msl = extract_var(Dataset=echam_data, varname="slp", units="hPa")
    
labels_pi = ["AWI-ESM-1-1-LR", "CESM2", "CESM2-WACCM-FV2", "EC-Earth3-LR", "GISS-E2-1-G", 
             "HadGEM3-GC31-LL", "MIROC-ES2L", "MPI-ESM1-2-LR", "NorESM1-F"]

data_paths = [awi_path, cesm_path, cesm_waccm_path, ec_earth_path, giss_path, hadGEM_path, 
              miroc_path, mpi_esm_path, norESM_path] 


# for i,model in enumerate(labels_pi):
#     perform_eof_and_store(apply_varimax=False, filesname= model + "_standard_eof_pi", path_to_data=data_paths[i],
#                           vmax=20, vmin=-20)
#     perform_eof_and_store(apply_varimax=True, filesname= model + "_varimax_eof_pi", path_to_data=data_paths[i])
    


    
perform_eof_and_store(apply_varimax=False, filesname="ECHAM5-wiso" + "_standard_eof_pi", data=pi_msl,
                      vmax=20, vmin=-20)
perform_eof_and_store(apply_varimax=True, filesname="ECHAM5-wiso" + "_varimax_eof_pi", data=pi_msl,
                      vmax=20, vmin=-20)