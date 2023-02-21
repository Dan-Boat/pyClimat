# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 11:21:17 2023

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

main_path_lgm = "D:/Datasets/CMIP6/PMIP/postprocessed/LGM"

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/LGM"
path_to_files = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/LGM"
main_path_echam = "D:/Datasets/Model_output_pst/LGM"


awi_path = os.path.join(main_path_lgm, "AWI-ESM-1-1-LR")
cesm_waccm_path = os.path.join(main_path_lgm, "CESM2-WACCM-FV2")
inm_cm_path = os.path.join(main_path_lgm, "INM-CM4-8")
miroc_path = os.path.join(main_path_lgm, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_lgm, "MPI-ESM1-2-LR")

def perform_eof_and_store(apply_varimax, filesname, path_to_data=None, vmax=15, vmin=-15, data=None):
    
    if data is None:
        data = read_from_path(path=path_to_data, filename="psl_monthly.nc", 
                              varname="psl", decode=True) / 100 #Pa --> hPa
    
    
    
    data_pcs = extract_eofs_data(data=data, figname=filesname, 
                                 units="hPa", variable="Mean Sea Level Pressure", vmax=vmax, vmin=vmin,
                                 path_to_plots=path_to_plots, apply_varimax=apply_varimax, save_files=True,
                                 filename=filesname, path_to_files=path_to_files)


echam_data = read_from_path(main_path_echam, "LGM_1003_1017_monthly.nc", decode=True) 
lgm_msl = extract_var(Dataset=echam_data, varname="slp", units="hPa") 

labels_lgm = ["MPI-ESM1-2-LR"]#["AWI-ESM-1-1-LR", "CESM2-WACCM-FV2", "INM-CM4-8", "MIROC-ES2L", "MPI-ESM1-2-LR"] 


data_paths = [mpi_esm_path]#[awi_path, cesm_waccm_path, inm_cm_path, miroc_path, mpi_esm_path]

# for i,model in enumerate(labels_lgm):
#     perform_eof_and_store(apply_varimax=False, filesname= model + "_standard_eof_lgm", path_to_data=data_paths[i],
#                           vmax=20, vmin=-20)
#     perform_eof_and_store(apply_varimax=True, filesname= model + "_varimax_eof_lgm", path_to_data=data_paths[i])
    


    
perform_eof_and_store(apply_varimax=False, filesname="ECHAM5-wiso" + "_standard_eof_lgm", data=lgm_msl,
                      vmax=20, vmin=-20)
perform_eof_and_store(apply_varimax=True, filesname="ECHAM5-wiso" + "_varimax_eof_lgm", data=lgm_msl,
                      vmax=15, vmin=-15)