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

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"

awi_path = os.path.join(main_path_lgm, "AWI-ESM-1-1-LR")
cesm_waccm_path = os.path.join(main_path_lgm, "CESM2-WACCM-FV2")
inm_cm_path = os.path.join(main_path_lgm, "INM-CM4-8")
miroc_path = os.path.join(main_path_lgm, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_lgm, "MPI-ESM1-2-LR")

data = read_from_path(path=cesm_waccm_path, filename="psl_monthly.nc", 
                      varname="psl", decode=True) / 100 #Pa --> hPa



data_pcs = extract_eofs_data(data=data, figname="cesm_waccm_lgm_msl_varimax_wNAO", 
                             units="hPa", variable="Mean Sea Level Pressure", vmax=15, vmin=-15,
                             path_to_plots=path_to_plots, apply_varimax=False)