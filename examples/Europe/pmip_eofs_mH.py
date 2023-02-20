# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 10:54:03 2023

@author: dboateng
Extract the relevant teleconnections for the midHolocene using the PMIP4
model simulations
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
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"


awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
cesm_path = os.path.join(main_path_mh, "CESM2")
ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")




data = read_from_path(path=awi_path, filename="psl_monthly.nc", 
                      varname="psl", decode=True) / 100 #Pa --> hPa



data_pcs = extract_eofs_data(data=data, figname="awi_mh_msl_varimax_wNAO", 
                             units="hPa", variable="Mean Sea Level Pressure", vmax=15, vmin=-15,
                             path_to_plots=path_to_plots, apply_varimax=False)