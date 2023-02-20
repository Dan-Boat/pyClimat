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

main_path_mh = "D:/Datasets/CMIP6/PMIP/postprocessed/PLIO"

cesm_path = os.path.join(main_path_plio, "CESM2")
ec_earth_path = os.path.join(main_path_plio, "EC-Earth3-LR")
giss_path = os.path.join(main_path_plio, "GISS-E2-1-G")
ipsl_path = os.path.join(main_path_plio, "IPSL-CM6A-LR")
norESM_path = os.path.join(main_path_plio, "NorESM1-F")


data = read_from_path(path=awi_path, filename="psl_monthly.nc", 
                      varname="psl", decode=True) / 100 #Pa --> hPa



data_pcs = extract_eofs_data(data=data, figname="awi_mh_msl_varimax_wNAO", 
                             units="hPa", variable="Mean Sea Level Pressure", vmax=10, vmin=-10)