# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 16:39:21 2023

@author: dboateng
"""

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


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"


# define paths 
main_path = "D:/Datasets/Model_output_pst/"
lgm_path = os.path.join(main_path, "LGM", "Inputs")
plio_path = os.path.join(main_path, "PLIO", "Inputs")
mh_path = os.path.join(main_path, "MH", "Inputs")
pi_path = os.path.join(main_path, "PI", "Inputs")

filename_sst = "SST.nc"
# read long-term means

PI_data = read_from_path(pi_path, filename_sst)
LGM_data = read_from_path(lgm_path, filename_sst)
PLIO_data = read_from_path(plio_path, filename_sst)
MH_data = read_from_path(mh_path, filename_sst)

def extract_vars_and_analysis(data, pi_data):
    
    sst = extract_var(Dataset=data , varname="sst", units="°C")
    
    
    sst_pi = extract_var(Dataset=pi_data , varname="sst", units="°C")
    
    
    
    #compute climatologies difference
    
    sst_diff = compute_lterm_diff(data_control=sst_pi, data_main=sst, time="annual")
    
    
    
    
    return_data = {"SST":sst_diff}
    
    return return_data


lgm_diff = extract_vars_and_analysis(data=LGM_data, pi_data=PI_data)