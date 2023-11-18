# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:50:30 2023

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
from pyClimat.data import read_from_path
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff, compute_spatial_means
from pyClimat.variables import extract_var

# define paths 
main_path = "D:/Datasets/Model_output_pst/"
lgm_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
plio_path = os.path.join(main_path, "PLIO", "MONTHLY_MEANS")
mh_path = os.path.join(main_path, "MH", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")


filename_lterm = "1003_1017_1m_mlterm.nc"
# read long-term means

PI_data = read_from_path(pi_path, filename_lterm)
LGM_data = read_from_path(lgm_path, filename_lterm)
PLIO_data = read_from_path(plio_path, filename_lterm)
MH_data = read_from_path(mh_path, filename_lterm)


# temp means (trials)
def print_global_means(dataset, wiso=None):
    glob_temp = compute_spatial_means(dataset=dataset, varname="temp2", units="°C")
    glob_temp_land = compute_spatial_means(dataset=dataset, varname="temp2", units="°C", land=True)
    glob_temp_ocean = compute_spatial_means(dataset=dataset, varname="temp2", units="°C", ocean=True)
    
    glob_prec = compute_spatial_means(dataset=dataset, varname="prec", units="mm/month")
    glob_prec_land = compute_spatial_means(dataset=dataset, varname="prec", units="mm/month", land=True)
    glob_prec_ocean = compute_spatial_means(dataset=dataset, varname="prec", units="mm/month", ocean=True)
    
    # glob_d18Op = compute_spatial_means(dataset=dataset, varname="d18op", units="per mil", 
    #                                    Dataset_wiso= wiso)
    # glob_d18Op_land = compute_spatial_means(dataset=dataset, varname="d18op", units="per mil",
    #                                         Dataset_wiso= wiso, land=True)
    # glob_d18Op_ocean = compute_spatial_means(dataset=dataset, varname="d18op", units="per mil", 
    #                                          Dataset_wiso= wiso, ocean=True)
    
    print(" Global temperature: {:.2f} , for land: {:.2f} , for ocean: {:.2f} (°C)".format(glob_temp, glob_temp_land,
                                                                                      glob_temp_ocean))
    
    print(" Global precipitation: {:.2f} , for land: {:.2f} , for ocean: {:.2f} (mm/month)".format(glob_prec, glob_prec_land,
                                                                                      glob_prec_ocean))
    
    # print(" Global d18Op: {:.2f} , for land: {:.2f} , for ocean: {:.2f} (per mil)".format(glob_d18Op, glob_d18Op_land,
    #                                                                                   glob_d18Op_ocean))
    
    



print_global_means(dataset=PI_data)
print_global_means(dataset=MH_data)
print_global_means(dataset=LGM_data)
print_global_means(dataset=PLIO_data)
