# -*- coding: utf-8 -*-
"""
Created on Sat May 20 17:27:36 2023

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
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, compute_lterm_diff


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

# experiments
W2E1_Mio278_filename = "a017_hpc-bw_e5w2.3_t159_MIO_W2E1_278ppm_t159l31.6h"
W2E1_Mio450_filename = "a016_hpc-bw_e5w2.3_t159_MIO_W2E1_450ppm_t159l31.6h"

W2E15_Mio278_filename = "a023_dkrz-levante_e5w2.3_t159_MIO_W2E1.5_278ppm_t159l31.6h"
W2E15_Mio450_filename = "a022_hpc-bw_e5w2.3_t159_MIO_W2E1.5_450ppm_t159l31.6h"


W2E0_Mio278_filename="a019_hpc-bw_e5w2.3_t159_MIO_W2E0_278ppm_t159l31.6h"
W2E0_Mio450_filename="a018_hpc-bw_e5w2.3_t159_MIO_W2E0_450ppm_t159l31.6h"

W2E2_Mio278_filename="a020_dkrz-levante_e5w2.3_t159_MIO_W2E2_278ppm_t159l31.6h"
W2E2_Mio450_filename="a021_dkrz-levante_e5w2.3_t159_MIO_W2E2_450ppm_t159l31.6h"


# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E1_278_data, W2E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_450_data, W2E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)

W2E15_278_data, W2E15_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E15_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E15_450_data, W2E15_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E15_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E0_278_data, W2E0_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_Mio278_filename, 
                                                    years=years, period=period, read_wiso=True)

W2E0_450_data, W2E0_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E2_278_data, W2E2_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio278_filename, 
                                                    years=years, period=period, read_wiso=True)

W2E2_450_data, W2E2_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


def extract_qvi_and_analysis(exp_data, exp_wiso, W1E1_data, W1E1_wiso):
    
    # extract prec
    exp_qvi_data = extract_var(Dataset=exp_data, varname="apmeb")
    w1e1_qvi_data = extract_var(Dataset=W1E1_data, varname="qvi")
    
    simulated_qvi_change_alt = compute_lterm_diff(data_control=w1e1_qvi_data, data_main=exp_qvi_data,
                                          time="annual")
    
    
    return_data_monthly = {"control_mon": w1e1_qvi_data, "simulated_mon": exp_qvi_data,}
    return_data_ltm = {"simulated_change_ltm": simulated_qvi_change_alt}
    
    return_data = return_data_monthly | return_data_ltm
    
    
    return return_data

W2E0_278_diff = extract_qvi_and_analysis(exp_data=W2E0_278_data, exp_wiso=W2E0_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso)