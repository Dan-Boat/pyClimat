#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:45:48 2023

@author: dboateng

This script generates the plot of precipitation and temperature for the Monsoon region using both ERA5 and ECHAM5-wiso (1980-200)
"""
import os 
import pandas as pd 

from pyClimat.data import read_ERA_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff

# define paths 

path_to_plots = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"
ERA5_path = "/home/dboateng/Datasets/ERA5/monthly_1950_2021/"

ERA5_tp_path = os.path.join(ERA5_path, "tp_monthly.nc")
ERA5_t2m_path = os.path.join(ERA5_path, "t2m_monthly.nc")
# ERA5_v10_path = os.path.join(ERA5_path, "v10_monthly.nc")
# ERA5_u10_path = os.path.join(ERA5_path, "u10_monthly.nc")

#define time of amip 
from1980to2000 = pd.date_range(start="1980-01-01", end="2000-12-31", freq="MS")


#read in postprocessed and analysed data 
ERA5_t2m = read_ERA_processed(path=ERA5_t2m_path, varname="t2m")   - 273.15 #°C
ERA5_tp = read_ERA_processed(path=ERA5_tp_path, varname="tp") * 1000 * 30  #mm/month
# ERA5_v10 = read_ERA_processed(path=ERA5_v10_path, varname="v10") #m/s
# ERA5_u10 = read_ERA_processed(path=ERA5_u10_path, varname="u10") #m/s

ERA5_t2m_alt = compute_lterm_mean(data=ERA5_t2m, time="month", month="JJAS", time_range=from1980to2000)
ERA5_tp_alt = compute_lterm_mean(data=ERA5_tp, time="month", month="JJAS", time_range=from1980to2000)
# ERA5_v10_alt = compute_lterm_mean(data=ERA5_v10, time="month", month="JJAS", time_range=from1980to2000)
# ERA5_u10_alt = compute_lterm_mean(data=ERA5_u10, time="month", month="JJAS", time_range=from1980to2000)



from wam_analysis import PD_t2m_alt, PD_prec_alt, PD_v10_alt, PD_u10_alt


# set up plot function for specific variable

