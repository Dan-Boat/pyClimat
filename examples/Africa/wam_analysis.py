#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:29:38 2022

@author: dboateng
extract variables and compute seasonal and monthly means using the read files from ./read_data.py
Pivot into functions to help run specific task as a time
"""

import os 
import xarray as xr
import pandas as pd 
import numpy as np 

#pyClimat models
from pyClimat.analysis import extract_var, compute_lterm_mean, extract_transect

# relative imports
from read_data import * # EXP_ID_data for surface variables and EXP_ID_plev_data for pressure variables

# define variables
t2m = "temp2"
prec = "prec"
v10 = "v10"
u10 = "u10"

# Pre_Industrial
PI_t2m = extract_var(Dataset=PI_data, varname=t2m, units="°C")
PI_prec = extract_var(Dataset=PI_data, varname=prec, units="mm/month")
PI_v10 = extract_var(Dataset=PI_data, varname=v10) # default in m/s
PI_u10 = extract_var(Dataset=PI_data, varname=u10) # default in m/s
PI_slp = extract_var(Dataset=PI_plev_data, varname="slp", units="hPa")
PI_v = extract_var(Dataset=PI_plev_data, varname="v", lev_units="hPa")
PI_u = extract_var(Dataset=PI_plev_data, varname="u", lev_units="hPa")
PI_omega = extract_var(Dataset=PI_plev_data, varname="omega", lev_units="hPa")





# select the monsoon months
PI_t2m_alt = compute_lterm_mean(data=PI_t2m, time="month", month="JJAS")
PI_prec_alt = compute_lterm_mean(data=PI_prec, time="month", month="JJAS")
PI_v10_alt = compute_lterm_mean(data=PI_v10, time="month", month="JJAS")
PI_u10_alt = compute_lterm_mean(data=PI_u10, time="month", month="JJAS")
PI_v_alt = compute_lterm_mean(data=PI_v, time="month", month="JJAS")
PI_u_alt = compute_lterm_mean(data=PI_u, time="month", month="JJAS")
PI_omega_alt = compute_lterm_mean(data=PI_omega, time="month", month="JJAS")


# extract the sections for monthly variability (Sahel: 10-20 N, 20W - 30E; coast of Guinea: 5-10N, 20W-30E, Sahara region: 20-30N, 20W-30E)

def weighted_mean(data):
    weights = np.cos(np.deg2rad(data.lat))
    weights.name = "weights"
    
    data_weighted = data.weighted(weights)
    
    data_mean = data_weighted.mean(dim=("lat", "lon"), skipna=True)
    return data_mean 

def weighted_std(data):
    weights = np.cos(np.deg2rad(data.lat))
    weights.name = "weights"
    
    data_weighted = data.weighted(weights)
    
    data_std = data_weighted.std(dim=("lat", "lon"), skipna=True)
    return data_std 


def compute_mean_std(data, save=False, path_to_save=None):
    data_mean = weighted_mean(data)
    data_std = weighted_std(data)
    
    import calendar
    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    df = pd.DataFrame(index = mnames, columns = ["mean","std"])
    df["mean"] = data_mean
    df["std"] = data_std
    if save == True:
        df.to_csv(path_to_save)
    return df


def extract_sections(data):
    
    max_lon, max_lat, min_lon, min_lat = 30, 20, -20, 10
    
    extract_sahel = extract_transect(data=data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat)
    
    df_sahel = compute_mean_std(extract_sahel)
    
    max_lon, max_lat, min_lon, min_lat = 30, 10, -20, 5
    
    extract_guinea = extract_transect(data=data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat)
    
    df_guinea = compute_mean_std(extract_guinea)
    
    max_lon, max_lat, min_lon, min_lat = 30, 30, -20, 20
    
    extract_sahara = extract_transect(data=data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat)
    
    df_sahara = compute_mean_std(extract_sahara)
    
    return df_sahara, df_sahel, df_guinea


 
PI_month_sahara_prec, PI_month_sahel_prec, PI_month_guinea_prec = extract_sections(data=PI_prec)
PI_month_sahara_t2m, PI_month_sahel_t2m, PI_month_guinea_t2m = extract_sections(data=PI_t2m)

# extract vertical section to show the ATJ, EAJ, WAM surface westerlies, and updraft and subsidence

