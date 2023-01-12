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
from pyClimat.analysis import extract_var, compute_lterm_mean, extract_transect, compute_lterm_diff
from pyClimat.analysis import extract_vertical_section

# relative imports
from read_data import * # EXP_ID_data for surface variables and EXP_ID_plev_data for pressure variables

# define variables
t2m = "temp2"
prec = "prec"
v10 = "v10"
u10 = "u10"


def extract_all_variables(data_surface, data_plev):
    
    data_t2m = extract_var(Dataset=data_surface, varname=t2m, units="°C")
    data_prec = extract_var(Dataset=data_surface, varname=prec, units="mm/month")
    data_v10 = extract_var(Dataset=data_surface, varname=v10) # default in m/s
    data_u10 = extract_var(Dataset=data_surface, varname=u10) # default in m/s
    data_slp = extract_var(Dataset=data_plev, varname="slp", units="hPa")
    data_v = extract_var(Dataset=data_plev, varname="v", lev_units="hPa")
    data_u = extract_var(Dataset=data_plev, varname="u", lev_units="hPa")
    data_omega = extract_var(Dataset=data_plev, varname="omega", lev_units="hPa")
    data_v850 = data_v.sel(lev=850)
    data_u850 = data_u.sel(lev=850)

    return data_t2m, data_prec, data_v10, data_u10, data_slp, data_v, data_u, data_omega, data_v850, data_u850

def extract_only_prec_temp_winds(data_surface):
    
    data_t2m = extract_var(Dataset=data_surface, varname=t2m, units="°C")
    
    data_prec = extract_var(Dataset=data_surface, varname=prec, units="mm/month")
    
    data_v10 = extract_var(Dataset=data_surface, varname=v10) # default in m/s
    data_u10 = extract_var(Dataset=data_surface, varname=u10) # default in m/s
    
    
    return data_prec, data_t2m, data_v10, data_u10


# Pre_Industrial - other experiments (function)
PI_t2m, PI_prec, PI_v10, PI_u10, PI_slp, PI_v, PI_u, PI_omega, PI_v850, PI_u850 = extract_all_variables(data_surface=PI_data, 
                                                                                                                            data_plev=PI_plev_data)

LGM_t2m, LGM_prec, LGM_v10, LGM_u10, LGM_slp, LGM_v, LGM_u, LGM_omega, LGM_v850, LGM_u850 = extract_all_variables(data_surface=LGM_data, 
                                                                                                                            data_plev=LGM_plev_data)

MH_t2m, MH_prec, MH_v10, MH_u10, MH_slp, MH_v, MH_u, MH_omega, MH_v850, MH_u850 = extract_all_variables(data_surface=MH_data, 
                                                                                                                            data_plev=MH_plev_data)

PLIO_t2m, PLIO_prec, PLIO_v10, PLIO_u10, PLIO_slp, PLIO_v, PLIO_u, PLIO_omega, PLIO_v850, PLIO_u850 = extract_all_variables(data_surface=PLIO_data, 
                                                                                                                            data_plev=PLIO_plev_data)

PD_prec, PD_t2m, PD_v10, PD_u10 = extract_only_prec_temp(PD_data)


# select the monsoon months then estimate annual means
#PD

PD_t2m_alt = compute_lterm_mean(data=PD_t2m, time="month", month="JJAS")
PD_prec_alt = compute_lterm_mean(data=PD_prec, time="month", month="JJAS")
PD_v10_alt = compute_lterm_mean(data=PD_v10, time="month", month="JJAS")
PD_u10_alt = compute_lterm_mean(data=PD_u10, time="month", month="JJAS")

#PI
PI_t2m_alt = compute_lterm_mean(data=PI_t2m, time="month", month="JJAS")
PI_prec_alt = compute_lterm_mean(data=PI_prec, time="month", month="JJAS")
PI_v10_alt = compute_lterm_mean(data=PI_v10, time="month", month="JJAS")
PI_u10_alt = compute_lterm_mean(data=PI_u10, time="month", month="JJAS")
PI_v_alt = compute_lterm_mean(data=PI_v, time="month", month="JJAS")
PI_u_alt = compute_lterm_mean(data=PI_u, time="month", month="JJAS")
PI_omega_alt = compute_lterm_mean(data=PI_omega, time="month", month="JJAS")
PI_slp_alt = compute_lterm_mean(data=PI_slp, time="month", month="JJAS")
PI_u850_alt = compute_lterm_mean(data=PI_u850, time="month", month="JJAS")
PI_v850_alt = compute_lterm_mean(data=PI_v850, time="month", month="JJAS")

#LGM
LGM_t2m_alt = compute_lterm_mean(data=LGM_t2m, time="month", month="JJAS")
LGM_prec_alt = compute_lterm_mean(data=LGM_prec, time="month", month="JJAS")
LGM_v10_alt = compute_lterm_mean(data=LGM_v10, time="month", month="JJAS")
LGM_u10_alt = compute_lterm_mean(data=LGM_u10, time="month", month="JJAS")
LGM_v_alt = compute_lterm_mean(data=LGM_v, time="month", month="JJAS")
LGM_u_alt = compute_lterm_mean(data=LGM_u, time="month", month="JJAS")
LGM_omega_alt = compute_lterm_mean(data=LGM_omega, time="month", month="JJAS")
LGM_slp_alt = compute_lterm_mean(data=LGM_slp, time="month", month="JJAS")
LGM_u850_alt = compute_lterm_mean(data=LGM_u850, time="month", month="JJAS")
LGM_v850_alt = compute_lterm_mean(data=LGM_v850, time="month", month="JJAS")

#MH
MH_t2m_alt = compute_lterm_mean(data=MH_t2m, time="month", month="JJAS")
MH_prec_alt = compute_lterm_mean(data=MH_prec, time="month", month="JJAS")
MH_v10_alt = compute_lterm_mean(data=MH_v10, time="month", month="JJAS")
MH_u10_alt = compute_lterm_mean(data=MH_u10, time="month", month="JJAS")
MH_v_alt = compute_lterm_mean(data=MH_v, time="month", month="JJAS")
MH_u_alt = compute_lterm_mean(data=MH_u, time="month", month="JJAS")
MH_omega_alt = compute_lterm_mean(data=MH_omega, time="month", month="JJAS")
MH_slp_alt = compute_lterm_mean(data=MH_slp, time="month", month="JJAS")
MH_u850_alt = compute_lterm_mean(data=MH_u850, time="month", month="JJAS")
MH_v850_alt = compute_lterm_mean(data=MH_v850, time="month", month="JJAS")

#PLIO
PLIO_t2m_alt = compute_lterm_mean(data=PLIO_t2m, time="month", month="JJAS")
PLIO_prec_alt = compute_lterm_mean(data=PLIO_prec, time="month", month="JJAS")
PLIO_v10_alt = compute_lterm_mean(data=PLIO_v10, time="month", month="JJAS")
PLIO_u10_alt = compute_lterm_mean(data=PLIO_u10, time="month", month="JJAS")
PLIO_v_alt = compute_lterm_mean(data=PLIO_v, time="month", month="JJAS")
PLIO_u_alt = compute_lterm_mean(data=PLIO_u, time="month", month="JJAS")
PLIO_omega_alt = compute_lterm_mean(data=PLIO_omega, time="month", month="JJAS")
PLIO_slp_alt = compute_lterm_mean(data=PLIO_slp, time="month", month="JJAS")
PLIO_u850_alt = compute_lterm_mean(data=PLIO_u850, time="month", month="JJAS")
PLIO_v850_alt = compute_lterm_mean(data=PLIO_v850, time="month", month="JJAS")

# compute the long-term mean difference



LGM_t2m_alt_diff = compute_lterm_diff(data_control=PI_t2m, data_main= LGM_t2m, time="month", month="JJAS")
LGM_prec_alt_diff = compute_lterm_diff(data_control=PI_prec, data_main= LGM_prec, time="month", month="JJAS")
LGM_slp_alt_diff = compute_lterm_diff(data_control=PI_slp, data_main=LGM_slp, time="month", month="JJAS")

MH_t2m_alt_diff = compute_lterm_diff(data_control=PI_t2m, data_main= MH_t2m, time="month", month="JJAS")
MH_prec_alt_diff = compute_lterm_diff(data_control=PI_prec, data_main= MH_prec, time="month", month="JJAS")
MH_slp_alt_diff = compute_lterm_diff(data_control=PI_slp, data_main=MH_slp, time="month", month="JJAS")


PLIO_t2m_alt_diff = compute_lterm_diff(data_control=PI_t2m, data_main= PLIO_t2m, time="month", month="JJAS")
PLIO_prec_alt_diff = compute_lterm_diff(data_control=PI_prec, data_main= PLIO_prec, time="month", month="JJAS")
PLIO_slp_alt_diff = compute_lterm_diff(data_control=PI_slp, data_main=PLIO_slp, time="month", month="JJAS")


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
LGM_month_sahara_prec, LGM_month_sahel_prec, LGM_month_guinea_prec = extract_sections(data=LGM_prec)
MH_month_sahara_prec, MH_month_sahel_prec, MH_month_guinea_prec = extract_sections(data=MH_prec)
PLIO_month_sahara_prec, PLIO_month_sahel_prec, PLIO_month_guinea_prec = extract_sections(data=PLIO_prec)


PI_month_sahara_t2m, PI_month_sahel_t2m, PI_month_guinea_t2m = extract_sections(data=PI_t2m)
LGM_month_sahara_t2m, LGM_month_sahel_t2m, LGM_month_guinea_t2m = extract_sections(data=LGM_t2m)
MH_month_sahara_t2m, MH_month_sahel_t2m, MH_month_guinea_t2m = extract_sections(data=MH_t2m)
PLIO_month_sahara_t2m, PLIO_month_sahel_t2m, PLIO_month_guinea_t2m = extract_sections(data=PLIO_t2m)

# extract vertical section to show the ATJ, EAJ, WAM surface westerlies, and updraft and subsidence (zonal, meridoinal, and omega)

minlat = 0
maxlat = 30
minlon = -20
maxlon = 30

PI_cross_section_u = extract_vertical_section(data=PI_u_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
MH_cross_section_u = extract_vertical_section(data=MH_u_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
LGM_cross_section_u = extract_vertical_section(data=LGM_u_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
PLIO_cross_section_u = extract_vertical_section(data=PLIO_u_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")

PI_cross_section_v = extract_vertical_section(data=PI_v_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
MH_cross_section_v = extract_vertical_section(data=MH_v_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
LGM_cross_section_v = extract_vertical_section(data=LGM_v_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
PLIO_cross_section_v = extract_vertical_section(data=PLIO_v_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")

PI_cross_section_omega = extract_vertical_section(data=PI_omega_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
MH_cross_section_omega = extract_vertical_section(data=MH_omega_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
LGM_cross_section_omega = extract_vertical_section(data=LGM_omega_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
PLIO_cross_section_omega = extract_vertical_section(data=PLIO_omega_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")


MH_PI_u_cross_section = MH_cross_section_u - PI_cross_section_u
LGM_PI_u_cross_section = LGM_cross_section_u - PI_cross_section_u
PLIO_PI_u_cross_section = PLIO_cross_section_u - PI_cross_section_u