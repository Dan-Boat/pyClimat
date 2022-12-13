#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 19:26:47 2022

@author: dboateng
"""

import os 
import xarray as xr
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt

#pyClimat models
from pyClimat.analysis import extract_var, compute_lterm_mean, extract_transect, compute_lterm_diff
from pyClimat.analysis import extract_vertical_section
from pyClimat.plot_utils import *

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


# Pre_Industrial - other experiments (function)
PI_t2m, PI_prec, PI_v10, PI_u10, PI_slp, PI_v, PI_u, PI_omega, PI_v850, PI_u850 = extract_all_variables(data_surface=PI_data, 
                                                                                             data_plev=PI_plev_data)

#PI
PI_t2m_alt = compute_lterm_mean(data=PI_t2m, time="month")
PI_prec_alt = compute_lterm_mean(data=PI_prec, time="month")
PI_v_alt = compute_lterm_mean(data=PI_v, time="month")
PI_u_alt = compute_lterm_mean(data=PI_u, time="month")
PI_omega_alt = compute_lterm_mean(data=PI_omega, time="month")
PI_slp_alt = compute_lterm_mean(data=PI_slp, time="month")
PI_u850_alt = compute_lterm_mean(data=PI_u850, time="month")
PI_v850_alt = compute_lterm_mean(data=PI_v850, time="month")


minlat = 0
maxlat = 30
minlon = -20
maxlon = 30

PI_cross_section_prec = extract_vertical_section(data=PI_prec_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat",
                                              lev=None)

data = PI_cross_section_prec
y = data.index.values
x = data.columns.values
     
X,Y = np.meshgrid(x,y)
Z = data.values

fig, ax = plt.subplots(nrows=1, ncols=1, sharex= False, sharey= False, figsize=(15, 10))

p = ax.contourf(X,Y,Z, cmap=YlGnBu)
c = ax.contour(X,Y,Z,color="black", linewidth=2, levels=10)
clb = ax.clabel(c, fmt="%2.0f", use_clabeltext=True)

plt.show()