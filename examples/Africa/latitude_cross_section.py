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
from pyClimat.plots import plot_hovmoller_space_time

# relative imports
from read_data import * # EXP_ID_data for surface variables and EXP_ID_plev_data for pressure variables

path_to_store = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

# define variables
prec = "prec"



def extract_all_variables():
    
    PI_prec = extract_var(Dataset=PI_data , varname=prec, units="mm/month")
    LGM_prec = extract_var(Dataset=LGM_data , varname=prec, units="mm/month")
    MH_prec = extract_var(Dataset=MH_data , varname=prec, units="mm/month")
    PLIO_prec = extract_var(Dataset=PLIO_data , varname=prec, units="mm/month")


    return PI_prec, MH_prec, LGM_prec, PLIO_prec


# Pre_Industrial - other experiments (function)
PI_prec, MH_prec, LGM_prec, PLIO_prec = extract_all_variables()

#PI
PI_prec_alt = compute_lterm_mean(data=PI_prec, time="month")
MH_prec_alt = compute_lterm_mean(data=MH_prec, time="month")
LGM_prec_alt = compute_lterm_mean(data=LGM_prec, time="month")
PLIO_prec_alt = compute_lterm_mean(data=PLIO_prec, time="month")


minlat = 0
maxlat = 30
minlon = -20
maxlon = 30

PI_tp = extract_vertical_section(data=PI_prec_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat",
                                              lev=None)
MH_tp = extract_vertical_section(data=MH_prec_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat",
                                              lev=None)
LGM_tp = extract_vertical_section(data=LGM_prec_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat",
                                              lev=None)
PLIO_tp = extract_vertical_section(data=PLIO_prec_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
 
                                   
 
apply_style(fontsize=22, style=None, linewidth=2)                                   
fig, ((ax1,ax2),(ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(24, 18), sharex=False, sharey=False)

plot_hovmoller_space_time(variable="Precipitation", data=PI_tp, cmap=YlGnBu, units="mm/month", plot_colorbar=True,
                          ax=ax1, bottom_labels=True, left_labels=True, fig=fig, vmax=400, vmin=0,
                                            levels=22, level_ticks=6, cbar_pos=[0.90, 0.35, 0.02, 0.35], title= "[A]  PI")

plot_hovmoller_space_time(variable="Precipitation", data=MH_tp, cmap=YlGnBu, units="mm/month", plot_colorbar=False,
                          ax=ax2, bottom_labels=True, left_labels=True, fig=fig, vmax=400, vmin=0,
                                            levels=22, level_ticks=6, title= "[A]  MH")

plot_hovmoller_space_time(variable="Precipitation", data=LGM_tp, cmap=YlGnBu, units="mm/month", plot_colorbar=False,
                          ax=ax3, bottom_labels=True, left_labels=True, fig=fig, vmax=400, vmin=0,
                                            levels=22, level_ticks=6, title= "[A]  LGM")

plot_hovmoller_space_time(variable="Precipitation", data=PLIO_tp, cmap=YlGnBu, units="mm/month", plot_colorbar=False,
                          ax=ax4, bottom_labels=True, left_labels=True, fig=fig, vmax=400, vmin=0,
                                            levels=22, level_ticks=6, title= "[A]  PLIO")

plt.tight_layout()
plt.subplots_adjust(left=0.02, right=0.86, top=0.98, bottom=0.03)

plt.savefig(os.path.join(path_to_store, "time_space_tp.svg"), format= "svg", bbox_inches="tight", dpi=300)



