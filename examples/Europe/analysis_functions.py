# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:46:22 2023

@author: dboateng

Note that the method for the EOF can be the standard one from Eof.xarray, or the xefos methods (sklearn),
    or the varimax and promax extentions

Perform all the require analysis (eg. EOFs, sorting of indixes, correlation, regional means)
1. Try the routine for extracting pcs and projecting X fields onto the PD eofs
2. Try the correlation scheme for both one time axis with grids and also two 3D arrays for spearman and pearson, 
get the pvalues and also the uncertainties

"""

import os 
import numpy as np 
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt 

from pyClimat.analysis import extract_var
from pyClimat.stats import EOF_standard
from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance

# analysis for ERA5 Dataset

from read_data import ERA5_msl

# initiate the eof instance 
ERA5_EOF = EOF_standard(data=ERA5_msl, weights=True, standardize=True, 
                      extract_region=True, extract_season=True, neofs=4)

# select the region of interest and season
ERA5_EOF.select_time_and_region(maxlon=60, minlon=-80, maxlat=80, minlat=20, time="season", 
                              season="DJF")

# calculate the anomalies and apply norm
ERA5_EOF.calculate_anomalies()

method = "Eof"

# fit the eof with the solver
ERA5_EOF.eof_solver(method=method, apply_promax=False, apply_varimax=False)

# extract the eofs (as the eigenvectors of the covariance matric of the X field)

eofs = ERA5_EOF.eofs(eofscaling=2) # 2 - multiply by the singular values or square root of the eigen values
pcs = ERA5_EOF.pcs(pscaling=1)
variance_ratio = ERA5_EOF.explained_variance_ratio()

# visualize
apply_style(fontsize=22, style=None, linewidth=2) 
projection = ccrs.PlateCarree()

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols=2, 
                                             figsize=(15, 8), subplot_kw={"projection": projection})


plot_eofsAsCovariance(variable= "slp", data=eofs.sel(mode=1), mode_var=variance_ratio[0], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=False,
                      ax=ax1, fig=fig, title="[A]", bottom_labels=False)

plt.show()
# analysis for echam PD


# analysis of projecting PI to PD eofs



# from read_data import PD_data, LGM_data, PI_data

# #extract variable from data 
# PD_slp = extract_var(PD_data, "slp", units="hPa")
# LGM_slp = extract_var(LGM_data, "slp", units="hPa")
# PI_slp = extract_var(PI_data, "slp", units="hPa")


# #applying EOF class

# PD_EOF = EOF_standard(data=PD_slp, weights=True, standardize=True, 
#                       extract_region=True, extract_season=True, neofs=4)

# PD_EOF.select_time_and_region(maxlon=60, minlon=-80, maxlat=80, minlat=20, time="season", 
#                               season="DJF", month="ONDJFM")
# PD_EOF.calculate_anomalies()
# #PD_EOF.eof_solver(method="xeofs", apply_promax=True)
# PD_EOF.eof_solver(method="Eof", apply_promax=True)
# eofs = PD_EOF.eofs(eofscaling=2)
# pcs = PD_EOF.pcs(pscaling=1)