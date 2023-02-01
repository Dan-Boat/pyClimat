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


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"
# analysis for ERA5 Dataset

from read_data import ERA5_msl

# initiate the eof instance 
ERA5_EOF = EOF_standard(data=ERA5_msl, weights=True, standardize=True, 
                      extract_region=True, extract_season=True, neofs=4)

# select the region of interest and season
ERA5_EOF.select_time_and_region(maxlon=60, minlon=-80, maxlat=80, minlat=20, time="season", 
                              season="JJA", month="JA") # month="AMJJAS"

# calculate the anomalies and apply norm
ERA5_EOF.calculate_anomalies()

method = "xeofs"

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
                                             figsize=(24, 13), subplot_kw={"projection": projection})

# loop through this !!

plot_eofsAsCovariance(variable= "slp", data=eofs.sel(mode=1), mode_var=variance_ratio[1], units="hPa", vmax=15, vmin=-15, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=False,
                      ax=ax1, fig=fig, title="[A]", bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=eofs.sel(mode=2), mode_var=variance_ratio[2], units="hPa", vmax=15, vmin=-15, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B]", bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=eofs.sel(mode=3), mode_var=variance_ratio[3], units="hPa", vmax=15, vmin=-15, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C]", bottom_labels=True)

plot_eofsAsCovariance(variable= "slp", data=eofs.sel(mode=4), mode_var=variance_ratio[4], units="hPa", vmax=15, vmin=-15, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D]", bottom_labels=True)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, "sNAO_EOF_standard.svg"), format= "svg", bbox_inches="tight", dpi=300)
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