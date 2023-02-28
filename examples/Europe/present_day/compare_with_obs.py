# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 20:10:19 2023

@author: dboateng

This script aim to compare the NAO and EA from ERA5 and ECHAM to obs records
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdate

from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance

#defining paths
main_path_to_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/PD"
main_path_to_obs = "C:/Users/dboateng/Desktop/Datasets/NAO"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"


EA_valencia_DJF_path = os.path.join(main_path_to_obs, "EA_Valencia_DJF")
EA_valencia_JJA_path = os.path.join(main_path_to_obs, "EA_Valencia_JJA")

NAO_Gilbraltar_JJA_path = os.path.join(main_path_to_obs, "NAO_Gilbraltar_JJA.csv")
NAO_Gilbraltar_DJF_path = os.path.join(main_path_to_obs, "NAO_Gilbraltar_DJF.csv")

NAO_CDC_JJA_path = os.path.join(main_path_to_obs, "NAO_CDC_JJA.csv")
NAO_CDC_DJF_path = os.path.join(main_path_to_obs, "NAO_CDC_DJF.csv")

EOFs_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_eofsAsCovariance.nc") # standard method
EOFs_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_eofsAsCovariance.nc")

PCs_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_pcs.csv") # standard method
PCs_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_pcs.csv")

Vars_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_pcs_variance.csv") # standard method
Vars_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_pcs_variance.csv")


EOFs_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_eofsAsCovariance.nc") # standard method
EOFs_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_eofsAsCovariance.nc")

PCs_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_pcs.csv") # standard method
PCs_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_pcs.csv")

Vars_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_pcs_variance.csv") # standard method
Vars_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_pcs_variance.csv")



# read all the required datasets (start with winter)

# read EOFs and plot with thier variance

eof_era_djf_data = xr.open_dataarray(EOFs_ERA_DJF_path)
eof_echam_djf_data = xr.open_dataarray(EOFs_ECHAM_DJF_path)
variance_era_djf = pd.read_csv(Vars_ERA_DJF_path, index_col=["mode"])
variance_echam_djf = pd.read_csv(Vars_ECHAM_DJF_path, index_col=["mode"])

# plot eofs

apply_style(fontsize=22, style=None, linewidth=2)
projection = ccrs.AlbersEqualArea(
    central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))

fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols=2, 
                                             figsize=(18, 15), subplot_kw={"projection": projection})

modes = [1,2]
units="hPa" 
variable="Mean Sea Level Pressure"
vmax=15
vmin=-15
figname ="NAO_era_echam_obs"
plot_eofsAsCovariance(variable= variable, data=eof_era_djf_data.sel(mode=modes[0]) * -1, mode_var=variance_era_djf.loc[modes[0]], units=units, vmax=vmax,
                      vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                      ax=ax1, fig=fig, title="[A] ERA5 (1959-2021)", bottom_labels=True)

plot_eofsAsCovariance(variable= variable, data=eof_echam_djf_data.sel(mode=modes[0]), mode_var=variance_echam_djf.loc[modes[0]], units=units, 
                      vmax=vmax, vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] ECHAM5-wiso (1979-2014)", bottom_labels=True, use_AlberEqualArea=True)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, figname + ".svg"), format= "svg", bbox_inches="tight", dpi=300)
