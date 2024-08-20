# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:55:38 2023

@author: dboateng
"""
import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance

from path_to_data_lm import *

main_path_to_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots/data"
path_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"


#path to covariance
CESM_DJF_EOFs = os.path.join(main_path_to_data, "CESM_DJF_eofsAsCorrelation.nc")
GISS_DJF_EOFs = os.path.join(main_path_to_data, "GISS_DJF_eofsAsCorrelation.nc")
ECHAM_DJF_EOFs = os.path.join(main_path_to_data, "ECHAM5_DJF_eofsAsCorrelation.nc")



# PCS
CESM_DJF_VARs = os.path.join(main_path_to_data, "CESM_DJF_pcs_variance.csv")
ECHAM_DJF_VARs = os.path.join(main_path_to_data, "ECHAM5_DJF_pcs_variance.csv")
GISS_DJF_VARs = os.path.join(main_path_to_data, "GISS_DJF_pcs_variance.csv")




def read_eof_var(eof_path, vars_path):
    eof = xr.open_dataarray(eof_path)
    vars_data = pd.read_csv(vars_path, index_col=["mode"])
    
    return eof, vars_data 


cesm_eof, cesm_vars = read_eof_var(eof_path=CESM_DJF_EOFs, vars_path=CESM_DJF_VARs)
giss_eof, giss_vars = read_eof_var(eof_path=GISS_DJF_EOFs , vars_path=GISS_DJF_VARs)
echam_eof, echam_vars = read_eof_var(eof_path=ECHAM_DJF_EOFs , vars_path=ECHAM_DJF_VARs)

apply_style(fontsize=23, style="seaborn-talk", linewidth=3,) 
#apply_style(fontsize=22, style=None, linewidth=2)
projection = ccrs.AlbersEqualArea(
    central_latitude=35, central_longitude=-35, standard_parallels=(20, 80))

fig, ((ax1,ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(18, 22), 
                                                                           subplot_kw={"projection": projection})

units="-" 
variable="Correlation"
vmax=1
vmin=-1

figname ="NAO_EA_spatial_fields"


plot_eofsAsCovariance(variable= variable, data=cesm_eof.sel(mode=1) * -1, mode_var=cesm_vars.loc[1], units=units, vmax=vmax,
                      vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.60, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=True,
                      ax=ax1, fig=fig, title="(a) iCESM  [NAO]", bottom_labels=False)

plot_eofsAsCovariance(variable= variable, data=giss_eof.sel(mode=1) * -1, mode_var=giss_vars.loc[1], units=units, vmax=vmax,
                      vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=True,
                      ax=ax3, fig=fig, title=" (c) GISS-E2-R", bottom_labels=False)

plot_eofsAsCovariance(variable= variable, data=echam_eof.sel(mode=1) * -1, mode_var=echam_vars.loc[1], units=units, vmax=vmax,
                      vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=True,
                      ax=ax5, fig=fig, title="(e) ECHAM5-wiso", bottom_labels=True)

plot_eofsAsCovariance(variable= variable, data=cesm_eof.sel(mode=3), mode_var=cesm_vars.loc[3], units=units, vmax=vmax,
                      vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=True,
                      ax=ax2, fig=fig, title=" (b )iCESM [EA]", bottom_labels=False)

plot_eofsAsCovariance(variable= variable, data=giss_eof.sel(mode=3), mode_var=giss_vars.loc[3], units=units, vmax=vmax,
                      vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=True,
                      ax=ax4, fig=fig, title=" (c) GISS-E2-R", bottom_labels=False)

plot_eofsAsCovariance(variable= variable, data=echam_eof.sel(mode=4), mode_var=echam_vars.loc[4], units=units, vmax=vmax,
                      vmin=vmin, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=True,
                      ax=ax6, fig=fig, title="(f) ECHAM5-wiso", bottom_labels=True)



fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.98, top=0.98, bottom=0.10,)
plt.savefig(os.path.join(path_plots, figname + ".pdf"), format= "pdf", bbox_inches="tight", dpi=600)
