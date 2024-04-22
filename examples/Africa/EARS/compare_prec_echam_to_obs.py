# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 12:31:12 2024

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl 
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_from_path, read_ERA_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"


# define paths
pd_echam_path = "D:/Datasets/Model_output_pst/PD/"
era_path = "D:/Datasets/ERA5/monthly_1950_2021/tp_monthly.nc"
cru_path = "E:/Datasets/CRU/"

from1901to2021 = pd.date_range(start="1901-01-01", end="2021-12-31", freq="MS")
from1981to2000 = pd.date_range(start="1981-01-01", end="2000-12-31", freq="MS")
from1980to2014 = pd.date_range(start="1980-01-01", end="2014-12-31", freq="MS")
#read data

data_echam = read_from_path(path=pd_echam_path, filename="PD_1980_2014_monthly.nc")
data_cru = read_from_path(path= cru_path, filename="pr_monthly.nc", decode=False, varname="pre")
data_cru["time"] = from1901to2021

ERA5_data = read_ERA_processed(path=era_path, varname="tp") * 1000 * 30  #mm/month




# compute long-term means
echam_prec = extract_var(Dataset=data_echam, varname="prec", units="mm/month")
echam_prec_alt = compute_lterm_mean(data=echam_prec, time="annual")

cru_prec_alt = compute_lterm_mean(data=data_cru, time="annual", time_range=from1981to2000)
cru_prec_alt = xr.where(cru_prec_alt > 1000, 0, cru_prec_alt)  # cut continent


era_prec_alt = compute_lterm_mean(data=ERA5_data, time="annual", time_range=from1980to2014)


# plotting
apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
            
projection = ccrs.Robinson(central_longitude=0, globe=None)

pprojection_trans = ccrs.PlateCarree()

fig,(ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(18,13), subplot_kw={"projection":projection})
plot_annual_mean(variable="Precipitation", data_alt=echam_prec_alt, ax=ax1,
                 cmap="YlGnBu", units="mm/month", vmax=300, vmin=50, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(a) ECHAM5-wiso", 
                domain="Africa", cbar_pos= [0.35, 0.01, 0.25, 0.02], center=False,
                orientation="horizontal")


plot_annual_mean(variable="Precipitation", data_alt=era_prec_alt, ax=ax2,
                 cmap="YlGnBu", units="mm/month", vmax=300, vmin=50, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) ERA5", 
                domain="Africa", center=False)

plot_annual_mean(variable="Precipitation", data_alt=cru_prec_alt, ax=ax3,
                 cmap="YlGnBu", units="mm/month", vmax=300, vmin=50, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) CRU", 
                domain="Africa", center=False)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "poster_fig2.pdf"), format= "pdf", bbox_inches="tight", dpi=600)

