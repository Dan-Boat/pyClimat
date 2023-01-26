#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 17:00:52 2023

@author: dboateng
"""
import os 
import xarray as xr
import pandas as pd
from pyClimat.data import read_from_path, read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean
from pyClimat.plots import plot_annual_mean
from pyClimat.plot_utils import *

path_to_store = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"
path_cru = "/media/dboateng/Boateng/Datasets/CRU"
path_chrip = "/media/dboateng/Boateng/Datasets/CHIRPS"

main_path = "/home/dboateng/Model_output_pst"
exp_name = "t004_dkrz-mistral_e5w2.3_AMIP_t159l31.6h"    # simulation with present-day simulation (not different from PI simulations)

years = "1980_2000"
period = "1m"

PD_data = read_ECHAM_processed(main_path=main_path , exp_name= exp_name, years=years,
                                                  period=period, read_wiso=False)

echam_prec = extract_var(Dataset=PD_data, varname="prec", units="mm/month")

#define time of amip 
from1901to2021 = pd.date_range(start="1901-01-01", end="2021-12-31", freq="MS")
from1981to2022 = pd.date_range(start="1981-01-01", end="2022-11-30", freq="MS")
from1981to2000 = pd.date_range(start="1981-01-01", end="2000-12-31", freq="MS")

data_cru = read_from_path(path= path_cru, filename="pr_monthly.nc", decode=False, varname="pre")
data_chirp = read_from_path(path= path_chrip, filename="pr_monthly.nc", decode=True, varname="precip")
data_cru["time"] = from1901to2021
data_chirp["time"] = from1981to2022

echam_prec_alt = compute_lterm_mean(data=echam_prec, time="month", month="JJAS")
cru_prec_alt = compute_lterm_mean(data=data_cru, time="month", month="JJAS", time_range=from1981to2000)
cru_data = xr.where(cru_prec_alt > 1000, 0, cru_prec_alt)  # cut continent

chirp_prec_alt = compute_lterm_mean(data=data_chirp, time="month", month="JJAS", time_range=from1981to2000)

apply_style(fontsize=22, style=None, linewidth=2) 

projection = ccrs.PlateCarree()
fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(15, 13), subplot_kw={"projection": projection})


plot_annual_mean(ax=ax1, fig=fig, variable="Precipitation", data_alt=echam_prec_alt, cmap=YlGnBu, units="mm/month", vmax=400, vmin=50, domain="West Africa", 
                  levels=22, level_ticks=6, title="ECHAM5-wiso", left_labels=True, bottom_labels=True, 
                  add_colorbar=True, cbar_pos = [0.35, 0.01, 0.25, 0.02], orientation= "horizontal")

plot_annual_mean(ax=ax2, fig=fig, variable="Precipitation", data_alt=cru_data, cmap=YlGnBu, units="mm/month", vmax=400, vmin=50, domain="West Africa", 
                  levels=22, level_ticks=6, title="CRU", left_labels=False, bottom_labels=True, 
                  add_colorbar=False)
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_plots, "compare_cru_and_echam.svg"), format= "svg", bbox_inches="tight", dpi=300)