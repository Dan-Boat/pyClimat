# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 16:38:38 2023

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"


# define paths 
main_path = "D:/Datasets/Model_output_pst/"
lgm_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
plio_path = os.path.join(main_path, "PLIO", "MONTHLY_MEANS")
mh_path = os.path.join(main_path, "MH", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")


filename_lterm = "1003_1017_1m_mlterm.nc"
# read long-term means

PI_data = read_from_path(pi_path, filename_lterm)
LGM_data = read_from_path(lgm_path, filename_lterm)
PLIO_data = read_from_path(plio_path, filename_lterm)
MH_data = read_from_path(mh_path, filename_lterm)


def extract_vars_and_analysis(data, pi_data):
    
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    
    temp2_pi = extract_var(Dataset=pi_data , varname="temp2", units="°C")
    prec_pi = extract_var(Dataset= pi_data , varname="prec", units="mm/month")
    
    
    #compute climatologies difference
    
    temp2_diff = compute_lterm_diff(data_control=temp2_pi, data_main=temp2, time="annual")
    prec_diff = compute_lterm_diff(data_control=prec_pi, data_main=prec, time="annual")
    
    
    
    return_data = {"temperature":temp2_diff, "precipitation":prec_diff}
    
    return return_data


lgm_diff = extract_vars_and_analysis(data=LGM_data, pi_data=PI_data)
mh_diff = extract_vars_and_analysis(data=MH_data, pi_data=PI_data)
mplio_diff = extract_vars_and_analysis(data=PLIO_data, pi_data=PI_data)


apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)

fig,((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(25,18), subplot_kw={"projection":projection})

plot_annual_mean(variable="Precipitation anomalies", data_alt=mh_diff.get("precipitation"), ax=ax1,
                 cmap="BrBG", units="‰", vmax=100, vmin=-100, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(a) MH - PI", 
                orientation="vertical",  cbar_pos= [0.95, 0.65, 0.02, 0.25])

plot_annual_mean(variable="Precipitation anomalies", data_alt=lgm_diff.get("precipitation"), ax=ax2,
                 cmap="BrBG", units="‰", vmax=100, vmin=-100, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) LGM - PI", 
                )

plot_annual_mean(variable="Precipitation anomalies", data_alt=mplio_diff.get("precipitation"), ax=ax3,
                 cmap="BrBG", units="‰", vmax=100, vmin=-100, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(c) mPLIO - PI", 
                )


plot_annual_mean(variable="Temperature anomalies", data_alt=mh_diff.get("temperature"), ax=ax4,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, orientation="vertical",  cbar_pos= [0.95, 0.25, 0.02, 0.25],
                plot_projection=projection, title="(d) MH - PI",)

plot_annual_mean(variable="Temperature anomalies", data_alt=lgm_diff.get("temperature"), ax=ax5,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, 
                plot_projection=projection, title="(e) LGM - PI",)

plot_annual_mean(variable="Temperature anomalies", data_alt=mplio_diff.get("temperature"), ax=ax6,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False,
                plot_projection=projection, title="(f) mPLIO - PI",)


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "global_anomalies.svg"), format= "svg", bbox_inches="tight", dpi=600)