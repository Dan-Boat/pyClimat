# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 20:27:38 2024

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
from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"



miocene_path= "D:/Datasets/Model_output_pst/a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"
miocene_ea_high_path = "D:/Datasets/Model_output_pst/a025_dkrz-levante_e5w2.3_t159_MIO_EA_high_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"



# read data
echam_filename = "1003_1017_1m_mlterm.nc"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM



def read_precipitation_uv_data(path_to_data, mpi_varname=None, echam=True, time_range=None):
    
    data = read_from_path(path=path_to_data, filename=echam_filename)
    prec = extract_var(Dataset=data , varname="prec", units="mm/month")
    u = extract_var(Dataset=data , varname="u10")
    v = extract_var(Dataset=data , varname="v10")
        
    return prec, u, v


def compute_anomalies(data_pi, data_exp, time="month", month="JJAS"):
    
    prec_diff_alt = compute_lterm_diff(data_control=data_pi, data_main=data_exp,
                                       time=time, month=month)
    
    return prec_diff_alt



# miocene
mio_prec_data, mio_u_data, mio_v_data = read_precipitation_uv_data(path_to_data=miocene_path)

mio_ea_prec_data, mio_ea_u_data, mio_ea_v_data = read_precipitation_uv_data(path_to_data=miocene_ea_high_path)


mio_prec_diff_j = compute_anomalies(data_pi=mio_prec_data, data_exp=mio_ea_prec_data, 
                             time="month", month="JJAS")

mio_u_diff_j = compute_anomalies(data_pi=mio_u_data, data_exp=mio_ea_u_data, 
                             time="month", month="JJAS")

mio_v_diff_j = compute_anomalies(data_pi=mio_v_data, data_exp=mio_ea_v_data, 
                             time="month", month="JJAS")


mio_prec_diff_m = compute_anomalies(data_pi=mio_prec_data, data_exp=mio_ea_prec_data, 
                             time="month", month="MAM")

mio_u_diff_m = compute_anomalies(data_pi=mio_u_data, data_exp=mio_ea_u_data, 
                             time="month", month="MAM")

mio_v_diff_m = compute_anomalies(data_pi=mio_v_data, data_exp=mio_ea_v_data, 
                             time="month", month="MAM")


# plotting
apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
            
projection = ccrs.Robinson(central_longitude=0, globe=None)

projection_trans = ccrs.PlateCarree()


fig,(ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(8,16), subplot_kw={"projection":projection})

plot_annual_mean(variable="Precipitation anomalies", data_alt=mio_prec_diff_j, ax=ax1,
                 cmap="BrBG", units="mm/month", vmax=100, vmin=-100, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, plot_projection=projection, title="(a) High EARS - CTL (JJAS)", 
                orientation="horizontal",  cbar_pos= [0.20, 0.05, 0.55, 0.02], sea_land_mask=mio_slm,
                domain="West Africa", data_v=mio_v_diff_j, data_u=mio_u_diff_j,
                show_arrow_scale=True, plot_winds=True, wind_scale=40)


plot_annual_mean(variable="Precipitation anomalies", data_alt=mio_prec_diff_m, ax=ax2,
                 cmap="BrBG", units="mm/month [JJAS]", vmax=100, vmin=-100, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, plot_projection=projection, title="(b) High EARS - CTL (MAM)", 
                orientation="horizontal",  cbar_pos= [0.20, 0.05, 0.55, 0.02], sea_land_mask=mio_slm,
                domain="East Africa", data_v=mio_v_diff_m, data_u=mio_u_diff_m,
                show_arrow_scale=True, plot_winds=True, wind_scale=40)


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "poster_map_ears_uplift.pdf"), format= "pdf", bbox_inches="tight", dpi=600)