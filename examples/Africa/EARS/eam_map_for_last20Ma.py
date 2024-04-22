# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 22:42:30 2024

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


# define paths
pi_path = "D:/Datasets/Model_output_pst/PI/MONTHLY_MEANS/"
miocene_path= "D:/Datasets/Model_output_pst/a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"
pliocene_path= "D:/Datasets/Model_output_pst/PLIO/MONTHLY_MEANS/"
holocene_path = "D:/Datasets/Model_output_pst/MH/MONTHLY_MEANS/"
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


def compute_anomalies(data_pi, data_exp, time="month", month="MAM"):
    
    prec_diff_alt = compute_lterm_diff(data_control=data_pi, data_main=data_exp,
                                       time=time, month=month)
    
    return prec_diff_alt


#pi 
pi_prec_data, pi_u_data, pi_v_data = read_precipitation_uv_data(path_to_data=pi_path)


pi_prec_alt = compute_lterm_mean(data=pi_prec_data, time="month", month="MAM")


# miocene
mio_prec_data, mio_u_data, mio_v_data = read_precipitation_uv_data(path_to_data=miocene_path)
mio_prec_diff = compute_anomalies(data_pi=pi_prec_data, data_exp=mio_prec_data, 
                             time="month", month="MAM")

mio_u_diff = compute_lterm_mean(data=mio_u_data, 
                             time="month", month="MAM")

mio_v_diff = compute_lterm_mean(data=mio_v_data, 
                             time="month", month="MAM")


# pliocene
plio_prec_data, plio_u_data, plio_v_data = read_precipitation_uv_data(path_to_data=pliocene_path)
plio_prec_diff = compute_anomalies(data_pi=pi_prec_data, data_exp=plio_prec_data,
                              time="month", month="MAM")

plio_u_diff = compute_lterm_mean(data=plio_u_data, 
                             time="month", month="MAM")

plio_v_diff = compute_lterm_mean(data=plio_v_data, 
                             time="month", month="MAM")

# holocene
holocene_prec_data, holocene_u_data, holocene_v_data = read_precipitation_uv_data(path_to_data=holocene_path)
holocene_prec_diff = compute_anomalies(data_pi=pi_prec_data, data_exp=holocene_prec_data,
                                  time="month", month="MAM")

holocene_u_diff = compute_lterm_mean(data=holocene_u_data, 
                             time="month", month="MAM")

holocene_v_diff = compute_lterm_mean(data=holocene_v_data, 
                             time="month", month="MAM")



# plotting
apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
            
projection = ccrs.Robinson(central_longitude=0, globe=None)

projection_trans = ccrs.PlateCarree()


fig,(ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(26,13), subplot_kw={"projection":projection})

plot_annual_mean(variable="Precipitation anomalies", data_alt=mio_prec_diff, ax=ax1,
                 cmap="BrBG", units="mm/month [MAM]", vmax=120, vmin=-120, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, plot_projection=projection, title="(a) Mid-Miocene - PI", 
                orientation="horizontal",  cbar_pos= [0.20, 0.05, 0.25, 0.02], sea_land_mask=mio_slm,
                domain="East Africa", data_v=mio_v_diff, data_u=mio_u_diff,
                show_arrow_scale=True, plot_winds=True)

# plot_annual_mean(variable="Precipitation", data_alt=mio_ea_high_diff, ax=ax2,
#                  cmap="BrBG", units="mm/month", vmax=200, vmin=-200, 
#                 levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=True,
#                 left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) Mid-Miocene (high EARS) - PI", 
#                 domain="West Africa", center=True, sea_land_mask=mio_slm,)

plot_annual_mean(variable="Precipitation anomalies", data_alt=plio_prec_diff, ax=ax2,
                 cmap="BrBG", units="mm/month", vmax=120, vmin=-120, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=True,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) Mid-Pliocene - PI", 
                domain="East Africa", data_v=plio_v_diff, data_u=plio_u_diff,
                show_arrow_scale=False, plot_winds=True)

plot_annual_mean(variable="Precipitation anomalies", data_alt=holocene_prec_diff, ax=ax3,
                 cmap="BrBG", units="mm/month", vmax=120, vmin=-120, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=True,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(c) Mid-Holocene - PI", 
                domain="East Africa", data_v=holocene_v_diff, data_u=holocene_u_diff,
                show_arrow_scale=False, plot_winds=True)


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "poster_map_eam.pdf"), format= "pdf", bbox_inches="tight", dpi=600)