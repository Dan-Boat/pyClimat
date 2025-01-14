# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 13:11:23 2024

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/armelle"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

path_miomip_temp = "D:/Datasets/zenodo/MioMIP/MioMIP1.nc"

cosmos_278 = "COSMOS Mid Miocene 278ppm"
cosmos_450 = "COSMOS Mid Miocene 450ppm"

data_miomip = xr.open_dataset(path_miomip_temp)
cosmos_278_data = data_miomip.TS_ANO.sel(exp=cosmos_278)
cosmos_450_data = data_miomip.TS_ANO.sel(exp=cosmos_450)


mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"
W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

# reading data 
# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_PI_data, W1E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)


# extract values and compute long-term means

def extract_vars_and_analysis(data, wiso, pi_data, pi_wiso):
    
    temp2 = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    
    temp2_pi = extract_var(Dataset=pi_data , varname="temp2", units="°C")
    prec_pi = extract_var(Dataset= pi_data , varname="prec", units="mm/month")
    d18op_pi = extract_var(Dataset=pi_data , varname="d18op", units="per mil", Dataset_wiso= pi_wiso)
    
    #compute climatologies difference
    
    temp2_diff = compute_lterm_diff(data_control=temp2_pi, data_main=temp2, time="annual")
    prec_diff = compute_lterm_diff(data_control=prec_pi, data_main=prec, time="annual")
    d18op_diff = compute_lterm_diff(data_control=d18op_pi, data_main=d18op, time="annual")
    
    
    return_data = {"temperature":temp2_diff, "precipitation":prec_diff, "d18Op":d18op_diff,}
    
    return return_data

# compute area-weighted means

def compute_spatial_mean(data):
    weights = np.cos(np.deg2rad(data.lat)) #np.sqrt(np.abs(np.cos(data_raw.lat*np.pi/180)))

    weights.name = "weights"

    means = data.weighted(weights).mean(dim=("lon", "lat"), skipna=True)
    
    return means


    
# read data 
Mio278_data = extract_vars_and_analysis(data=W1E1_278_data, wiso=W1E1_278_wiso, pi_data=W1E1_PI_data, 
                                        pi_wiso=W1E1_PI_wiso)
Mio450_data = extract_vars_and_analysis(data=W1E1_450_data, wiso=W1E1_450_wiso, pi_data=W1E1_PI_data, 
                                        pi_wiso=W1E1_PI_wiso)



apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)

fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(18,16), subplot_kw={"projection":projection})

mean_echam_278 = compute_spatial_mean(Mio278_data.get("temperature"))
plot_annual_mean(variable="Temperature anomalies", data_alt=Mio278_data.get("temperature"), ax=ax1,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(a) ECHAM5-wiso MIO 278ppm - PI ({:.2f})".format(mean_echam_278), orientation="horizontal",
                cbar_pos= [0.35, 0.05, 0.35, 0.02])

mean_cosmos_278 = compute_spatial_mean(cosmos_278_data)
plot_annual_mean(variable="Temperature anomalies", data_alt=cosmos_278_data, ax=ax2,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(b) COSMOS MIO 450ppm - PI ({:.2f})".format(mean_cosmos_278))

mean_echam_450 = compute_spatial_mean(Mio450_data.get("temperature"))
plot_annual_mean(variable="Temperature anomalies", data_alt=Mio450_data.get("temperature"), ax=ax3,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(c) ECHAM5-wiso MIO 450ppm - PI ({:.2f})".format(mean_echam_450))

mean_cosmos_450 = compute_spatial_mean(cosmos_450_data)
plot_annual_mean(variable="Temperature anomalies", data_alt=cosmos_450_data, ax=ax4,
                 cmap="RdBu_r", units="°C", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, title="(d) COSMOS MIO 450ppm - PI ({:.2f})".format(mean_cosmos_450))


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05,)
plt.savefig(os.path.join(path_to_plots, "global_anomalies_cosmos.pdf"), format= "pdf", bbox_inches="tight", dpi=600)
plt.savefig(os.path.join(path_to_plots, "global_anomalies_cosmos.png"), format= "png", bbox_inches="tight", dpi=600)
plt.show()


