# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 15:15:00 2024

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
#lgm_path= "D:/Datasets/Model_output_pst/LGM/MONTHLY_MEANS/"
holocene_path = "D:/Datasets/Model_output_pst/MH/MONTHLY_MEANS/"


mpi_pi_path = "D:/Datasets/CMIP6/PMIP/postprocessed/PI/MPI-ESM1-2-LR/"
mpi_cmip_85_path = "D:/Datasets/CMIP6/CMIP/SSP585/MPI-ESMI-2-LR/monthly/"
mpi_cmip_26_path = "D:/Datasets/CMIP6/CMIP/SSP126/MPI-ESMI-2-LR/monthly/"


#proxies
path_to_bartlein_mh = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/map_delta_06ka_ALL_grid_2x2_global.csv"
#path_to_bartlein_lgm = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/map_delta_21ka_ALL_grid_2x2.csv"
path_to_tierney = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/tierney_data.csv"
path_to_plio = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/41467_2022_28814_MOESM4_ESM.csv"


# read data
echam_filename = "1003_1017_1m_mlterm.nc"
mpi_esm_filename = "tp_monthly.nc"


mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


from2070to2100 = pd.date_range(start='2070-01-01', end='2100-12-31', freq='MS')




#proxies extracts

df_b_mh = pd.read_csv(path_to_bartlein_mh)

mh_b_wetter = df_b_mh[df_b_mh["map_anm_mean"] > 0]
mh_b_drier = df_b_mh[df_b_mh["map_anm_mean"] < 0]


#df_b_lgm = pd.read_csv(path_to_bartlein_lgm)

# lgm_b_wetter = df_b_lgm[df_b_lgm["map_anm_mean"] > 0]
# lgm_b_drier = df_b_lgm[df_b_lgm["map_anm_mean"] < 0]

df_tierney = pd.read_csv(path_to_tierney)

df_plio_data = pd.read_csv(path_to_plio)

plio_wetter = df_plio_data[df_plio_data["Interpretation"] == "Wetter"]
plio_drier = df_plio_data[df_plio_data["Interpretation"] == "Drier"]
plio_no_change = df_plio_data[df_plio_data["Interpretation"] == "No Change"]


def read_precipitation_data(path_to_data, mpi_varname=None, echam=True, time_range=None):
    
    if echam:
        data = read_from_path(path=path_to_data, filename=echam_filename)
        prec = extract_var(Dataset=data , varname="prec", units="mm/month")
    else:
        prec = read_from_path(path=path_to_data, filename=mpi_esm_filename,
                              varname=mpi_varname, decode=True, )*60*60*24*365  #mm/month
        
        if time_range is not None:
            prec = prec.sel(time=time_range)
        
    return prec


def compute_anomalies(data_pi, data_exp, time="annual"):
    
    prec_diff_alt = compute_lterm_diff(data_control=data_pi, data_main=data_exp)
    
    return prec_diff_alt


#pi 
pi_echam_data = read_precipitation_data(path_to_data=pi_path)
pi_mpi_esm_data = read_precipitation_data(path_to_data=mpi_pi_path, mpi_varname="pr",
                                          echam=False)

pi_echam_alt = compute_lterm_mean(data=pi_echam_data, time="annual")
# miocene
mio_data = read_precipitation_data(path_to_data=miocene_path)
mio_diff = compute_anomalies(data_pi=pi_echam_data, data_exp=mio_data)

# pliocene
plio_data = read_precipitation_data(path_to_data=pliocene_path)
plio_diff = compute_anomalies(data_pi=pi_echam_data, data_exp=plio_data)

# holocene
holocene_data = read_precipitation_data(path_to_data=holocene_path)
holocene_diff = compute_anomalies(data_pi=pi_echam_data, data_exp=holocene_data)


#future
cmip_26_data = read_precipitation_data(path_to_data=mpi_cmip_26_path, mpi_varname="tp",
                                       echam=False, time_range=from2070to2100)
cmip_26_diff = compute_anomalies(data_pi=pi_mpi_esm_data, data_exp=cmip_26_data)

cmip_85_data = read_precipitation_data(path_to_data=mpi_cmip_85_path, mpi_varname="tp",
                                       echam=False, time_range=from2070to2100)
cmip_85_diff = compute_anomalies(data_pi=pi_mpi_esm_data, data_exp=cmip_85_data)



# plotting
apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
            
projection = ccrs.Robinson(central_longitude=0, globe=None)

projection_trans = ccrs.PlateCarree()


fig,(ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(nrows=1, ncols=6, figsize=(32,13), subplot_kw={"projection":projection})

plot_annual_mean(variable="Precipitation anomalies", data_alt=mio_diff, ax=ax1,
                 cmap="BrBG", units="mm/month", vmax=120, vmin=-120, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(a) Mid-Miocene - PI", 
                orientation="horizontal",  cbar_pos= [0.20, 0.05, 0.25, 0.02], sea_land_mask=mio_slm,
                domain="Africa")

plot_annual_mean(variable="Precipitation anomalies", data_alt=plio_diff, ax=ax2,
                 cmap="BrBG", units="mm/month", vmax=120, vmin=-120, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) Mid-Pliocene - PI", 
                domain="Africa")

plot_annual_mean(variable="Precipitation anomalies", data_alt=holocene_diff, ax=ax3,
                 cmap="BrBG", units="mm/month", vmax=120, vmin=-120, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(c) Mid-Holocene - PI", 
                domain="Africa")

plot_annual_mean(variable="Precipitation", data_alt=pi_echam_alt, ax=ax4,
                 cmap="YlGnBu", units="mm/month", vmax=300, vmin=50, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(d) Pre-Industrial", 
                domain="Africa", cbar_pos= [0.65, 0.05, 0.25, 0.02], center=False)

plot_annual_mean(variable="Precipitation anomalies", data_alt=cmip_26_diff, ax=ax5,
                 cmap="BrBG", units="mm/month", vmax=120, vmin=-120, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(e) CMIP6 SSP26 [2070-2100] - PI", 
                domain="Africa")

plot_annual_mean(variable="Precipitation anomalies", data_alt=cmip_85_diff, ax=ax6,
                 cmap="BrBG", units="mm/month", vmax=120, vmin=-120, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(f) CMIP6 SSP85 [2070-2100] - PI", 
                domain="Africa")




#plot mid-Holocene Bartlien data
sc = ax3.scatter(mh_b_wetter['lon'], mh_b_wetter['lat'], s=120, color="#004d00",
                edgecolor="k", transform = projection_trans)

sc = ax3.scatter(mh_b_drier['lon'], mh_b_drier['lat'], s=120, color="#663300",
                edgecolor="k", transform = projection_trans)

sc = ax3.scatter(df_tierney['lon'], df_tierney['lat'], s=120, color="#004d00",
                edgecolor="red", transform = projection_trans)


# plot lgm Bartlien data 
# sc = ax2.scatter(lgm_b_wetter['lon'], lgm_b_wetter['lat'], s=120, color="#004d00",
#                 edgecolor="k", transform = projection_trans)

# sc = ax2.scatter(lgm_b_drier['lon'], lgm_b_drier['lat'], s=120, color="#663300",
#                 edgecolor="k", transform = projection_trans)


# plot plio from Rang
sc = ax2.scatter(plio_wetter['Lon'], plio_wetter['Lat '], s=120, color="#004d00",
                edgecolor="k", transform = projection_trans)

sc = ax2.scatter(plio_drier['Lon'], plio_drier['Lat '], s=120, color="#663300",
                edgecolor="k", transform = projection_trans)
    
    
    


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "poster_fig1.pdf"), format= "pdf", bbox_inches="tight", dpi=600)





