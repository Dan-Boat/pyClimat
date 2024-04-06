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

path_to_bartlein_mh = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/map_delta_06ka_ALL_grid_2x2_global.csv"
path_to_bartlein_lgm = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/map_delta_21ka_ALL_grid_2x2.csv"
path_to_tierney = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/tierney_data.csv"
path_to_plio = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/41467_2022_28814_MOESM4_ESM.csv"


df_b_mh = pd.read_csv(path_to_bartlein_mh)

mh_b_wetter = df_b_mh[df_b_mh["map_anm_mean"] > 0]
mh_b_drier = df_b_mh[df_b_mh["map_anm_mean"] < 0]


df_b_lgm = pd.read_csv(path_to_bartlein_lgm)

lgm_b_wetter = df_b_lgm[df_b_lgm["map_anm_mean"] > 0]
lgm_b_drier = df_b_lgm[df_b_lgm["map_anm_mean"] < 0]

df_tierney = pd.read_csv(path_to_tierney)

df_plio_data = pd.read_csv(path_to_plio)

plio_wetter = df_plio_data[df_plio_data["Interpretation"] == "Wetter"]
plio_drier = df_plio_data[df_plio_data["Interpretation"] == "Drier"]
plio_no_change = df_plio_data[df_plio_data["Interpretation"] == "No Change"]


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





def plot_global():
    apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    projection_trans = ccrs.PlateCarree()
    
    
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
    
    
    
    #plot mid-Holocene Bartlien data
    sc = ax1.scatter(mh_b_wetter['lon'], mh_b_wetter['lat'], s=120, color="#004d00",
                    edgecolor="k", transform = projection_trans)
    
    sc = ax1.scatter(mh_b_drier['lon'], mh_b_drier['lat'], s=120, color="#663300",
                    edgecolor="k", transform = projection_trans)
    
    sc = ax1.scatter(df_tierney['lon'], df_tierney['lat'], s=120, color="#004d00",
                    edgecolor="k", transform = projection_trans, marker="s")
    
    
    # plot lgm Bartlien data 
    sc = ax2.scatter(lgm_b_wetter['lon'], lgm_b_wetter['lat'], s=120, color="#004d00",
                    edgecolor="k", transform = projection_trans)
    
    sc = ax2.scatter(lgm_b_drier['lon'], lgm_b_drier['lat'], s=120, color="#663300",
                    edgecolor="k", transform = projection_trans)
    
    
    # plot plio from Rang
    sc = ax3.scatter(plio_wetter['Lon'], plio_wetter['Lat '], s=120, color="#004d00",
                    edgecolor="k", transform = projection_trans)
    
    sc = ax3.scatter(plio_drier['Lon'], plio_drier['Lat '], s=120, color="#663300",
                    edgecolor="k", transform = projection_trans)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, "global_anomalies.png"), format= "png", bbox_inches="tight", dpi=600)
    
    
    
apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)

projection_trans = ccrs.PlateCarree()


fig,((ax1, ax2, ax3)) = plt.subplots(nrows=1, ncols=3, figsize=(25,15), subplot_kw={"projection":projection})

plot_annual_mean(variable="Precipitation anomalies", data_alt=mh_diff.get("precipitation"), ax=ax1,
                 cmap="BrBG", units="mm/month", vmax=100, vmin=-100, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(a) MH - PI", 
                orientation="horizontal",  cbar_pos= [0.35, 0.25, 0.25, 0.02], domain="Africa")

plot_annual_mean(variable="Precipitation anomalies", data_alt=lgm_diff.get("precipitation"), ax=ax2,
                 cmap="BrBG", units="mm/month", vmax=100, vmin=-100, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) LGM - PI", 
                domain="Africa")

plot_annual_mean(variable="Precipitation anomalies", data_alt=mplio_diff.get("precipitation"), ax=ax3,
                 cmap="BrBG", units="mm/month", vmax=100, vmin=-100, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(c) mPLIO - PI", 
                domain="Africa")


#plot mid-Holocene Bartlien data
sc = ax1.scatter(mh_b_wetter['lon'], mh_b_wetter['lat'], s=250, color="#004d00",
                edgecolor="k", transform = projection_trans, linewidth=3)

sc = ax1.scatter(mh_b_drier['lon'], mh_b_drier['lat'], s=250, color="#663300",
                edgecolor="k", transform = projection_trans, linewidth=3)

sc = ax1.scatter(df_tierney['lon'], df_tierney['lat'], s=250, color="#004d00",
                edgecolor="red", transform = projection_trans, linewidth=3)


# plot lgm Bartlien data 
sc = ax2.scatter(lgm_b_wetter['lon'], lgm_b_wetter['lat'], s=250, color="#004d00",
                edgecolor="k", transform = projection_trans, linewidth=3)

sc = ax2.scatter(lgm_b_drier['lon'], lgm_b_drier['lat'], s=250, color="#663300",
                edgecolor="k", transform = projection_trans, linewidth=3)


# plot plio from Rang
sc = ax3.scatter(plio_wetter['Lon'], plio_wetter['Lat '], s=250, color="#004d00",
                edgecolor="k", transform = projection_trans, linewidth=3)

sc = ax3.scatter(plio_drier['Lon'], plio_drier['Lat '], s=250, color="#663300",
                edgecolor="k", transform = projection_trans, linewidth=3)


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "africa_anomalies_with_proxies.pdf"), format= "pdf", bbox_inches="tight", dpi=600)