# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:54:17 2023

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xr
import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import scatter_plot_laspe_rate
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, extract_transect, extract_profile


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

proxies_path_low = "C:/Users/dboateng/Desktop/Datasets/Alps_d18op/low_elevation.csv"
proxies_path_high = "C:/Users/dboateng/Desktop/Datasets/Alps_d18op/high_elevation.csv"


proxy_data_low = pd.read_csv(proxies_path_low)
proxy_data_high = pd.read_csv(proxies_path_high)

# read data

W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"


W2E1_Mio278_filename = "a017_hpc-bw_e5w2.3_t159_MIO_W2E1_278ppm_t159l31.6h"
W2E1_Mio450_filename = "a016_hpc-bw_e5w2.3_t159_MIO_W2E1_450ppm_t159l31.6h"


W2E0_Mio278_filename="a019_hpc-bw_e5w2.3_t159_MIO_W2E0_278ppm_t159l31.6h"
W2E0_Mio450_filename="a018_hpc-bw_e5w2.3_t159_MIO_W2E0_450ppm_t159l31.6h"


W2E2_Mio278_filename="a020_dkrz-levante_e5w2.3_t159_MIO_W2E2_278ppm_t159l31.6h"
W2E2_Mio450_filename="a021_dkrz-levante_e5w2.3_t159_MIO_W2E2_450ppm_t159l31.6h"


W2E15_Mio278_filename = "a023_dkrz-levante_e5w2.3_t159_MIO_W2E1.5_278ppm_t159l31.6h"
W2E15_Mio450_filename = "a022_hpc-bw_e5w2.3_t159_MIO_W2E1.5_450ppm_t159l31.6h"


W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E1_278_data, W2E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_450_data, W2E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)

W2E15_278_data, W2E15_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E15_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E15_450_data, W2E15_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E15_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E0_278_data, W2E0_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_Mio278_filename, 
                                                    years=years, period=period, read_wiso=True)

W2E0_450_data, W2E0_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E2_278_data, W2E2_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio278_filename, 
                                                    years=years, period=period, read_wiso=True)

W2E2_450_data, W2E2_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_PI_data, W1E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)

def extract_vars_and_analysis(data, wiso):

    # d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    d18op_carb = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso,
                                  use_PDB=True)
    
    d18op_carb_alt = compute_lterm_mean(data=d18op_carb, time="annual")
    
    lon_d18op_carb = extract_profile(d18op_carb_alt, maxlon=25, minlon=-5, maxlat=47, minlat=45, dim="lon",
                                     to_pandas=True)
    lat_d18op_carb = extract_profile(d18op_carb_alt, maxlon=11, minlon=9, maxlat=54, minlat=40, dim="lat", 
                                     to_pandas=True)
    
    return_data = {"lon": lon_d18op_carb, "lat": lat_d18op_carb}
    
    return return_data




def plot_profiles_all(varname, units, data_mio278, data_mio450, ax=None, path_to_store=None, filename=None,
                       colors=None, xlabel=True, ylabel=True, title=None, ax_legend=True,
                       ymin=None, ymax=None, dim="lon", proxy=None, control_data=None,
                       marker="v", proxy_label="low elevation samples", vmax=25, vmin=5, fig=None):
    
    topo_names = ["W1E1", "W2E0", "W2E1", "W2E1.5", "W2E2"]
    
    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
    if colors is None:
        colors = ['#21130d', '#1f77b4', '#2ca02c', '#d62728', purple]
        
    
    data_mio278.plot(ax=ax, linestyle=":", color=colors, linewidth=3) 
    data_mio450.plot(ax=ax, linestyle="--", color=colors, linewidth=3) #marker="o", markersize=7
    
    if control_data is not None:
        control_data.plot(ax=ax, linestyle="-", color=black, linewidth=3)
    
    if proxy is not None:
        p = ax.scatter(x=proxy[dim], y=proxy["d18op"], s=500, c=proxy["age"], cmap="winter",
                   marker=marker, vmax=vmax, vmin=vmin, alpha=0.95)
        
        ax.errorbar(x=proxy[dim], y=proxy["d18op"], yerr=proxy["d18op_error"], xerr=None,
                    ls="none", color="black", capthick=3, capsize=4, label=proxy_label)
        
        cbar_pos = [0.94, 0.30, 0.03, 0.45]
        
        cbar_ax = fig.add_axes(cbar_pos)
        
        plt.colorbar(p, label="age (Ma)", shrink=0.3, cax=cbar_ax)
            
    if ylabel:
        ax.set_ylabel(varname + " [" + units + "]", fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_yticklabels([])
    
    if xlabel is not None:
        if dim == "lon":
             ax.set_xlabel("Longitude [E°]", fontsize=28, fontweight="bold")
        elif dim == "lat":
             ax.set_xlabel("Latitude [N°]", fontsize=28, fontweight="bold")
        else:
            raise ValueError("Define dim as lat or lon")
        
        
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_xticklabels([])
        
    if all(parameter is not None for parameter in [ymax, ymin]):
        ax.set_ylim(ymin, ymax)
        
    if ax_legend:
        ax.legend(frameon=True, fontsize=24, bbox_to_anchor=(0.01, 1.05, 1, 0.102,), loc=3, borderaxespad=0,
        ncol=4)
    else:
        ax.legend([],[], frameon=False)
        
       
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="center")
        
    
        
        
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    
    if path_to_store is not None:
        plt.savefig(os.path.join(path_to_store, filename), bbox_inches="tight", format= "svg")
# extract all

extracts_278 = {}
extracts_450 = {}

labels = ["W1E1", "W2E0", "W2E1", "W2E1.5", "W2E2"]  
exp_data_278 = [W1E1_278_data, W2E0_278_data, W2E1_278_data, W2E15_278_data, W2E2_278_data]
exp_wiso_278 = [W1E1_278_wiso, W2E0_278_wiso, W2E1_278_wiso, W2E15_278_wiso, W2E2_278_wiso] 

exp_data_450 = [W1E1_450_data, W2E0_450_data, W2E1_450_data, W2E15_450_data, W2E2_450_data]
exp_wiso_450 = [W1E1_450_wiso, W2E0_450_wiso, W2E1_450_wiso, W2E15_450_wiso, W2E2_450_wiso] 

for i,topo in enumerate(labels):
    extracts_278[topo] = extract_vars_and_analysis(data=exp_data_278[i], wiso=exp_wiso_278[i])
    extracts_450[topo] = extract_vars_and_analysis(data=exp_data_450[i], wiso=exp_wiso_450[i])
    
    if i == 0:
        df_mio278_lon = pd.DataFrame(index=extracts_278.get(topo)["lon"].index.values, columns=labels)
        df_mio450_lon = pd.DataFrame(index=extracts_450.get(topo)["lon"].index.values, columns=labels)
        
        df_mio278_lat = pd.DataFrame(index=extracts_278.get(topo)["lat"].index.values, columns=labels)
        df_mio450_lat = pd.DataFrame(index=extracts_450.get(topo)["lat"].index.values, columns=labels)
        
    df_mio278_lon[topo] = extracts_278.get(topo)["lon"]
    df_mio278_lat[topo] = extracts_278.get(topo)["lat"]
    df_mio450_lon[topo] = extracts_450.get(topo)["lon"]
    df_mio450_lat[topo] = extracts_450.get(topo)["lat"]
    
    
    
PI_profile = extract_vars_and_analysis(data=W1E1_PI_data, wiso=W1E1_PI_wiso)
df_PI_lon = pd.DataFrame(index=PI_profile.get("lon").index.values, columns=["PI(W1E1)"]) 
df_PI_lon["PI(W1E1)"] = PI_profile.get("lon") 

df_PI_lat = pd.DataFrame(index=PI_profile.get("lat").index.values, columns=["PI(W1E1)"]) 
df_PI_lat["PI(W1E1)"] = PI_profile.get("lat")



def plot_profile_low():
    apply_style(fontsize=28, style="seaborn-talk", linewidth=3,)
    
    fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(30, 15), sharey=False)
    
    
    plot_profiles_all(varname="$\delta^{18}$Op vs SMOW", units="‰", data_mio278=df_mio278_lon, 
                      data_mio450=df_mio450_lon, ax=ax1, ax_legend=True,ymax=-2, ymin=-20, dim="lon", 
                      control_data=df_PI_lon, proxy=proxy_data_low, fig=fig, 
                      proxy_label="Low elevation proxies", marker="v",  vmax=20, vmin=10)
    
    plot_profiles_all(varname="$\delta^{18}$Op vs SMOW", units="‰", data_mio278=df_mio278_lat, 
                      data_mio450=df_mio450_lat, ax=ax2, ax_legend=False,ymax=-2, ymin=-20, dim="lat", 
                      control_data=df_PI_lat, proxy=proxy_data_low, ylabel=True, fig=fig, 
                      proxy_label="Low elevation proxies",  marker="v", vmax=20, vmin=10)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.90, top=0.97, bottom=0.05, wspace=0.12)
    plt.savefig(os.path.join(path_to_plots, "d18Op_profile_with_low.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
    
def plot_profile_high():
    apply_style(fontsize=28, style="seaborn-talk", linewidth=3,)
    
    fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(30, 15), sharey=False)
    
    
    plot_profiles_all(varname="$\delta^{18}$Op vs SMOW", units="‰", data_mio278=df_mio278_lon, 
                      data_mio450=df_mio450_lon, ax=ax1, ax_legend=True,ymax=-2, ymin=-20, dim="lon", 
                      control_data=df_PI_lon, proxy=proxy_data_high, fig=fig, 
                      proxy_label="High elevation proxies", marker="o",  vmax=20, vmin=10)
    
    plot_profiles_all(varname="$\delta^{18}$Op vs SMOW", units="‰", data_mio278=df_mio278_lat, 
                      data_mio450=df_mio450_lat, ax=ax2, ax_legend=False,ymax=-2, ymin=-20, dim="lat", 
                      control_data=df_PI_lat, proxy=proxy_data_high, ylabel=True, fig=fig, 
                      proxy_label="High elevation proxies",  marker="o", vmax=20, vmin=10)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.90, top=0.97, bottom=0.05, wspace=0.12)
    plt.savefig(os.path.join(path_to_plots, "d18Op_profile_with_high.svg"), format= "svg", bbox_inches="tight", dpi=600)


plot_profile_low()
plot_profile_high()