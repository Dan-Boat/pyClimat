# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 15:00:08 2023

@author: dboateng
Analyse the latitudinal and longitudinal gradient of MAP, MAT, and isotopes accross the Alps. 
Compare the results with proxies (check for new compilations or better rely on current one)
"""

"""
Created on Tue Mar  7 14:59:50 2023

@author: dboateng
1. extract the spatial distribution of t2m, prec, d18op across the Alps (box plot)
Sections: low elevation, high elevation, Northern Alps, Eastern Alps, West-Central Alps, Mediterranean 
Also check thier seasonal differences (JJA and DJF, and annual)
"""
# import models
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


# read data

CTL_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

W2E1_PI_filename = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
W2E1_Mio278_filename = "a017_hpc-bw_e5w2.3_t159_MIO_W2E1_278ppm_t159l31.6h"
W2E1_Mio450_filename = "a016_hpc-bw_e5w2.3_t159_MIO_W2E1_450ppm_t159l31.6h"

W2E0_PI_filename="a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
W2E0_Mio278_filename="a019_hpc-bw_e5w2.3_t159_MIO_W2E0_278ppm_t159l31.6h"
W2E0_Mio450_filename="a018_hpc-bw_e5w2.3_t159_MIO_W2E0_450ppm_t159l31.6h"

W2E2_PI_filename="t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"
W2E2_Mio278_filename="a020_dkrz-levante_e5w2.3_t159_MIO_W2E2_278ppm_t159l31.6h"
W2E2_Mio450_filename="a021_dkrz-levante_e5w2.3_t159_MIO_W2E2_450ppm_t159l31.6h"

years = "1003_1017"

years_not_complete="1003_1010"


period = "1m"

CTL_data, CTL_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=CTL_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E1_PI_data, W2E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_278_data, W2E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_450_data, W2E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)

W2E0_PI_data, W2E0_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_PI_filename, years=years, 
                                          period=period, read_wiso=True)

W2E0_278_data, W2E0_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_Mio278_filename, 
                                                    years=years_not_complete, period=period, read_wiso=True)

W2E0_450_data, W2E0_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E0_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E2_PI_data, W2E2_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_PI_filename, years=years, 
                                          period=period, read_wiso=True)

W2E2_278_data, W2E2_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio278_filename, 
                                                    years=years, period=period, read_wiso=True)

W2E2_450_data, W2E2_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E2_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)

def extract_and_compute(W1E1, W2E0, W2E1, W2E2, varname, units, W1E1_wiso=None, W2E0_wiso=None, 
                        W2E1_wiso=None, W2E2_wiso=None):
    
    # extract variables 
    var_data_w1e1 = extract_var(Dataset=W1E1 , varname=varname, units=units, Dataset_wiso=W1E1_wiso)
    var_data_w2e0 = extract_var(Dataset=W2E0 , varname=varname, units=units, Dataset_wiso=W2E0_wiso)
    var_data_w2e1 = extract_var(Dataset=W2E1 , varname=varname, units=units, Dataset_wiso=W2E1_wiso)
    var_data_w2e2 = extract_var(Dataset=W2E2 , varname=varname, units=units, Dataset_wiso=W2E2_wiso)
    
    #compute the means
    w1e1_alt = compute_lterm_mean(data=var_data_w1e1, time="annual")
    w2e0_alt = compute_lterm_mean(data=var_data_w2e0, time="annual")
    w2e1_alt = compute_lterm_mean(data=var_data_w2e1, time="annual")
    w2e2_alt = compute_lterm_mean(data=var_data_w2e2, time="annual")
    
    return w1e1_alt, w2e0_alt, w2e1_alt, w2e2_alt


def extract_profile_topos(varname, units, dim, W1E1, W2E0, W2E1, W2E2, maxlon, minlon, maxlat, minlat, W1E1_wiso=None,
                          W2E0_wiso=None, W2E1_wiso=None, W2E2_wiso=None,):
    
    
    # extract var and compute
    w1e1_alt, w2e0_alt, w2e1_alt, w2e2_alt = extract_and_compute(W1E1, W2E0, W2E1, W2E2, varname, units,
                                                                 W1E1_wiso=W1E1_wiso, W2E0_wiso=W2E0_wiso, 
                                                                 W2E1_wiso=W2E1_wiso, W2E2_wiso=W2E2_wiso)
    
    
    profile_w1e1 = extract_profile(data = w1e1_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                   dim=dim, to_pandas=True)
    
    profile_w2e0 = extract_profile(data = w2e0_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                         dim=dim, to_pandas=True)
    
    profile_w2e1 = extract_profile(data = w2e1_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                         dim=dim, to_pandas=True)
    
    profile_w2e2 = extract_profile(data = w2e2_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                         dim=dim, to_pandas=True)
    
    
    profile_values = [profile_w1e1, profile_w2e0, profile_w2e1, profile_w2e2]
    
    return profile_values

def extract_profile_for_all_topo(varname, units, maxlon, minlon,
                                      maxlat, minlat, dim):
    
    topo_names = ["W1E1", "W2E0", "W2E1", "W2E2"]
    
    
    # extract variables 
    profile_values_pi = extract_profile_topos(varname=varname, units=units, 
                                                W1E1=CTL_data, W2E0=W2E0_PI_data, 
                                                W2E1=W2E1_PI_data, W2E2=W2E2_PI_data, maxlon=maxlon,
                                                minlon=minlon, maxlat=maxlat, minlat=minlat,
                                                W1E1_wiso=CTL_wiso, W2E0_wiso=W2E0_PI_wiso, 
                                                W2E1_wiso=W2E1_PI_wiso, W2E2_wiso=W2E2_PI_wiso,dim=dim)
    
    profile_values_mio278 = extract_profile_topos(varname=varname, units=units, 
                                                W1E1=W1E1_278_data, W2E0=W2E0_278_data, 
                                                W2E1=W2E1_278_data, W2E2=W2E2_278_data, maxlon=maxlon,
                                                minlon=minlon, maxlat=maxlat, minlat=minlat,
                                                W1E1_wiso=W1E1_278_wiso, W2E0_wiso=W2E0_278_wiso, 
                                                W2E1_wiso=W2E1_278_wiso, W2E2_wiso=W2E2_278_wiso,
                                                dim=dim)
    
    profile_values_mio450 = extract_profile_topos(varname=varname, units=units, 
                                                W1E1=W1E1_450_data, W2E0=W2E0_450_data, 
                                                W2E1=W2E1_450_data, W2E2=W2E2_450_data, maxlon=maxlon,
                                                minlon=minlon, maxlat=maxlat, minlat=minlat,
                                                W1E1_wiso=W1E1_450_wiso, W2E0_wiso=W2E0_450_wiso, 
                                                W2E1_wiso=W2E1_450_wiso, W2E2_wiso=W2E2_450_wiso,
                                                dim=dim)
    
    
    topo_names = ["W1E1", "W2E0", "W2E1", "W2E2"]
    
    for i,topo in enumerate(topo_names):
        if i ==0:
            
            df_pi = pd.DataFrame(index=profile_values_pi[i].index.values, columns=topo_names)
            df_mio278 = pd.DataFrame(index=profile_values_mio278[i].index.values, columns=topo_names)
            df_mio450 = pd.DataFrame(index=profile_values_mio450[i].index.values, columns=topo_names)
        
        df_pi[topo] = profile_values_pi[i]
        df_mio278[topo] = profile_values_mio278[i]
        df_mio450[topo] = profile_values_mio450[i]
        
    return df_pi, df_mio278, df_mio450


def plot_profiles_all(varname, units, data_pi, data_mio278, data_mio450, ax=None, path_to_store=None, filename=None,
                       colors=None, xlabel=True, ylabel=True, title=None, ax_legend=True,
                       ymin=None, ymax=None, dim="lon"):
    
    topo_names = ["W1E1", "W2E0", "W2E1", "W2E2"]
    
    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
    if colors is None:
        colors = ['#21130d', '#1f77b4', '#2ca02c', '#d62728']
    data_pi.plot(ax=ax, linestyle="-", color=colors, linewidth=3)    #marker="o", markersize=7
    data_mio278.plot(ax=ax, linestyle=":", color=colors, linewidth=3) 
    data_mio450.plot(ax=ax, linestyle="--", color=colors, linewidth=3)
            
            
    if ylabel:
        ax.set_ylabel(varname + " [" + units + "]", fontweight="bold", fontsize=20)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_yticklabels([])
    
    if xlabel is not None:
        if dim == "lon":
             ax.set_xlabel("Longitude [E°]", fontsize=22, fontweight="bold")
        elif dim == "lat":
             ax.set_xlabel("Latitude [N°]", fontsize=22, fontweight="bold")
        else:
            raise ValueError("Define dim as lat or lon")
        
        
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_xticklabels([])
        
    if all(parameter is not None for parameter in [ymax, ymin]):
        ax.set_ylim(ymin, ymax)
        
    if ax_legend:
        ax.legend(frameon=True, fontsize=22, loc="lower right")
    else:
        ax.legend([],[], frameon=False)
        
       
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="center")
        
    
        
        
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    
    if path_to_store is not None:
        plt.savefig(os.path.join(path_to_store, filename), bbox_inches="tight", format= "svg")
    
    
df_pi_d18op_lon, df_mio278_d18op_lon, df_mio450_18op_lon = extract_profile_for_all_topo(varname="d18op", units="per mil", 
                                            maxlon=25, minlon=-5, maxlat=47, minlat=45, dim="lon")

df_pi_prec_lon, df_mio278_prec_lon, df_mio450_prec_lon = extract_profile_for_all_topo(varname="prec", units="mm/month", 
                                            maxlon=25, minlon=-5, maxlat=47, minlat=45, dim="lon")

df_pi_temp_lon, df_mio278_temp_lon, df_mio450_temp_lon = extract_profile_for_all_topo(varname="temp2", units="°C", 
                                            maxlon=25, minlon=-5, maxlat=47, minlat=45, dim="lon")


df_pi_d18op_lat, df_mio278_d18op_lat, df_mio450_18op_lat = extract_profile_for_all_topo(varname="d18op", units="per mil", 
                                            maxlon=11, minlon=9, maxlat=54, minlat=40, dim="lat")

df_pi_prec_lat, df_mio278_prec_lat, df_mio450_prec_lat = extract_profile_for_all_topo(varname="prec", units="mm/month", 
                                            maxlon=11, minlon=9, maxlat=54, minlat=40, dim="lat")

df_pi_temp_lat, df_mio278_temp_lat, df_mio450_temp_lat = extract_profile_for_all_topo(varname="temp2", units="°C", 
                                            maxlon=11, minlon=9, maxlat=54, minlat=40, dim="lat")



apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
fig, ((ax1,ax2, ax3),( ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(36, 20), sharey=False)


plot_profiles_all(varname="$\delta^{18}$Op vs SMOW", units="‰", data_pi=df_pi_d18op_lon, data_mio278=df_mio278_d18op_lon, 
                  data_mio450=df_mio450_18op_lon, ax=ax1, ax_legend=True,ymax=-2, ymin=-20, dim="lon")

plot_profiles_all(varname="Temperature", units="°C", data_pi=df_pi_temp_lon, data_mio278=df_mio278_temp_lon, 
                  data_mio450=df_mio450_temp_lon, ax=ax2, ax_legend=False,ymax=25, ymin=-10, dim="lon")

plot_profiles_all(varname="Precipitation", units="mm/month", data_pi=df_pi_prec_lon, data_mio278=df_mio278_prec_lon, 
                  data_mio450=df_mio450_prec_lon, ax=ax3, ax_legend=False,ymax=250, ymin=20, dim="lon")


plot_profiles_all(varname="$\delta^{18}$Op vs SMOW", units="‰", data_pi=df_pi_d18op_lat, data_mio278=df_mio278_d18op_lat, 
                  data_mio450=df_mio450_18op_lat, ax=ax4, ax_legend=False,ymax=-2, ymin=-20, dim="lat")

plot_profiles_all(varname="Temperature", units="°C", data_pi=df_pi_temp_lat, data_mio278=df_mio278_temp_lat, 
                  data_mio450=df_mio450_temp_lat, ax=ax5, ax_legend=False,ymax=25, ymin=-10, dim="lat")

plot_profiles_all(varname="Precipitation", units="mm/month", data_pi=df_pi_prec_lat, data_mio278=df_mio278_prec_lat, 
                  data_mio450=df_mio450_prec_lat, ax=ax6, ax_legend=False,ymax=250, ymin=20, dim="lat")


plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
plt.savefig(os.path.join(path_to_plots, "profiles_d18op_temp_prec.svg"), format= "svg", bbox_inches="tight", dpi=600)