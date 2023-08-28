# -*- coding: utf-8 -*-
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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import scatter_plot_laspe_rate
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, extract_transect, linregression
from pyClimat.variables import extract_var


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


def extract_region_topos(varname, units, W1E1, W2E0, W2E1, W2E2, maxlon, minlon, maxlat, minlat, minelev=None, 
                         maxelev=None, land_mask=False, W1E1_wiso=None, W2E0_wiso=None, 
                                                 W2E1_wiso=None, W2E2_wiso=None, cal_mean=False,
                                                 lapse_rate= -0.002):
    
    
    # extract var and compute
    w1e1_alt, w2e0_alt, w2e1_alt, w2e2_alt = extract_and_compute(W1E1, W2E0, W2E1, W2E2, varname, units,
                                                                 W1E1_wiso=W1E1_wiso, W2E0_wiso=W2E0_wiso, 
                                                                 W2E1_wiso=W2E1_wiso, W2E2_wiso=W2E2_wiso)
    
    
    region_w1e1 = extract_transect(data = w1e1_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                         sea_land_mask=land_mask, minelev=minelev, maxelev=maxelev, Dataset=W1E1)
    
    region_w2e0 = extract_transect(data = w2e0_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                         sea_land_mask=land_mask, minelev=minelev, maxelev=maxelev, Dataset=W2E0)
    
    region_w2e1 = extract_transect(data = w2e1_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                         sea_land_mask=land_mask, minelev=minelev, maxelev=maxelev, Dataset=W2E1)
    
    region_w2e2 = extract_transect(data = w2e2_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                         sea_land_mask=land_mask, minelev=minelev, maxelev=maxelev, Dataset=W2E2)
    
    if cal_mean == True:
        region_values = [region_w1e1.mean().values.ravel() * lapse_rate, region_w2e0.mean().values.ravel()*lapse_rate, 
                         region_w2e1.mean().values.ravel()*lapse_rate, region_w2e2.mean().values.ravel()*lapse_rate]
       
    else:
        region_values = [region_w1e1.values.ravel(), region_w2e0.values.ravel(), 
                         region_w2e1.values.ravel(), region_w2e2.values.ravel()]
    
    return region_values



def extract_var_for_section_all_topo(varname, units, maxlon, minlon,
                                      maxlat, minlat, minelev=None, maxelev=None, land_mask=False,
                                      cal_mean=False,):
    
    topo_names = ["W1E1", "W2E0", "W2E1", "W2E2"]
    
    
    # extract variables 
    region_values_pi = extract_region_topos(varname=varname, units=units, 
                                                W1E1=CTL_data, W2E0=W2E0_PI_data, 
                                                W2E1=W2E1_PI_data, W2E2=W2E2_PI_data, maxlon=maxlon,
                                                minlon=minlon, maxlat=maxlat, minlat=minlat,
                                                minelev=minelev, maxelev=maxelev, land_mask=land_mask,
                                                W1E1_wiso=CTL_wiso, W2E0_wiso=W2E0_PI_wiso, 
                                                W2E1_wiso=W2E1_PI_wiso, W2E2_wiso=W2E2_PI_wiso, 
                                                cal_mean=cal_mean)
    
    region_values_mio278 = extract_region_topos(varname=varname, units=units, 
                                                W1E1=W1E1_278_data, W2E0=W2E0_278_data, 
                                                W2E1=W2E1_278_data, W2E2=W2E2_278_data, maxlon=maxlon,
                                                minlon=minlon, maxlat=maxlat, minlat=minlat,
                                                minelev=minelev, maxelev=maxelev, land_mask=land_mask,
                                                W1E1_wiso=W1E1_278_wiso, W2E0_wiso=W2E0_278_wiso, 
                                                W2E1_wiso=W2E1_278_wiso, W2E2_wiso=W2E2_278_wiso,
                                                cal_mean=cal_mean)
    
    region_values_mio450 = extract_region_topos(varname=varname, units=units, 
                                                W1E1=W1E1_450_data, W2E0=W2E0_450_data, 
                                                W2E1=W2E1_450_data, W2E2=W2E2_450_data, maxlon=maxlon,
                                                minlon=minlon, maxlat=maxlat, minlat=minlat,
                                                minelev=minelev, maxelev=maxelev, land_mask=land_mask,
                                                W1E1_wiso=W1E1_450_wiso, W2E0_wiso=W2E0_450_wiso, 
                                                W2E1_wiso=W2E1_450_wiso, W2E2_wiso=W2E2_450_wiso,
                                                cal_mean=cal_mean)
    
    
    for i,topo in enumerate(topo_names):
        if i ==0:
            len_pi = len(region_values_pi[i])
            len_mio278 = len(region_values_mio278[i])
            len_mio450 = len(region_values_mio450[i])
            
            df_pi = pd.DataFrame(index=np.arange(len_pi), columns=topo_names).assign(Paleoclimate="PI")
            df_mio278 = pd.DataFrame(index=np.arange(len_mio278), columns=topo_names).assign(Paleoclimate="MIO 278ppm")
            df_mio450 = pd.DataFrame(index=np.arange(len_mio450), columns=topo_names).assign(Paleoclimate="MIO 450ppm")
        
        df_pi[topo] = region_values_pi[i]
        df_mio278[topo] = region_values_mio278[i]
        df_mio450[topo] = region_values_mio450[i]
    
   
    cdf = pd.concat([df_pi, df_mio278, df_mio450])
    mdf = pd.melt(cdf, id_vars=["Paleoclimate"],var_name=varname)
    
    return mdf

def plot_cross_section(varname, units, data, hue, ax=None, path_to_store=None, filename=None,
                       colors=None, xlabel=True, ylabel=True, title=None, ax_legend=True,
                       ymin=None, ymax=None, points_data=None, point_hue=None):
    
    topo_names = ["W1E1", "W2E0", "W2E1", "W2E2"]
    apply_style(fontsize=22, style=None, linewidth=2)
    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
    boxplot = sns.violinplot(data=data, x="Paleoclimate", y="value", saturation=0.8, ax=ax,
                          hue=hue, scale="width") # count, width, area
    
    
    if points_data is not None:
        
        pointplot = sns.pointplot(data=points_data, x="Paleoclimate", y="value", ax=ax,
                              hue=point_hue, markers="x", linestyles="", scale=2,
                              dodge=True)
    if colors is not None:
        
        for patch, color in zip(boxplot["boxes"], colors):
            patch.set_facecolor(color)
            
            
    if ylabel:
        ax.set_ylabel(varname + " [" + units + "]", fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_yticklabels([])
    
    if xlabel is not None:
        ax.set_xlabel("Experiment", fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_xticklabels([])
        
    if all(parameter is not None for parameter in [ymax, ymin]):
        ax.set_ylim(ymin, ymax)
        
    if ax_legend:
        ax.legend(frameon=True, fontsize=24,
                  bbox_to_anchor=(0.01, 1.05, 1, 0.102,), loc=3, borderaxespad=0,
                  ncol=4)
    else:
        ax.legend([],[], frameon=False)
        
       
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 28, "fontweight":"bold"}, loc="center")
        
    
        
        
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    
    if path_to_store is not None:
        plt.savefig(os.path.join(path_to_store, filename), bbox_inches="tight", format= "svg")


# extract
def estimate_all_transect(varname, units, cal_mean=False):        
    mdf_low = extract_var_for_section_all_topo(varname=varname, units=units, maxlon=19, minlon=-2, maxlat=51, minlat=42,
                                     maxelev=500, land_mask=True, cal_mean=cal_mean)
    mdf_high = extract_var_for_section_all_topo(varname=varname, units=units, maxlon=19, minlon=-2, maxlat=51, minlat=42,
                                     minelev=1000, land_mask=True, cal_mean=cal_mean)
    
    
    # mdf_west = extract_var_for_section_all_topo(varname=varname, units=units, maxlon=8, minlon=1, maxlat=47, minlat=44,
    #                                  land_mask=True)
    # mdf_north = extract_var_for_section_all_topo(varname=varname, units=units, maxlon=16, minlon=5, maxlat=51, minlat=46.5,
    #                                  land_mask=True)
    # mdf_south = extract_var_for_section_all_topo(varname=varname, units=units, maxlon=16, minlon=7.5, maxlat=47, minlat=43,
    #                                  land_mask=True)
    
    mdfs = {"low_elev":mdf_low, "high_elev":mdf_high}
    
    return mdfs



# extract for all variables 
d18op_mdfs = estimate_all_transect(varname="d18op", units="per mil")

elev_mdfs = estimate_all_transect(varname="elev", units="m", cal_mean=True)


# t2m_mdfs = estimate_all_transect(varname="temp2", units="°C")
# prec_mdfs = estimate_all_transect(varname="prec", units="mm/month")

# plotting 

apply_style(fontsize=28, style="seaborn-paper", linewidth=2,)

fig,(ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(24, 15))

plot_cross_section(varname="$\delta^{18}$Op vs SMOW", units="‰", data=d18op_mdfs.get("low_elev"), 
                   hue="d18op", xlabel=True, title="Low elevation", ax_legend=True, ymin=-24, ymax=0,
                   ax=ax1, points_data=None,)

plot_cross_section(varname="$\delta^{18}$Op vs SMOW", units="‰", data=d18op_mdfs.get("high_elev"), 
                   hue="d18op", xlabel=True, title="High elevation", ax_legend=False, ymin=-24 , ymax=0,
                   ax=ax2, points_data=None, )

ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='y', which='major', length=8, width=3)
ax1.tick_params(axis='y', which='minor', length=4, width=1.5)

ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='y', which='major', length=8, width=3)
ax2.tick_params(axis='y', which='minor', length=4, width=1.5)


plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
plt.savefig(os.path.join(path_to_plots, "transect_low_high_d18op.svg"), format= "svg", bbox_inches="tight", dpi=600)

        
