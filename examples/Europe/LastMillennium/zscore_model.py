# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 16:03:59 2024

@author: dboateng
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs
import seaborn as sns

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_transect, compute_lterm_diff, compute_lterm_mean
from pyClimat.variables import extract_var
from pyClimat.utils import extract_region


from path_to_data_lm import *

main_path = "D:/Datasets/iGCM_datasets/"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"

filenames = ["CESM", "ECHAM5", "GISS", "HADCM3", "CCSM"]



def extract_climatologies(filename, time="annual", season=None):
    
    filename_data = filename + "_wiso_vars.nc"
    
    
    d18O = read_from_path(path=main_path, filename=filename_data, decode=True,
                      varname="d18O")  
    
    data  = d18O.sel(time= slice("0950-01", "1950-01"))
    
    mean = data.mean(dim="time", skipna=True)
    std = data.std(dim="time", skipna=True)
    data = (data - mean) / std
    
    
    date_range_mca = d18O.sel(time= slice("0950-01", "1250-01")).time
    date_range_lia = d18O.sel(time= slice("1650-01", "1850-01")).time
    date_range_pi = d18O.sel(time= slice("1851-01", "1950-01")).time
        
        
    
    d18O_mca = compute_lterm_mean(data=data, time=time, time_range=date_range_mca, 
                                  season=season)
    d18O_lia = compute_lterm_mean(data=data, time=time, time_range=date_range_lia, 
                                  season=season)
    
    d18O_pi = compute_lterm_mean(data=data, time=time, time_range=date_range_pi, 
                                  season=season)
    
    return d18O_mca, d18O_lia, d18O_pi


# test function 
# cesm_mca, cesm_lia, cesm_diff = extract_climatologies(filename="CESM", time="annual")
# giss_mca, giss_lia, giss_diff = extract_climatologies(filename="GISS", time="annual")
data_mca, data_lia, data_pi = extract_climatologies(filename="HADCM3", time="annual")

regions = ["Greenland", "South EU", "North EU"]

def extract_regions_from_data(data):
    greenland = extract_transect(data=data, maxlon=-15, minlon=-90, maxlat=85, minlat=60)
    south_eu = extract_transect(data=data, maxlon=35, minlon=-15, maxlat=50, minlat=30)
    north_eu = extract_transect(data=data, maxlon=35, minlon=-15, maxlat=85, minlat=48)

    return [greenland.values.ravel(), south_eu.values.ravel(), north_eu.values.ravel()]


def extract_all(data_mca, data_lia, data_pi):
    
    mca_regions = extract_regions_from_data(data_mca)
    lia_regions = extract_regions_from_data(data_lia)
    pi_regions = extract_regions_from_data(data_pi)
    
    
    max_len = max(len(m) for m in mca_regions)

    df_mca = pd.DataFrame(index=np.arange(max_len), columns=regions).assign(Paleoclimate="MCA")

    for i,r in enumerate(regions):
        current_len = len(mca_regions[i])
        df_mca.loc[:current_len - 1, r] = mca_regions[i]
        
        
    
    max_len = max(len(l) for l in lia_regions)

    df_lia = pd.DataFrame(index=np.arange(max_len), columns=regions).assign(Paleoclimate="LIA")

    for i,r in enumerate(regions):
        current_len = len(lia_regions[i])
        df_lia.loc[:current_len - 1, r] = lia_regions[i]
    
   
    
    max_len = max(len(p) for p in pi_regions)

    df_pi = pd.DataFrame(index=np.arange(max_len), columns=regions).assign(Paleoclimate="PI")

    for i,r in enumerate(regions):
        current_len = len(pi_regions[i])
        df_pi.loc[:current_len - 1, r] = pi_regions[i]
        
    
    cdf = pd.concat([df_mca, df_lia, df_pi])
    mdf = pd.melt(cdf, id_vars=["Paleoclimate"],var_name="d18O")
    mdf["value"] = mdf["value"].astype("float")
    
    mdf_mean = mdf.groupby(["Paleoclimate", "d18O"]).mean()
    
    return mdf, mdf_mean 
    
    
def plot_distribution(varname, units, data, hue, ax=None, path_to_plots=None, filename=None,
                       colors=None, xlabel=True, ylabel=True, title=None, ax_legend=True,
                       ymin=None, ymax=None, points_data=None, violin_plot=False):
    
    regions = ["Greenland", "North EU", "South EU"]
    
    apply_style(fontsize=25, style=None, linewidth=3)
    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
    if violin_plot:    
        boxplot = sns.violinplot(data=data, x="Paleoclimate", y="value", saturation=0.8, ax=ax,
                          hue=hue, density_norm="width", fill=False, linewidth=3) # count, width, area
        
    else:
        boxplot = sns.boxplot(data=data, x="Paleoclimate", y="value", saturation=0.7, ax=ax,
                              hue=hue, linewidth=3, ) # count, width, area
    
    if points_data is not None:
        
        pointplot = sns.pointplot(data=points_data, x="Paleoclimate", y="value", ax=ax,
                              hue=hue, markers="x", linestyles="", scale=2,
                              dodge=True, legend=False,)
    
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
        ax.set_xlabel("Paleoclimate", fontweight="bold", fontsize=28)
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
    
    if path_to_plots is not None:
        plt.savefig(os.path.join(path_to_plots, filename), bbox_inches="tight", format= "png")
    
   
    
    



mdf, mdf_mean = extract_all(data_mca, data_lia, data_pi)  

apply_style(fontsize=25, style=None, linewidth=3)

plot_distribution(varname='$\delta^{18}$Op VSMOW', units="‰", data=mdf, hue="d18O", ymax=0.25, ymin=-0.1,
                  ax_legend=True,xlabel=True, title="iHadCM3 (0950-01-01 to 1950-12-31)",
                  path_to_plots=path_to_plots, filename="model_regions.png", points_data=None,
                  )  
    
    