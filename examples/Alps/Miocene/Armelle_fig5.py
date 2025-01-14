# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 14:30:15 2024

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
from pyClimat.data import read_ECHAM_processed, read_from_path
from pyClimat.analysis import compute_lterm_mean, extract_transect
from pyClimat.variables import extract_var


main_path = "D:/Datasets/Model_output_pst/PD"

path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/armelle"


#load datasets 
PD_data = read_from_path(main_path, "PD_1980_2014_monthly.nc", decode=True)
PD_wiso = read_from_path(main_path, "PD_1980_2014_monthly_wiso.nc", decode=True)


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"
#W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"

# reading data 
# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)

def weighted_mean(data):
    weights = np.cos(np.deg2rad(data.lat))
    weights.name = "weights"
    
    data_weighted = data.weighted(weights)
    
    data_mean = data_weighted.mean(dim=("lat", "lon"), skipna=True)
    return data_mean

def compute_mean(data, save=False, path_to_save=None):
    data_mean = weighted_mean(data)
    
    
    import calendar
    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    df = pd.DataFrame(index = mnames, columns = ["mean"])
    df["mean"] = data_mean
    
    if save == True:
        df.to_csv(path_to_save)
    return df


# # extract values and compute long-term means

def extract_vars_and_analysis(data, wiso):
    
    d18Op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    
    t2m = extract_var(Dataset=data, varname="temp2", units="°C")
    prec = extract_var(Dataset=data, varname="prec", units="mm/month")
    
    d18Op_alt = compute_lterm_mean(data=d18Op, time="month")
    prec_alt = compute_lterm_mean(data=prec, time="month") 
    t2m_alt = compute_lterm_mean(data=t2m, time="month")
    
    return d18Op_alt, prec_alt, t2m_alt



    


def extract_sections(data):
    
    max_lon, max_lat, min_lon, min_lat = 9.80, 47.4, 6.66, 46.6
    
    extract_nafb = extract_transect(data=data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat)
    
    df_nafb = compute_mean(extract_nafb)
    
    max_lon, max_lat, min_lon, min_lat = 7.4, 44.7, 4.9, 42.7
    
    extract_dvb = extract_transect(data=data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat)
    
    df_dvb = compute_mean(extract_dvb)
    
    
    return df_nafb, df_dvb



def plot_monthly_period_per_section(PD, Mio278, Mio450, ax, title, ymax, ymin, varname=None):
    
    ax.plot(PD["mean"], "-", color="black", label="PD", linewidth=4)
    ax.plot(Mio278["mean"], "-", color="blue", label="MIO 278 ppm", linewidth=4)
    ax.plot(Mio450["mean"], "-", color="red", label="MIO 450 ppm", linewidth=4)
    
    
    if varname is not None:
        ax.set_ylabel(varname, fontweight="bold", fontsize=22)
    
    #ax.grid(False, linestyle="--", color=grey, alpha=0.8)
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc= "left")
    ax.axes.tick_params(which="both", labelsize=25)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
    
# read data 
Mio278_d18Op, Mio278_prec, Mio278_t2m = extract_vars_and_analysis(data=W1E1_278_data, wiso=W1E1_278_wiso)
Mio450_d18Op, Mio450_prec,  Mio450_t2m = extract_vars_and_analysis(data=W1E1_450_data, wiso=W1E1_450_wiso)
PD_d18Op, PD_prec, PD_t2m = extract_vars_and_analysis(data=PD_data, wiso=PD_wiso)



Mio278_d18Op_nafb, Mio278_d18Op_dvb = extract_sections(Mio278_d18Op)
Mio450_d18Op_nafb, Mio450_d18Op_dvb = extract_sections(Mio450_d18Op)
PD_d18Op_nafb, PD_d18Op_dvb = extract_sections(PD_d18Op)

apply_style(fontsize=24, style="seaborn-poster", linewidth=3) 
fig, ((ax1,ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(15, 18), sharey=False, sharex=True)

plot_monthly_period_per_section(PD=PD_d18Op_nafb, Mio278=Mio278_d18Op_nafb, Mio450=Mio450_d18Op_nafb,
                                ax=ax1, title="(a) NAFB", ymax=-2, ymin=-14, 
                                varname="$\delta^{18}$Op VSMOW (‰)")

plot_monthly_period_per_section(PD=PD_d18Op_dvb, Mio278=Mio278_d18Op_dvb, Mio450=Mio450_d18Op_dvb,
                                ax=ax2, title="(b) DVB", ymax=-2, ymin=-14, 
                                varname=None)


Mio278_prec_nafb, Mio278_prec_dvb = extract_sections(Mio278_prec)
Mio450_prec_nafb, Mio450_prec_dvb = extract_sections(Mio450_prec)
PD_prec_nafb, PD_prec_dvb = extract_sections(PD_prec)


plot_monthly_period_per_section(PD=PD_prec_nafb, Mio278=Mio278_prec_nafb, Mio450=Mio450_prec_nafb,
                                ax=ax3, title="(c)", ymax=180, ymin=0, 
                                varname="Precipitation [mm/month]")

plot_monthly_period_per_section(PD=PD_prec_dvb, Mio278=Mio278_prec_dvb, Mio450=Mio450_prec_dvb,
                                ax=ax4, title="(d)", ymax=180, ymin=0, 
                                varname=None)

Mio278_t2m_nafb, Mio278_t2m_dvb = extract_sections(Mio278_t2m)
Mio450_t2m_nafb, Mio450_t2m_dvb = extract_sections(Mio450_t2m)
PD_t2m_nafb, PD_t2m_dvb = extract_sections(PD_t2m)


plot_monthly_period_per_section(PD=PD_t2m_nafb, Mio278=Mio278_t2m_nafb, Mio450=Mio450_t2m_nafb,
                                ax=ax5, title="(e)", ymax=30, ymin=2, 
                                varname="Temperature (°C)")

plot_monthly_period_per_section(PD=PD_t2m_dvb, Mio278=Mio278_t2m_dvb, Mio450=Mio450_t2m_dvb,
                                ax=ax6, title="(f)", ymax=30, ymin=2, 
                                varname=None)
    
axes = [ax1,ax2, ax3, ax4, ax5, ax6]
for ax in axes:
    ax.grid(True, linestyle="--", color="grey")
plt.tight_layout()
plt.savefig(os.path.join(path_to_plots, "monthly_sections.pdf"), format= "pdf", bbox_inches="tight", dpi=300)
