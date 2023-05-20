# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:35:57 2023

@author: dboateng
"""

# import models
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
import seaborn as sns


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import scatter_plot_laspe_rate
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, extract_transect, linregression


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"


# read data

W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
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

W1E1_PI_data, W1E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_PI_filename, years=years, 
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



def extract_vars_and_analysis(data, wiso, maxlon=5, maxlat=46, minlon=-1, minlat=43):

    d18op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= wiso)
    elev = extract_var(Dataset=data , varname="elev", units="m")

    # compute annual means
    d18op_alt = compute_lterm_mean(data=d18op, time="annual")
    elev_alt = compute_lterm_mean(data=elev, time="annual")
    
    elev_transect = extract_transect(data=elev_alt, maxlon=maxlon, minlon=minlon, 
                                  maxlat=maxlat, minlat=minlat, sea_land_mask=True, 
                                  Dataset=data)
    
    
    
    d18op_transect = extract_transect(data=d18op_alt, maxlon=maxlon, minlon=minlon, 
                                  maxlat=maxlat, minlat=minlat, sea_land_mask=True, 
                                  Dataset=data)
   
   
    reg, df = linregression(data_x=elev_transect, data_y=d18op_transect, return_yhat=True)
    
    
    return_data = {"reg":reg, "df":df}
    
    return return_data

def extract_label(data):
    # labels
    slope_pi = data.get("reg").slope*1000
    coef_err_pi = data.get("reg").stderr*1000
    r2_pi = data.get("reg").rvalue*-1
    
    label = "{:.2f} +/-{:.2f} [‰/km], r²={:.2f}".format(slope_pi, coef_err_pi,r2_pi)
    
    return label


def plot_section_lapse_rate(data_pi, data_278, data_450, ax=None, xlabel=True, ylabel=True, title=None, 
                            ax_legend=True, ymin=None, ymax=None,):
    
   
    

    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
    
                                                              
    # plot the regression line and points
    label_pi = extract_label(data_pi)
    sns.regplot(data=data_pi.get("df"), x="X", y="Y", marker="*", 
               scatter_kws={"color":"black", "s":200}, color="black",ax=ax,
               label=label_pi)
    
    label_278 = extract_label(data_278)
    sns.regplot(data=data_278.get("df"), x="X", y="Y", marker="D", 
               scatter_kws={"color":"red", "s":200},
               color ="red", label=label_278, ax=ax)
    
    
    label_450 = extract_label(data_450)
    sns.regplot(data=data_450.get("df"), x="X", y="Y", marker="^", 
               scatter_kws={"color":"green", "s":200},
               color ="green", label=label_450, ax=ax)
    
    
    # add the labels (with pvalues estimates)
            
            
    if ylabel:
        ax.set_ylabel(u'$\delta^{18}$O ‰ vs SMOW (‰)', fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_yticklabels([])
    
    if xlabel:
        ax.set_xlabel("Elevation (m)", fontweight="bold", fontsize=28)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_xticklabels([])
        
    if all(parameter is not None for parameter in [ymax, ymin]):
        ax.set_ylim(ymin, ymax)
        
    if ax_legend:
        ax.legend(frameon=True, fontsize=28,
                  loc="lower left", borderaxespad=0,
                  ncol=1) #bbox_to_anchor=(0.01, 1.05, 1, 0.102,)
    else:
        ax.legend([],[], frameon=False)
        
    ax.set_box_aspect(1)
    
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 28, "fontweight":"bold"}, loc="left")
        
    
        
        
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    


def plot_west_transect():

    W1E1_PI_stats = extract_vars_and_analysis(data=W1E1_PI_data, wiso=W1E1_PI_wiso, 
                                            maxlon=8, maxlat=47, minlon=1, minlat=44)
    W1E1_278_stats = extract_vars_and_analysis(data= W1E1_278_data, wiso=W1E1_278_wiso)
    W1E1_450_stats = extract_vars_and_analysis(data= W1E1_450_data, wiso=W1E1_450_wiso)
    
    W2E0_PI_stats = extract_vars_and_analysis(data=W2E0_PI_data, wiso=W2E0_PI_wiso, 
                                            maxlon=8, maxlat=47, minlon=1, minlat=44)
    W2E0_278_stats = extract_vars_and_analysis(data= W2E0_278_data, wiso=W2E0_278_wiso)
    W2E0_450_stats = extract_vars_and_analysis(data= W2E0_450_data, wiso=W2E0_450_wiso)
    
    
    W2E1_PI_stats = extract_vars_and_analysis(data=W2E1_PI_data, wiso=W2E1_PI_wiso, 
                                            maxlon=8, maxlat=47, minlon=1, minlat=44)
    W2E1_278_stats = extract_vars_and_analysis(data= W2E1_278_data, wiso=W2E1_278_wiso)
    W2E1_450_stats = extract_vars_and_analysis(data= W2E1_450_data, wiso=W2E1_450_wiso)
    
    
    W2E2_PI_stats = extract_vars_and_analysis(data=W2E2_PI_data, wiso=W2E2_PI_wiso, 
                                            maxlon=8, maxlat=47, minlon=1, minlat=44)
    W2E2_278_stats = extract_vars_and_analysis(data= W2E2_278_data, wiso=W2E2_278_wiso)
    W2E2_450_stats = extract_vars_and_analysis(data= W2E2_450_data, wiso=W2E2_450_wiso)
    
    
    apply_style(fontsize=28, style="seaborn-paper", linewidth=2,)
    
    fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(24, 22))
    
    
    plot_section_lapse_rate(data_pi=W1E1_PI_stats, data_278=W1E1_278_stats, data_450=W1E1_450_stats,
                            ax=ax1, title="(a) W1E1")
    
    plot_section_lapse_rate(data_pi=W2E0_PI_stats, data_278=W2E0_278_stats, data_450=W2E0_450_stats,
                            ax=ax2, title="(b) W2E0")
    
    plot_section_lapse_rate(data_pi=W2E1_PI_stats, data_278=W2E1_278_stats, data_450=W2E1_450_stats,
                            ax=ax3, title="(c) W2E1")
    
    plot_section_lapse_rate(data_pi=W2E2_PI_stats, data_278=W2E2_278_stats, data_450=W2E2_450_stats,
                            ax=ax4, title="(d) W2E2")
    
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_plots, "lapse_rate_west.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
    
def plot_north_transect():

    W1E1_PI_stats = extract_vars_and_analysis(data=W1E1_PI_data, wiso=W1E1_PI_wiso, 
                                            maxlon=16, maxlat=50, minlon=5, minlat=46.5)
    W1E1_278_stats = extract_vars_and_analysis(data= W1E1_278_data, wiso=W1E1_278_wiso,
                                               maxlon=15, maxlat=49, minlon=4, minlat=45.5)
    W1E1_450_stats = extract_vars_and_analysis(data= W1E1_450_data, wiso=W1E1_450_wiso,
                                               maxlon=15, maxlat=49, minlon=4, minlat=45.5)
    
    W2E0_PI_stats = extract_vars_and_analysis(data=W2E0_PI_data, wiso=W2E0_PI_wiso, 
                                            maxlon=16, maxlat=50, minlon=5, minlat=46.5)
    W2E0_278_stats = extract_vars_and_analysis(data= W2E0_278_data, wiso=W2E0_278_wiso,
                                               maxlon=15, maxlat=49, minlon=4, minlat=45.5)
    W2E0_450_stats = extract_vars_and_analysis(data= W2E0_450_data, wiso=W2E0_450_wiso,
                                               maxlon=15, maxlat=49, minlon=4, minlat=45.5)
    
    
    W2E1_PI_stats = extract_vars_and_analysis(data=W2E1_PI_data, wiso=W2E1_PI_wiso, 
                                            maxlon=16, maxlat=50, minlon=5, minlat=46.5)
    W2E1_278_stats = extract_vars_and_analysis(data= W2E1_278_data, wiso=W2E1_278_wiso,
                                               maxlon=15, maxlat=49, minlon=4, minlat=45.5)
    W2E1_450_stats = extract_vars_and_analysis(data= W2E1_450_data, wiso=W2E1_450_wiso,
                                               maxlon=15, maxlat=49, minlon=4, minlat=45.5)
    
    
    W2E2_PI_stats = extract_vars_and_analysis(data=W2E2_PI_data, wiso=W2E2_PI_wiso, 
                                            maxlon=16, maxlat=50, minlon=5, minlat=46.5)
    W2E2_278_stats = extract_vars_and_analysis(data= W2E2_278_data, wiso=W2E2_278_wiso,
                                               maxlon=15, maxlat=49, minlon=4, minlat=45.5)
    W2E2_450_stats = extract_vars_and_analysis(data= W2E2_450_data, wiso=W2E2_450_wiso,
                                               maxlon=15, maxlat=49, minlon=4, minlat=45.5)
    
    
    apply_style(fontsize=28, style="seaborn-paper", linewidth=2,)
    
    fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(24, 22))
    
    
    plot_section_lapse_rate(data_pi=W1E1_PI_stats, data_278=W1E1_278_stats, data_450=W1E1_450_stats,
                            ax=ax1, title="(a) W1E1")
    
    plot_section_lapse_rate(data_pi=W2E0_PI_stats, data_278=W2E0_278_stats, data_450=W2E0_450_stats,
                            ax=ax2, title="(b) W2E0")
    
    plot_section_lapse_rate(data_pi=W2E1_PI_stats, data_278=W2E1_278_stats, data_450=W2E1_450_stats,
                            ax=ax3, title="(c) W2E1")
    
    plot_section_lapse_rate(data_pi=W2E2_PI_stats, data_278=W2E2_278_stats, data_450=W2E2_450_stats,
                            ax=ax4, title="(d) W2E2")
    
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_plots, "lapse_rate_north.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
    
def plot_south_transect():

    W1E1_PI_stats = extract_vars_and_analysis(data=W1E1_PI_data, wiso=W1E1_PI_wiso, 
                                            maxlon=15, maxlat=47, minlon=7.5, minlat=43)
    W1E1_278_stats = extract_vars_and_analysis(data= W1E1_278_data, wiso=W1E1_278_wiso,
                                               maxlon=12.5, maxlat=46, minlon=5, minlat=42)
    W1E1_450_stats = extract_vars_and_analysis(data= W1E1_450_data, wiso=W1E1_450_wiso,
                                               maxlon=12.5, maxlat=46, minlon=5, minlat=42)
    
    W2E0_PI_stats = extract_vars_and_analysis(data=W2E0_PI_data, wiso=W2E0_PI_wiso, 
                                            maxlon=15, maxlat=47, minlon=7.5, minlat=43)
    W2E0_278_stats = extract_vars_and_analysis(data= W2E0_278_data, wiso=W2E0_278_wiso,
                                               maxlon=12.5, maxlat=46, minlon=5, minlat=42)
    W2E0_450_stats = extract_vars_and_analysis(data= W2E0_450_data, wiso=W2E0_450_wiso,
                                               maxlon=12.5, maxlat=46, minlon=5, minlat=42)
    
    
    W2E1_PI_stats = extract_vars_and_analysis(data=W2E1_PI_data, wiso=W2E1_PI_wiso, 
                                            maxlon=15, maxlat=47, minlon=7.5, minlat=43)
    W2E1_278_stats = extract_vars_and_analysis(data= W2E1_278_data, wiso=W2E1_278_wiso,
                                               maxlon=12.5, maxlat=46, minlon=5, minlat=42)
    W2E1_450_stats = extract_vars_and_analysis(data= W2E1_450_data, wiso=W2E1_450_wiso,
                                               maxlon=12.5, maxlat=46, minlon=5, minlat=42)
    
    
    W2E2_PI_stats = extract_vars_and_analysis(data=W2E2_PI_data, wiso=W2E2_PI_wiso, 
                                            maxlon=15, maxlat=47, minlon=7.5, minlat=43)
    W2E2_278_stats = extract_vars_and_analysis(data= W2E2_278_data, wiso=W2E2_278_wiso,
                                               maxlon=12.5, maxlat=46, minlon=5, minlat=42)
    W2E2_450_stats = extract_vars_and_analysis(data= W2E2_450_data, wiso=W2E2_450_wiso,
                                               maxlon=12.5, maxlat=46, minlon=5, minlat=42)
    
    
    apply_style(fontsize=28, style="seaborn-paper", linewidth=2,)
    
    fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(24, 22))
    
    
    plot_section_lapse_rate(data_pi=W1E1_PI_stats, data_278=W1E1_278_stats, data_450=W1E1_450_stats,
                            ax=ax1, title="(a) W1E1")
    
    plot_section_lapse_rate(data_pi=W2E0_PI_stats, data_278=W2E0_278_stats, data_450=W2E0_450_stats,
                            ax=ax2, title="(b) W2E0")
    
    plot_section_lapse_rate(data_pi=W2E1_PI_stats, data_278=W2E1_278_stats, data_450=W2E1_450_stats,
                            ax=ax3, title="(c) W2E1")
    
    plot_section_lapse_rate(data_pi=W2E2_PI_stats, data_278=W2E2_278_stats, data_450=W2E2_450_stats,
                            ax=ax4, title="(d) W2E2")
    
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_plots, "lapse_rate_south.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
if __name__ == "__main__":
    #plot_west_transect()
    plot_north_transect()
    #plot_south_transect()