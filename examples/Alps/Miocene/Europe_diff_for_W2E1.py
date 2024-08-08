# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 16:33:16 2024

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
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var



path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

# experiments
W2E1_Mio278_filename = "a017_hpc-bw_e5w2.3_t159_MIO_W2E1_278ppm_t159l31.6h"
W2E1_Mio450_filename = "a016_hpc-bw_e5w2.3_t159_MIO_W2E1_450ppm_t159l31.6h"



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


def extract_and_analysis(exp_data, exp_wiso, W1E1_data, W1E1_wiso, return_control=True):
    
    # extract prec
    exp_prec_data = extract_var(Dataset=exp_data, varname="prec", units="mm/month")
    w1e1_prec_data = extract_var(Dataset=W1E1_data, varname="prec", units="mm/month")
    simulated_prec_change_alt = compute_lterm_diff(data_control=w1e1_prec_data, data_main=exp_prec_data,
                                          time="annual")
    
    # extract d18op
    exp_d18op_data = extract_var(Dataset=exp_data, varname="d18op", units="per mil", Dataset_wiso=exp_wiso)
    w1e1_d18op_data = extract_var(Dataset=W1E1_data, varname="d18op", units="per mil", Dataset_wiso=W1E1_wiso)
    simulated_d18op_change_alt = compute_lterm_diff(data_control=w1e1_d18op_data, data_main=exp_d18op_data,
                                          time="annual")
    
    # extract temp
    exp_temp_data = extract_var(Dataset=exp_data, varname="temp2", units="°C")
    w1e1_temp_data = extract_var(Dataset=W1E1_data, varname="temp2", units="°C")
    simulated_temp_change_alt = compute_lterm_diff(data_control=w1e1_temp_data, data_main=exp_temp_data,
                                          time="annual")
    
    
   
    return_data= {"d18O": simulated_d18op_change_alt, "temp": simulated_temp_change_alt,
                  "prec": simulated_prec_change_alt}
    
    return_data_exp_mon = {"d18O_mon": exp_d18op_data, "prec_mon": exp_prec_data,
                           "temp_mon": exp_temp_data}
    if return_control:
        return_data_control_mon = {"d18O_mon": w1e1_d18op_data, "prec_mon": w1e1_prec_data,
                               "temp_mon": w1e1_temp_data}
        return return_data, return_data_exp_mon, return_data_control_mon
    else:
        return return_data, return_data_exp_mon
    
    

W2E1_278_diff, W2E1_278_mon, W1E1_278_mon = extract_and_analysis(exp_data=W2E1_278_data, exp_wiso=W2E1_278_wiso, W1E1_data=W1E1_278_data,
                                            W1E1_wiso=W1E1_278_wiso, return_control=True)
W2E1_450_diff, W2E1_450_mon, W1E1_450_mon = extract_and_analysis(exp_data=W2E1_450_data, exp_wiso=W2E1_450_wiso, W1E1_data=W1E1_450_data,
                                            W1E1_wiso=W1E1_450_wiso, return_control=True)


apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)

fig,((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(22,25), 
                                                                subplot_kw={"projection":projection})

plot_annual_mean(variable="$\delta^{18}$Op VSMOW difference", data_alt=W2E1_278_diff.get("d18O"), 
                 ax=ax1, cmap=GryBr_r, units="‰", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.95, 0.65, 0.02, 0.20], 
                orientation="vertical", plot_coastlines=False, bottom_labels=False,
                left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, domain="Europe", compare_data1=W2E1_278_mon.get("d18O_mon"), 
                compare_data2=W1E1_278_mon.get("d18O_mon"), max_pvalue=0.1, plot_stats=True,
                hatches=".", title="(a) MIO 278 ppm (W2E1-W1E1)")

plot_annual_mean(variable="$\delta^{18}$Op VSMOW difference", data_alt=W2E1_450_diff.get("d18O"), 
                 ax=ax2, cmap=GryBr_r, units="‰", vmax=10, vmin=-10, 
                levels=22, level_ticks=11, add_colorbar=False, 
                plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, domain="Europe", compare_data1=W2E1_450_mon.get("d18O_mon"), 
                compare_data2=W1E1_450_mon.get("d18O_mon"), max_pvalue=0.1, plot_stats=True,
                hatches=".", title="(b) MIO 450 ppm (W2E1-W1E1)")

plot_annual_mean(variable="Temperature difference", data_alt=W2E1_278_diff.get("temp"), 
                 ax=ax3, cmap=RdBu_r, units="‰", vmax=3, vmin=-3, 
                levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.95, 0.35, 0.02, 0.20], 
                orientation="vertical", plot_coastlines=False, bottom_labels=False,
                left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, domain="Europe", compare_data1=W2E1_278_mon.get("temp_mon"), 
                compare_data2=W1E1_278_mon.get("temp_mon"), max_pvalue=0.1, plot_stats=True,
                hatches=".", title="(c)")

plot_annual_mean(variable="Temperature difference", data_alt=W2E1_450_diff.get("temp"), 
                 ax=ax4, cmap=RdBu_r, units="‰", vmax=3, vmin=-3, 
                levels=22, level_ticks=11, add_colorbar=False, 
                plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, domain="Europe", compare_data1=W2E1_450_mon.get("temp_mon"), 
                compare_data2=W1E1_450_mon.get("temp_mon"), max_pvalue=0.1, plot_stats=True,
                hatches=".", title="(d)")

plot_annual_mean(variable="Precipitation difference", data_alt=W2E1_278_diff.get("prec"), 
                 ax=ax5, cmap=PrecAno, units="‰", vmax=60, vmin=-60, 
                levels=22, level_ticks=11, add_colorbar=True, cbar_pos= [0.95, 0.05, 0.02, 0.20], 
                orientation="vertical", plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, domain="Europe", compare_data1=W2E1_278_mon.get("prec_mon"), 
                compare_data2=W1E1_278_mon.get("prec_mon"), max_pvalue=0.1, plot_stats=True,
                hatches=".", title="(e)")

plot_annual_mean(variable="Precipitation difference", data_alt=W2E1_450_diff.get("prec"), 
                 ax=ax6, cmap=PrecAno, units="‰", vmax=60, vmin=-60, 
                levels=22, level_ticks=11, add_colorbar=False, 
                plot_coastlines=False, bottom_labels=True,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                plot_projection=projection, domain="Europe", compare_data1=W2E1_450_mon.get("prec_mon"), 
                compare_data2=W1E1_450_mon.get("prec_mon"), max_pvalue=0.1, plot_stats=True,
                hatches=".", title="(f)")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
plt.savefig(os.path.join(path_to_plots, "europe_diff.svg"), format= "svg", bbox_inches="tight", dpi=600)


