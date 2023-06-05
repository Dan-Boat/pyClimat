# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 20:00:07 2023

@author: dboateng
"""

# import modules

import os
import xarray as xr 
import numpy as np
import pandas as pd 


from pyClimat.analysis import extract_var, compute_lterm_diff, compute_lterm_mean
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pyClimat.data import read_from_path
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean

main_path = "D:/Datasets/Model_output_pst"

path_to_store = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

lgm_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
plio_path = os.path.join(main_path, "PLIO", "MONTHLY_MEANS")
mh_path = os.path.join(main_path, "MH", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")
pd_path = os.path.join(main_path, )

filename_lterm = "1003_1017_1m_mlterm.nc"

PI_data = read_from_path(pi_path, filename_lterm)
LGM_data = read_from_path(lgm_path, filename_lterm)
PLIO_data = read_from_path(plio_path, filename_lterm)
MH_data = read_from_path(mh_path, filename_lterm)

filename_lterm = "1003_1017_wiso_1m_mlterm.nc"

PI_wiso = read_from_path(pi_path, filename_lterm)
LGM_wiso = read_from_path(lgm_path, filename_lterm)
PLIO_wiso = read_from_path(plio_path, filename_lterm)
MH_wiso = read_from_path(mh_path, filename_lterm)


def extract_d18Op_anomalies_data(control_data, control_wiso, main_data, main_wiso):
    
    d18Op_control = extract_var(Dataset=control_data, varname="d18op", units="per mil", Dataset_wiso=control_wiso)
    
    d18Op_main = extract_var(Dataset=main_data, varname="d18op", units="per mil", Dataset_wiso=main_wiso)
    
    diff = compute_lterm_diff(data_control=d18Op_control, data_main=d18Op_main,
                              time="month", month="JJAS")
    
    return_data={"control_mon":d18Op_control, "main_mon": d18Op_main, "diff": diff}
    
    return return_data

lgm_diff = extract_d18Op_anomalies_data(control_data=PI_data, control_wiso=PI_wiso, main_data=LGM_data, 
                                        main_wiso=LGM_wiso)
mh_diff = extract_d18Op_anomalies_data(control_data=PI_data, control_wiso=PI_wiso, main_data=MH_data, 
                                        main_wiso=MH_wiso)
plio_diff = extract_d18Op_anomalies_data(control_data=PI_data, control_wiso=PI_wiso, main_data=PLIO_data, 
                                        main_wiso=PLIO_wiso)


apply_style(fontsize=22, style=None, linewidth=2) 

projection = ccrs.PlateCarree()
fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 13), subplot_kw={"projection":
                                                                                                                  projection})
plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=mh_diff.get("diff") , cmap="PiYG", units="‰", ax=ax1, 
                 fig=fig, vmax=6, vmin=-6,levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                           title= ["(a)  MH - PI"], bottom_labels=True, left_labels=True, 
                           compare_data1=mh_diff.get("control_mon"),
                           compare_data2=mh_diff.get("main_mon"), max_pvalue=0.05, plot_stats=True, time="JJAS",
                           hatches=".")

plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=lgm_diff.get("diff") , cmap="PiYG", units="‰", ax=ax2, 
                 fig=fig, vmax=6, vmin=-6,levels=22, domain="West Africa", level_ticks=11, 
                           title= ["(b)  LGM - PI"], bottom_labels=True, left_labels=True, 
                           compare_data1=lgm_diff.get("control_mon"),
                           compare_data2=lgm_diff.get("main_mon"), max_pvalue=0.05, plot_stats=True, time="JJAS",
                           hatches=".", add_colorbar=True, cbar_pos = [0.35, 0.25, 0.25, 0.02], orientation= "horizontal")


plot_annual_mean(variable="$\delta^{18}$Op vs SMOW anomalies", data_alt=plio_diff.get("diff") , cmap="PiYG", units="‰", ax=ax3, 
                 fig=fig, vmax=6, vmin=-6,levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                           title= ["(c)  PLIO - PI"], bottom_labels=True, left_labels=True, 
                           compare_data1=plio_diff.get("control_mon"),
                           compare_data2=plio_diff.get("main_mon"), max_pvalue=0.05, plot_stats=True, time="JJAS",
                           hatches=".")


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.15)
plt.savefig(os.path.join(path_to_store, "d18Op_anomalies.svg"), format= "svg", bbox_inches="tight", dpi=300)

    
    