# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 10:13:05 2023

@author: dboateng
"""

import os
import xarray as xr 
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates
import numpy as np

from pyClimat.plot_utils import *
from pyClimat.plots import plot_causal_statistics

from paths_to_data import *


filename = "NAO_caused_by_d18op_test.nc"

filenames_NAO = ["NAO_caused_by_d18op.nc", "NAO_caused_by_t2m.nc", "NAO_caused_by_prec.nc",
                 "NAO_caused_by_d18op_reverse.nc", "NAO_caused_by_t2m_reverse.nc", 
                 "NAO_caused_by_prec_reverse.nc"]

labels_NAO = ["(a) Y($\delta^{18}$Op) to X(NAO)", "(b) Y(Temperature) to X(NAO)", "(c) Y(Precipitation) to X(NAO)",
              "(d) Y(NAO) to X($\delta^{18}$Op)", "(e) Y(NAO) to X(Temperature)", "(f) Y(NAO) to X(Precipitation)"]



filenames_EA = ["EA_caused_by_d18op.nc", "EA_caused_by_t2m.nc", "EA_caused_by_prec.nc",
                 "EA_caused_by_d18op_reverse.nc", "EA_caused_by_t2m_reverse.nc", 
                 "EA_caused_by_prec_reverse.nc"]

labels_EA = ["(a) Y($\delta^{18}$Op) to X(EA)", "(b) Y(Temperature) to X(EA)", "(c) Y(Precipitation) to X(EA)",
              "(d) Y(EA) to X($\delta^{18}$Op)", "(e) Y(EA) to X(Temperature)", "(f) Y(EA) to X(Precipitation)"]


filenames_NAO_EA = ["d18op_caused_by_NAO_EA.nc", "t2m_caused_by_NAO_EA.nc", "prec_caused_by_NAO_EA.nc"]

labels_NAO_EA = ["(a) Y(NAO,EA) to X($\delta^{18}$Op)", "(b) Y(NAO,EA) to X(Temperature)", "(c) Y(NAO,EA) to X(Precipitation)"]

def plot_causal(filenames, labels, figname):

    apply_style(fontsize=22, style=None, linewidth=2) 
    
    projection = ccrs.PlateCarree()
    
    if len(labels) ==6:
    
        fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(28, 22), 
                                                                                  subplot_kw={"projection": projection})
        axes = [ax1,ax2, ax3, ax4, ax5, ax6]
        
    elif len(labels) == 3:
        fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(20, 12), 
                                                                                  subplot_kw={"projection": projection})
        axes = [ax1,ax2, ax3]
        
    
    for i,filename in enumerate(filenames):
        
        data = xr.open_dataarray(os.path.join(main_path_to_data, filename))
        
        if i == 0:
            plot_causal_statistics(variable="p-value", data=data, ax=axes[i], cmap="Blues", levels=20,
                                   domain = "Europe Wide", vmin=0, vmax=1, level_ticks=8, cbar=True, 
                                   cbar_orientation="horizontal", cbar_pos=[0.40, 0.05, 0.25, 0.02],fig=fig,
                                   title=labels[i])
            
        else:
            plot_causal_statistics(variable="p-value", data=data, ax=axes[i], cmap="Blues", levels=20,
                                   domain = "Europe Wide", vmin=0, vmax=1, level_ticks=8, cbar=False, 
                                   fig=fig, title=labels[i])
            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
    plt.savefig(os.path.join(path_to_plots, figname +".svg"), format= "svg", bbox_inches="tight", dpi=300)
    


plot_causal(filenames=filenames_NAO, labels=labels_NAO, figname="causal_NAO_climate")
plot_causal(filenames=filenames_EA, labels=labels_EA, figname="causal_EA_climate")
plot_causal(filenames=filenames_NAO_EA, labels=labels_NAO_EA, figname="causal_NAO_EA_climate")