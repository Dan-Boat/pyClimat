# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:57:16 2024

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

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean, plot_correlation


path_to_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots/data"
path_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"

cor_cesm_path = os.path.join(path_to_data, "CESM_corr_data.nc")
cor_giss_path = os.path.join(path_to_data, "GISS_corr_data.nc")
# read data

cesm_data = xr.open_dataset(cor_cesm_path)
giss_data = xr.open_dataset(cor_giss_path)

apply_style(fontsize=28, style=None, linewidth=2.5) 
        
projection = ccrs.Robinson(central_longitude=0, globe=None)
def plot_d18O():
    fig,((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(24,26), subplot_kw={"projection":projection})
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_d18O"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=True, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax1, fig=fig, cbar_pos= [0.35, 0.01, 0.45, 0.02],
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_d18O"], bottom_labels=True, 
                     title=" (a) cor(NAO-$\delta^{18}$Op VSMOW) (All): CESM", plot_projection=projection,
                     extend="both")
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_d18O"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax2, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_d18O"], bottom_labels=True, 
                     title=" (b) cor(NAO-$\delta^{18}$Op VSMOW) (All): GISS", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_d18O_EQ"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax3, fig=fig,
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_d18O_EQ"], bottom_labels=True, 
                     title=" (c) cor(NAO-$\delta^{18}$Op VSMOW) (EQ): CESM", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_d18O_EQ"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax4, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_d18O_EQ"], bottom_labels=True, 
                     title=" (d) cor(NAO-$\delta^{18}$Op VSMOW) (EQ): GISS", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_d18O_OP"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax5, fig=fig,
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_d18O_OP"], bottom_labels=True, 
                     title=" (e) cor(NAO-$\delta^{18}$Op VSMOW) (OP): CESM", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_d18O_OP"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax6, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_d18O_OP"], bottom_labels=True, 
                     title=" (f) cor(NAO-$\delta^{18}$Op VSMOW) (OP): GISS", plot_projection=projection,)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.10, wspace=0.08)
    plt.savefig(os.path.join(path_plots, "correlation_d18O.png"), format= "png", bbox_inches="tight", dpi=600)


def plot_prec():
    fig,((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(24,26), subplot_kw={"projection":projection})
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_prec"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=True, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax1, fig=fig, cbar_pos= [0.35, 0.01, 0.45, 0.02],
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_prec"], bottom_labels=True, 
                     title=" (a) cor(NAO-Precipitation) (All): CESM", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_prec"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax2, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_prec"], bottom_labels=True, 
                     title=" (b) cor(NAO-Precipitation) (All): GISS", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_prec_EQ"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax3, fig=fig,
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_prec_EQ"], bottom_labels=True, 
                     title=" (c) cor(NAO-Precipitation) (EQ): CESM", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_prec_EQ"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax4, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_prec_EQ"], bottom_labels=True, 
                     title=" (d) cor(NAO-Precipitation) (EQ): GISS", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_prec_OP"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax5, fig=fig,
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_prec_OP"], bottom_labels=True, 
                     title=" (e) cor(NAO-Precipitation) (OP): CESM", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_prec_OP"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax6, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_prec_OP"], bottom_labels=True, 
                     title=" (f) cor(NAO-Precipitation) (OP): GISS", plot_projection=projection,)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.10, wspace=0.08)
    plt.savefig(os.path.join(path_plots, "correlation_prec.png"), format= "png", bbox_inches="tight", dpi=600)


def plot_temp():
    fig,((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(24,26), subplot_kw={"projection":projection})
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_temp"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=True, cbar_orientation="horizontal",
                     level_ticks=7, ax=ax1, fig=fig, cbar_pos= [0.35, 0.01, 0.45, 0.02],
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_temp"], bottom_labels=True, 
                     title=" (a) cor(NAO-Temperature) (All): CESM", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_temp"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax2, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_temp"], bottom_labels=True, 
                     title=" (b) cor(NAO-Temperature) (All): GISS", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_temp_EQ"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax3, fig=fig,
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_temp_EQ"], bottom_labels=True, 
                     title=" (c) cor(NAO-Temperature) (EQ): CESM", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_temp_EQ"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax4, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_temp_EQ"], bottom_labels=True, 
                     title=" (d) cor(NAO-Temperature) (EQ): GISS", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=cesm_data["nao_reg_temp_OP"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax5, fig=fig,
                     plot_pvalues=True, pvalue_data=cesm_data["nao_sig_temp_OP"], bottom_labels=True, 
                     title=" (e) cor(NAO-Temperature) (OP): CESM", plot_projection=projection,)
    
    plot_correlation(variable="Spearman Coefficients", data=giss_data["nao_reg_temp_OP"], units="-", vmax=0.8,
                     vmin=-0.8, cmap="PRGn", domain="NH", levels=22,cbar=False,
                     level_ticks=7, ax=ax6, fig=fig,
                     plot_pvalues=True, pvalue_data=giss_data["nao_sig_temp_OP"], bottom_labels=True, 
                     title=" (f) cor(NAO-Temperature) (OP): GISS", plot_projection=projection,)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.10, wspace=0.08)
    plt.savefig(os.path.join(path_plots, "correlation_temp.png"), format= "png", bbox_inches="tight", dpi=600)
    
    
plot_d18O()
plot_temp()
plot_prec()