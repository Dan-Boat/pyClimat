# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 20:10:57 2023

@author: dboateng
1. compare correlation with echam d18Op and dD and NAO and EA index, 
extract regions to check distribution of correlation, and perform correlation with
regional means. Perferfom the analysis with ERA and GNIP datasets
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_var, extract_transect
from pyClimat.utils import extract_region


from paths_to_data import *
echam_pd_data_path = "D:/Datasets/Model_output_pst/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021"


def plot_correlation(variable, data, cmap = None, levels=None, units=None, ax=None, domain=None, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None, cbar = None, cbar_orientation=None, 
                     cbar_pos = None,plot_pvalues=False, pvalue_data=None,
                     fig=None, vmax=None, vmin=None, left_labels= True, bottom_labels=True,
                     ):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
    from cartopy.mpl.ticker import (LatitudeLocator, LongitudeLocator) 
    from cartopy.util import add_cyclic_point
    
    norm = MidpointNormalize(midpoint=0)
    projection = ccrs.PlateCarree()
    
    # defining axis
    if ax is None:
    
        fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":projection})
             
    if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
        ticks = np.linspace(vmin, vmax, level_ticks)
        if cbar==True:
            if cbar_pos is None:
                cbar_pos = [0.90, 0.30, 0.03, 0.45]
                
            
            cbar_ax = fig.add_axes(cbar_pos)   # axis for subplot colorbar # left, bottom, width, height
            if cbar_orientation == "vertical":
                cbar_ax.get_xaxis().set_visible(False)
                cbar_ax.yaxis.set_ticks_position('right')
            else:
                cbar_ax.get_yaxis().set_visible(False)
                cbar_ax.xaxis.set_ticks_position('bottom')
                
            cbar_ax.set_yticklabels([])
            cbar_ax.tick_params(size=0)
            
            p = data.plot.imshow(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, transform = projection,
                                 cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": cbar_orientation, 
                                               "shrink": 0.7, "format": "%.1f", "ticks":ticks}, extend= "both", 
                                 add_colorbar=True, cbar_ax=cbar_ax,
                                 add_labels=False)
            
            
            p.colorbar.set_label(label=variable + " [" + units + "]", size= 22, fontweight="bold")
            p.colorbar.ax.tick_params(labelsize=22, size=0,)
            
        elif cbar == False:
            p = data.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, transform = projection, 
                                 add_colorbar=False, add_labels=False)
                                 
    else:
        p = data.plot.imshow(ax =ax, cmap=cmap, transform = projection, 
                                 cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": cbar_orientation, 
                                               "shrink": 0.50, "format": "%.1f", "ticks":ticks}, extend= "both")
    
        
    

    plot_background(p, domain= domain, bottom_labels=bottom_labels, left_labels=left_labels)
    
    if plot_pvalues:
        if pvalue_data is not None:
            pvalue_data.plot.contourf(ax=ax, colors="none", hatches=["."], add_colorbar=False,
                                      add_labels=False)
    
    if title is not None:
        ax.set_title(title , fontsize=22, weight="bold", loc="left")
                

        
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
        
        
        
# read all the required datsets (for winter)
df_era_pcs = pd.read_csv(PCs_ERA_DJF_path , parse_dates=["time"])
df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

#select the node ands convert to xarray

nao_index_era = xr.DataArray(df_era_pcs[str(2)] * -1, dims="time", coords={"time": df_era_pcs["time"]})
nao_index_echam = xr.DataArray(df_echam_pcs[str(1)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})

ea_index_era = xr.DataArray(df_era_pcs[str(3)], dims="time", coords={"time": df_era_pcs["time"]})
ea_index_echam = xr.DataArray(df_echam_pcs[str(2)], dims="time", coords={"time": df_echam_pcs["time"]})



def perform_correlation(atmos_index, varname, units=None, sig=0.05, lev=None, lev_units=None):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data = extract_var(Dataset=echam_data, varname=varname, units=units, Dataset_wiso=echam_wiso,
                       lev=lev, lev_units=lev_units)
    
    data_season = extract_region(data=data, maxlon=40, minlon=-30, maxlat=80, minlat=30, time="season", 
                                 season="DJF")
    
    sval, pval, sig = StatCorr(x=data_season, y=atmos_index, dim="time",
                                                 return_sig=True, sig=0.05)
    
    return sval, pval, sig
    
d18op_nao_sval, d18op_nao_pval, d18op_nao_sig = perform_correlation(atmos_index=nao_index_echam, 
                                                                    varname="d18op", units="per mil",
                                                                    )
d18op_ea_sval, d18op_ea_pval, d18op_ea_sig = perform_correlation(atmos_index=ea_index_echam, 
                                                                    varname="d18op", units="per mil",
                                                                    )


t2m_nao_sval, t2m_nao_pval, t2m_nao_sig = perform_correlation(atmos_index=nao_index_echam, 
                                                                    varname="temp2", units="°C",
                                                                    )
t2m_ea_sval, t2m_ea_pval, t2m_ea_sig = perform_correlation(atmos_index=ea_index_echam, 
                                                                    varname="temp2", units="°C",
                                                                    )

prec_nao_sval, prec_nao_pval, prec_nao_sig = perform_correlation(atmos_index=nao_index_echam, 
                                                                    varname="prec", units="mm/month",
                                                                    )
prec_ea_sval, prec_ea_pval, prec_ea_sig = perform_correlation(atmos_index=ea_index_echam, 
                                                                    varname="prec", units="mm/month",
                                                                    )

     
apply_style(fontsize=22, style=None, linewidth=2)
projection = ccrs.PlateCarree()
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols=3, 
                                                 figsize=(28, 18), subplot_kw={"projection": projection})

plot_correlation(variable="Spearman Corr", data=d18op_nao_sval, units="-", vmax=1,
                 vmin=-1, cmap=Seismic, domain="Europe", levels=22,cbar=True, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax1, fig=fig, cbar_pos= [0.35, 0.01, 0.30, 0.02],
                 plot_pvalues=True, pvalue_data=d18op_nao_sig, bottom_labels=False, title="[A] NAO-$\delta^{18}$Op")

plot_correlation(variable="Spearman Corr", data=d18op_ea_sval, units="-", vmax=1,
                 vmin=-1, cmap=Seismic, domain="Europe", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax4, fig=fig, plot_pvalues=True, pvalue_data=d18op_ea_sig,
                 title="[D] EA-$\delta^{18}$Op")

plot_correlation(variable="Spearman Corr", data=t2m_nao_sval, units="-", vmax=1,
                 vmin=-1, cmap=Seismic, domain="Europe", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax2, fig=fig, plot_pvalues=True, pvalue_data=t2m_nao_sig, bottom_labels=False,
                 left_labels=False, title="[B] NAO-t2m")

plot_correlation(variable="Spearman Corr", data=t2m_ea_sval, units="-", vmax=1,
                 vmin=-1, cmap=Seismic, domain="Europe", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax5, fig=fig, plot_pvalues=True, pvalue_data=t2m_ea_sig, left_labels=False,
                 title="[E] EA-t2m")

plot_correlation(variable="Spearman Corr", data=prec_nao_sval, units="-", vmax=1,
                 vmin=-1, cmap=Seismic, domain="Europe", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax3, fig=fig, plot_pvalues=True, pvalue_data=prec_nao_sig, bottom_labels=False,
                 left_labels=False, title="[C] NAO-prec")

plot_correlation(variable="Spearman Corr", data=prec_ea_sval, units="-", vmax=1,
                 vmin=-1, cmap=Seismic, domain="Europe", levels=22,cbar=False, cbar_orientation="horizontal",
                 level_ticks=7, ax=ax6, fig=fig, plot_pvalues=True, pvalue_data=prec_ea_sig,
                 left_labels=False, title="[F] EA-prec")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, "corr_nao_ea_d18op.svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.show()
 