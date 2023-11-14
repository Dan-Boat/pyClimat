# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:20:22 2023

@author: dboateng
"""


import os
import numpy as np
import pandas as pd 
import xarray as xr
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import cartopy.crs as ccrs 
import cartopy.feature as cfeature
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)

import dask.dataframe as dd
import dask.array as da

# path to datasets
main_path_to_data_oco2 = "D:/Datasets/OCO2/Data/csv_files"

main_path_to_data_oco3 = "D:/Datasets/OCO3/Data/csv_files"

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/scratch/OCO2_OCO3"

path_to_plots = "F:\scratch"


def combine_datasets(path1, path2, glob_name,):
    
    data1 = dd.read_csv(path1 + "/" + glob_name, parse_dates=["DateTime"])
    
    data2 = dd.read_csv(path2 + "/" + glob_name, parse_dates=["DateTime"])
    
    #df = dd.concat([data1, data2], axis=0, interleave_partitions=True)
    
    
    return data1.compute(), data2.compute()

def plot_glob_datasets_traj(data, x_name="Longitude", y_name="Latitude", z_name="Xco2", vmax=None, vmin=None, 
                       levels=None, plot_projection=None, ax=None,
                       add_colorbar=True,cmap=None, path_to_store=None, output_name=None,  
                       orientation="horizontal", use_colorbar_default=False, title=None,
                       cbar_pos=None, fig=None, norm=None, color=None):
    
    projection = ccrs.PlateCarree()
    
    if plot_projection is None:
        plot_projection = ccrs.PlateCarree()
        
    #generating plot using geoaxis predefined or from default
    if ax is None:
        fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":plot_projection})
        
    if add_colorbar == True:    
        if cbar_pos is None:
            cbar_pos = [0.90, 0.30, 0.03, 0.45]
        
        if use_colorbar_default == False:
            
            
            cbar_ax = fig.add_axes(cbar_pos)   # axis for subplot colorbar # left, bottom, width, height
            
            if orientation == "vertical":
                cbar_ax.get_xaxis().set_visible(False)
                cbar_ax.yaxis.set_ticks_position('right')
                cbar_ax.set_yticklabels([])
                cbar_ax.tick_params(size=0)
            else:
                cbar_ax.get_yaxis().set_visible(False)
                cbar_ax.xaxis.set_ticks_position('bottom')
                cbar_ax.set_xticklabels([])
                cbar_ax.tick_params(size=0)
                
       
    sc = ax.scatter(data[x_name], data[y_name], cmap=cmap, 
               c=color, norm=norm, transform = projection, s=1, )
    
    
    if add_colorbar == True:
    
        cbar = plt.colorbar(mappable= sc, ax=ax, cax=cbar_ax, label= "XCO2 (ppm)",
                            shrink=0.5, extend="both", pad=0.05)
        cbar.ax.tick_params(labelsize=10)
        
    
    ax.coastlines(resolution='50m', linewidth=1.5, color="black")
    ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth = 0.3)
    ax.set_global()
    gl= ax.gridlines(crs = projection, draw_labels = True, linewidth = 1,
                     edgecolor = "gray", linestyle = "--", color="gray", alpha=0.5)
    
    gl.xformatter = LongitudeFormatter()     # axis formatter
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}   #axis style 
    gl.ylabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}
    
    ax.set_title(title, fontsize=20, weight="bold", loc="left")



for i in np.arange(1, 13):
    if i <10:
        id_num = "0"+str(i)
    else:
        id_num = str(i)
    
    files_glob_oco2 = "_xco2_2020" + id_num + "*.csv"
    
    df1, df2 = combine_datasets(path1=main_path_to_data_oco3, path2=main_path_to_data_oco2, 
                          glob_name=files_glob_oco2)

    
    plot_projection = ccrs.Robinson(central_longitude=0, globe=None)
    projection = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18,15), subplot_kw={"projection":plot_projection})
    
    plot_glob_datasets_traj(data=df1, vmax=430, vmin=390, levels=20, plot_projection=None, ax=ax, fig=fig,
                            add_colorbar=False,cmap=None, path_to_store=None, color="black", 
                            output_name=None, norm=None,
                            orientation="vertical", use_colorbar_default=False, title="OCO2 (red), OCO3(black) month="
                            + id_num + "(2020)",)
    
    plot_glob_datasets_traj(data=df2, vmax=430, vmin=390, levels=20, plot_projection=None, ax=ax, fig=fig,
                            add_colorbar=False,cmap=None, path_to_store=None, color="red", 
                            output_name=None, norm=None,
                            orientation="vertical", use_colorbar_default=False, title="OCO2 (red), OCO3(black) month="
                            + id_num + "(2020)",)
                            
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, id_num + ".png"), format= "png", bbox_inches="tight", dpi=600)