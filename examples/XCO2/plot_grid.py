# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 19:30:47 2023

@author: dboateng
"""
import xarray as xr
import os
import glob
import  pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator, LongitudeLocator)
from cartopy.util import add_cyclic_point
import calendar
import matplotlib.colors as mcolors

def plot_background(p, domain=None, use_AlbersEqualArea=False,ax=None, left_labels=True,
                    bottom_labels=True, plot_borders=False, coast_resolution=None):
  
    p.axes.set_global()                    # setting global axis 
    

    if coast_resolution is None:
        p.axes.coastlines(resolution = "50m")  # add coastlines outlines to the current axis (110m, 50m , 10m)
        
    else:
        p.axes.coastlines(resolution = coast_resolution, linewidth=1.5, color="black")  # add coastlines outlines to the current axis
    
    if plot_borders == True:
        p.axes.add_feature(cfeature.BORDERS, edgecolor="black", linewidth = 0.3) #adding country boarder lines
    
    #setting domain size
    if domain is not None: 
        if domain == "Europe":   # Europe
            minLon = -20
            maxLon = 35
            minLat = 34
            maxLat = 65
        elif domain == "South America":   # South America
            minLon = -83
            maxLon = -32
            minLat = -35
            maxLat = 5
      
        
        elif domain == "Africa":  # Africa
            minLon = -30
            maxLon = 55
            minLat = -35
            maxLat = 40
       
            
        elif domain == "West Africa":
            minLon = -25
            maxLon = 40
            minLat = -5
            maxLat = 35
        
            
        else:
            print("ERROR: invalid geographical domain passed in options")
        p.axes.set_extent([minLon, maxLon, minLat, maxLat], ccrs.PlateCarree())
    if domain is None: 
        p.axes.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        
    
    # adding gridlines    
    gl= p.axes.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1,
                 edgecolor = "gray", linestyle = "--", color="gray", alpha=0.5)
    
    gl.top_labels = False                  # labesl at top
    gl.right_labels = False
    
    if left_labels == True:
        gl.left_labels = True
    else:
        gl.left_labels = False
    
    if bottom_labels == True:
        gl.bottom_labels =True
    else:
        gl.bottom_labels = False
        
    gl.xformatter = LongitudeFormatter()     # axis formatter
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}   #axis style 
    gl.ylabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}


def plot_XCO2(variable, data, cmap, units, ax=None, vmax=None, vmin=None, levels=None, domain=None, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None,  left_labels= True, bottom_labels=True, add_colorbar=True, 
                     fig=None, cbar_pos=None, use_colorbar_default=False, 
                     orientation = "horizontal", time=None, plot_projection=None, coast_resolution=None, plot_borders=False,):
  
   
    projection = ccrs.PlateCarree()
    
    if plot_projection is None:
        plot_projection = ccrs.PlateCarree()
        
    #generating plot using geoaxis predefined or from default
    if ax is None:
        fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":plot_projection})
        
    if add_colorbar == True:    
        if cbar_pos is None:
            cbar_pos = [0.90, 0.30, 0.03, 0.40]
        
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
        
    
    
    if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
        ticks = np.linspace(vmin, vmax, level_ticks)
        
        if add_colorbar == True:
            
            if use_colorbar_default == True:
                p = data.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                     levels=levels, transform = projection, 
                                     cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": orientation, 
                                                   "shrink": 0.70, "ticks":ticks}, extend= "both",
                                     add_colorbar=True, add_labels=False)
            else:
                
                p = data.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                     levels=levels, transform = projection, 
                                     cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": orientation, 
                                                   "shrink": 0.70, "ticks":ticks}, extend= "both",
                                     add_colorbar=True, cbar_ax = cbar_ax, add_labels=False)
        else:
            p = data.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                 levels=levels, transform = projection, add_colorbar=False, add_labels=False)
                                     
    # when limits are not defined for the plot            
    else:
        p = data.plot.imshow(ax =ax, cmap=cmap, transform = projection, 
                                 cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": orientation, 
                                               "shrink": 0.70,}, extend= "neither", 
                                 add_colorbar=True, cbar_ax = cbar_ax, add_labels=False)
    
    
    if add_colorbar == True:
        
        p.colorbar.set_label(label=variable + " [" + units + "]", size= 20, fontweight="bold")
        p.colorbar.ax.tick_params(labelsize=20, size=0)
       
    
    # ploting background extent
    plot_background(p, domain= domain, left_labels=left_labels, bottom_labels=bottom_labels,
                    coast_resolution=coast_resolution, plot_borders=plot_borders)
    
    if title is not None:
        ax.set_title(title, fontsize=20, weight="bold", loc="left")
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
   
def count_files_in_directory(path_to_data, glob_pattern="*.nc"):
    files = glob.glob(path_to_data + "/" + glob_pattern)
    return files        



main_path_to_data_tccon = "D:/Datasets/OCO2/TCCON/"

files = count_files_in_directory(main_path_to_data_tccon)
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/scratch/OCO2_OCO3/OCO2_OCO3_grid_1"

ds = xr.open_dataset("D:/Datasets/OCO3/monthly_gridded/raw_OCO2_OCO3_gridded_monthly_1.0_grid.nc")

# path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/scratch/OCO2"

# ds = xr.open_dataset("D:/Datasets/OCO3/monthly_gridded/raw_OCO2_gridded_monthly_1.0_grid.nc")

for i in np.arange(ds.time.shape[0]):
    data = ds.XCO2[i]
    vmax=430
    vmin=400
    levels=30
    
    
    from pyClimat.plot_utils import apply_style
    apply_style(fontsize=25, style="seaborn-talk", linewidth=4,)
    # plotting
    plot_projection = ccrs.Robinson(central_longitude=0, globe=None)
    projection = ccrs.PlateCarree()
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18,15), subplot_kw={"projection":plot_projection})
    
    cmap_levels =np.linspace(vmin, vmax, levels)
    cmap_adj = plt.cm.get_cmap("jet", len(cmap_levels)-1)
    norm_adj = mcolors.BoundaryNorm(boundaries=cmap_levels, ncolors=cmap_adj.N) 
    
    title = str(ds.time.data[i])[:10]
    
    plot_XCO2(variable="XCO2", data=data, cmap=cmap_adj, units="ppm", ax=ax, vmax=vmax, vmin=vmin, levels=levels, domain=None, output_name=None, 
                         output_format=None, level_ticks=9, title=title, path_to_store=None,  left_labels= True, bottom_labels=True, add_colorbar=True, 
                         fig=fig, cbar_pos=None, use_colorbar_default=False, 
                         orientation = "vertical", time=None, plot_projection=None, coast_resolution="50m", plot_borders=True,)
    
    for file in files:
        tccon_data  = xr.open_dataset(os.path.join(main_path_to_data_tccon, file))
        
        lat = tccon_data["lat"][0].data
        lon = tccon_data["long"][0].data
        
        
        df_xco2 = pd.DataFrame({
                                'Xco2': tccon_data['xco2'],
                                "time": pd.to_datetime(tccon_data["time"]),
                                "year": tccon_data["year"],
                                })
        
        df = df_xco2.set_index("time")
        
        #compute monthly
        df_month = df.resample("MS").mean()
        
        if ds.time.data[i] in df_month.index.values:
            mon_value = df_month.loc[ds.time.data[i]]["Xco2"]
            
            #ax.annotate(tccon_data.short_location, (lon, lat), transform=plot_projection)
            ax.scatter(x=lon, y=lat, c=mon_value, cmap=cmap_adj, norm=norm_adj, edgecolor="red", s= 140,
                       transform=projection, linewidth=2)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
    plt.savefig(os.path.join(path_to_plots, title + ".png"), format= "png", bbox_inches="tight", dpi=600)