# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 09:47:14 2023

@author: dboateng
"""
# import modules
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
main_path_to_data = "D:/Datasets/OCO2/Data/csv_files"
file = "_xco2_2015010100045007_.csv"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/scratch"

def count_files_in_directory(path_to_data, glob_pattern="*.nc"):
    files = glob.glob(path_to_data + "/" + glob_pattern)
    print(len(files))
    return files


def filter_grid_cells(grouped_df, min_data_points):
    if len(grouped_df) > min_data_points:
        return grouped_df['Xco2'].mean()
    else:
        return np.nan


# extract only monthly files (use glob)
def grid_monthly(main_path_to_csv, files_glob, grid_interval=1, min_data_points=20):
    
    data = dd.read_csv(main_path_to_csv + "/" + files_glob, parse_dates=["DateTime"])
    
    x_cell_edges = da.arange(-180, 181, grid_interval)
    y_cell_edges = da.arange(-90, 91, grid_interval)

    x_cell_labels = (x_cell_edges[:-1] + x_cell_edges[1:]) / 2
    y_cell_labels = (y_cell_edges[:-1] + y_cell_edges[1:]) / 2
    
    df = pd.DataFrame()
    df['x'] = pd.cut(data['Longitude'].values.compute(), bins=x_cell_edges.compute(), labels=x_cell_labels.compute())
    df['y'] = pd.cut(data['Latitude'].values.compute(), bins=y_cell_edges.compute(), labels=y_cell_labels.compute())
    df['Xco2'] = data['Xco2'].values.compute()
    
    grid = df.groupby(['x', 'y']).apply(filter_grid_cells, min_data_points=min_data_points).reset_index(name='mean').dropna()
    
    return grid 


# function for ploting
def plot_glob_datasets(data, vmax=None, vmin=None, levels=None, plot_projection=None, ax=None,
                       add_colorbar=True,cmap="cividis", path_to_store=None, output_name=None,  
                       orientation="horizontal", use_colorbar_default=False, title=None,
                       cbar_pos=None, fig=None, norm=None):
    
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
                
    #cmap_levels =np.linspace(data["mean"].min(), data["mean"].max(), levels)
    cmap_levels =np.linspace(vmin, vmax, levels)
    cmap_adj = plt.cm.get_cmap(cmap, len(cmap_levels)-1)
    norm_adj = mcolors.BoundaryNorm(boundaries=cmap_levels, ncolors=cmap_adj.N) 
           
    sc = ax.scatter(data["x"], data["y"], c = data["mean"], cmap=cmap, 
               norm=norm, transform = projection, s=10, )
    
    
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
    
    #plt.savefig(os.path.join(path_to_store, output_name + "." + "png"), format= output_format, bbox_inches="tight")

# function to read tccon file 
# def read_tccon_datasets(path_to_data, select_time=False, date_range=None):
    
#     data = xr.open_dataset(path_to_data)
#     df = pd.DataFrame({
#         'xco2': data['xco2'],
#         "time": data["time"],
#         "xco2_error": data["xco2_error"]})

#     return df        
    
    
    
if __name__=="__main__":
    # function that get file globe
    plot_projection = ccrs.Robinson(central_longitude=0, globe=None)
    projection = ccrs.PlateCarree()
    main_path_to_data_tccon = "D:/Datasets/OCO2/TCCON/"
    
    files = count_files_in_directory(main_path_to_data_tccon)
    
    vmax=406
    vmin=390
    levels=20
    
    #cmap_levels =np.linspace(data["mean"].min(), data["mean"].max(), levels)
    cmap_levels =np.linspace(vmin, vmax, levels)
    cmap_adj = plt.cm.get_cmap("jet", len(cmap_levels)-1)
    norm_adj = mcolors.BoundaryNorm(boundaries=cmap_levels, ncolors=cmap_adj.N) 
    
    
    for i in np.arange(1, 13):
        if i <10:
            id_num = "0"+str(i)
        else:
            id_num = str(i)
        files_glob = "_xco2_2015" + id_num + "*.csv"
        
        data = grid_monthly(main_path_to_csv=main_path_to_data, files_glob=files_glob,
                            grid_interval=2, min_data_points=20)
        
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18,15), subplot_kw={"projection":plot_projection})
        plot_glob_datasets(data, vmax=406, vmin=390, levels=20, plot_projection=None, ax=ax, fig=fig,
                                add_colorbar=True,cmap=cmap_adj, path_to_store=main_path_to_data, 
                                output_name=id_num, norm=norm_adj,
                                orientation="vertical", use_colorbar_default=False, title="XCO2 level 2 month=" + id_num)
        
        for file in files:
            tccon_data  = xr.open_dataset(os.path.join(main_path_to_data, file))
            
            lat = tccon_data["lat"][0].data
            lon = tccon_data["long"][0].data
        
        
            df_xco2 = pd.DataFrame({
                'Xco2': tccon_data['xco2'],
                "time": pd.to_datetime(tccon_data["time"]),
                "year": tccon_data["year"],
                })
        
            df = df_xco2.set_index("time")
        
            df_2015 = df.loc[df["year"]==2015]
            
            #df_2015["month"] = df_2015.index.month.values 
            
            
            if df_2015.empty == False:
                
                df_2015_month = df_2015.resample("M").mean()
                
                df_2015_month["month"] = df_2015_month.index.month.values 
                
                
                # check if the data has the whole year
                if i in df_2015_month.month.values:
                    mon_value = df_2015_month.loc[df_2015_month["month"]==i]["Xco2"].values
                    ax.scatter(x=lon, y=lat, c=mon_value, cmap=cmap_adj, norm=norm_adj, edgecolor="k", s= 140,
                               transform=projection)
                
                
            
            
        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout() 
        plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
        plt.savefig(os.path.join(path_to_plots, id_num + ".png"), format= "png", bbox_inches="tight", dpi=600)
            

        
        
        
        
    