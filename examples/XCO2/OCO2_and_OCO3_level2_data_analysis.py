# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 19:52:35 2023

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

def count_files_in_directory(path_to_data, glob_pattern="*.nc"):
    files = glob.glob(path_to_data + "/" + glob_pattern)
    print(len(files))
    return files


def combine_datasets(path1, path2, glob_name,):
    
    data1 = dd.read_csv(path1 + "/" + glob_name, parse_dates=["DateTime"])
    
    data2 = dd.read_csv(path2 + "/" + glob_name, parse_dates=["DateTime"])
    
    df = dd.concat([data1, data2], axis=0, interleave_partitions=True)
    
    
    return df.compute()

def filter_grid_cells(grouped_df, min_data_points):
    if len(grouped_df) > min_data_points:
        return grouped_df['Xco2'].mean()
    else:
        return np.nan


# extract only monthly files (use glob)
def grid_monthly(main_path_O2_to_csv, main_path_O3_to_csv, files_glob, grid_interval=1, min_data_points=20):
    
    data_O2 = dd.read_csv(main_path_O2_to_csv + "/" + files_glob, parse_dates=["DateTime"])
    data_O3 = dd.read_csv(main_path_O3_to_csv + "/" + files_glob, parse_dates=["DateTime"])
    
    data = dd.concat([data_O2, data_O3], axis=0, interleave_partitions=True)
    
    
    x_cell_edges = da.arange(-180, 181, grid_interval)
    y_cell_edges = da.arange(-90, 91, grid_interval)

    x_cell_labels = (x_cell_edges[:-1] + x_cell_edges[1:]) / 2
    y_cell_labels = (y_cell_edges[:-1] + y_cell_edges[1:]) / 2
    
    df = pd.DataFrame()
    df['Longitude'] = pd.cut(data['Longitude'].values.compute(), bins=x_cell_edges.compute(), labels=x_cell_labels.compute())
    df['Latitude'] = pd.cut(data['Latitude'].values.compute(), bins=y_cell_edges.compute(), labels=y_cell_labels.compute())
    df['Xco2'] = data['Xco2'].values.compute()
    
    grid = df.groupby(['Longitude', 'Latitude']).apply(filter_grid_cells, min_data_points=min_data_points).reset_index(name='Xco2').dropna()
    
    return grid

# function for ploting
def plot_glob_datasets(data, vmax=None, vmin=None, levels=None, plot_projection=None, ax=None,
                       add_colorbar=True,cmap="cividis", path_to_store=None, output_name=None,  
                       orientation="horizontal", use_colorbar_default=False, title=None,
                       cbar_pos=None, fig=None, norm=None, x_name="Longitude", y_name="Latitude", z_name="Xco2"):
    
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
           
    sc = ax.scatter(data[x_name], data[y_name], c = data[z_name], cmap=cmap, 
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

# function for ploting
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
    
def plot_trajectories():
    for i in np.arange(1, 13):
        if i <10:
            id_num = "0"+str(i)
        else:
            id_num = str(i)
        
        files_glob_oco2 = "_xco2_2020" + id_num + "*.csv"
        
        df = combine_datasets(path1=main_path_to_data_oco3, path2=main_path_to_data_oco2, 
                              glob_name=files_glob_oco2)
    
        
        plot_projection = ccrs.Robinson(central_longitude=0, globe=None)
        projection = ccrs.PlateCarree()
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18,15), subplot_kw={"projection":plot_projection})
        plot_glob_datasets_traj(data=df, vmax=406, vmin=390, levels=20, plot_projection=None, ax=ax, fig=fig,
                                add_colorbar=False,cmap=None, path_to_store=None, color="b", 
                                output_name=None, norm=None,
                                orientation="vertical", use_colorbar_default=False, title="OCO2, OCO3 level 2 trajectories month="
                                + id_num,)
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout() 
        plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
        plt.savefig(os.path.join(path_to_plots, id_num + ".png"), format= "png", bbox_inches="tight", dpi=600)
        
        
if __name__=="__main__":
    
    
    # function that get file globe
    plot_projection = ccrs.Robinson(central_longitude=0, globe=None)
    projection = ccrs.PlateCarree()
    main_path_to_data_tccon = "D:/Datasets/OCO2/TCCON/"
    
    files = count_files_in_directory(main_path_to_data_tccon)
    
    vmax=420
    vmin=400
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
        files_glob = "_xco2_2020" + id_num + "*.csv"
        
        data = grid_monthly(main_path_O2_to_csv=main_path_to_data_oco2, main_path_O3_to_csv=main_path_to_data_oco3,
                            files_glob=files_glob,
                            grid_interval=0.5, min_data_points=20)
        
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18,15), subplot_kw={"projection":plot_projection})
        plot_glob_datasets(data, vmax=406, vmin=390, levels=20, plot_projection=None, ax=ax, fig=fig,
                                add_colorbar=True,cmap=cmap_adj, path_to_store=path_to_plots, 
                                output_name=id_num, norm=norm_adj,
                                orientation="vertical", use_colorbar_default=False, title="XCO2 level 2 month=" + id_num)
        
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
        
            df_2015 = df.loc[df["year"]==2020]
            
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
        plt.savefig(os.path.join(path_to_plots, id_num + "OCO2_OC3_05_by_05_2020.png"), format= "png", bbox_inches="tight", dpi=600)
            
    
    
