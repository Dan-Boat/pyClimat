# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 13:45:02 2023

@author: dboateng
"""

# import os
# import numpy as np
# import pandas as pd
# import xarray as xr

# # Define the path to your CSV file
# path_to_csv = "C:/Users/dboateng/Desktop/Python_scripts/scratch/test.csv"

# # Create a DataFrame with latitude, longitude, and XCO2
# latitude = np.arange(-89.75, 90.25, 0.5)
# longitude = np.arange(-179.75, 180.25, 0.5)
# longitude, latitude = np.meshgrid(longitude, latitude)
# latitude = latitude.flatten()
# longitude = longitude.flatten()
# xco2 = [np.nan] * len(latitude)

# df1 = pd.DataFrame({'Latitude': latitude, 'Longitude': longitude, 'Xco2': xco2})

# # Load data from another CSV file
# df2 = pd.read_csv(path_to_csv, usecols=["Latitude", "Longitude", "Xco2"]).dropna()

# # Create a dictionary mapping (latitude, longitude) to new XCO2 values
# co2_mapping = dict(zip(zip(df2['Latitude'], df2['Longitude']), df2['Xco2']))

# # Use the map function to replace values in df1 based on the mapping
# df1['Xco2'] = df1.apply(lambda row: co2_mapping.get((row['Latitude'], row['Longitude']), row['Xco2']), axis=1)

# # Create an xarray dataset from the modified DataFrame
# ds = xr.Dataset({
#     'Xco2': ('latitude', 'longitude')}, df1['Xco2'].values.reshape(len(latitude), len(longitude))
#              , coords={'latitude': latitude, 'longitude': longitude})

# # Define the output NetCDF file path
# output_nc_file = "output_file.nc"

# # Save the xarray dataset to a NetCDF file
# ds.to_netcdf(output_nc_file)

# print(f"Saved the modified data to {output_nc_file}")


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

# Define the path to your CSV file
path_to_csv = "C:/Users/dboateng/Desktop/Python_scripts/scratch/test.csv"

# Create a DataFrame with latitude, longitude, and XCO2
latitude = np.arange(-89.75, 90.25, 0.5)
longitude = np.arange(-179.75, 180.25, 0.5)
longitude, latitude = np.meshgrid(longitude, latitude)
latitude = latitude.flatten()
longitude = longitude.flatten()
xco2 = [np.nan] * len(latitude)

df1 = pd.DataFrame({'Latitude': latitude, 'Longitude': longitude, 'Xco2': xco2})

# Load data from another CSV file
df2 = pd.read_csv(path_to_csv, usecols=["Latitude", "Longitude", "Xco2"]).dropna()

# Create a dictionary mapping (latitude, longitude) to new XCO2 values
co2_mapping = dict(zip(zip(df2['Latitude'], df2['Longitude']), df2['Xco2']))

# Update XCO2 values in df1 using the mapping
df1['Xco2'] = df1.apply(lambda row: co2_mapping.get((row['Latitude'], row['Longitude']), row['Xco2']), axis=1)


time = pd.date_range(start="2020-01-01", end="2020-12-31", freq="MS")


# Create an xarray dataset from the modified DataFrame
ds = xr.Dataset({
    'Xco2': (('latitude', 'longitude'), df1['Xco2'].values.reshape(360, 720))
}, coords={'latitude': np.arange(-89.75, 90.25, 0.5), 'longitude': np.arange(-179.75, 180.25, 0.5),
           'time': time[0]})
           
           
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

    
    p = data["Xco2"].plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                        levels=levels, transform = projection, norm=norm, 
                                        cbar_kwargs= {"pad":0.1, "drawedges": True, 
                                                      "shrink": 0.70}, extend= "both",
                                        add_colorbar=True, add_labels=False)
    
    p.colorbar.set_label(label="XCO2 (ppm)", size= 20, fontweight="bold")
    p.colorbar.ax.tick_params(labelsize=20, size=0,)
           
    # sc = ax.scatter(data[x_name], data[y_name], c = data[z_name], cmap=cmap, 
    #            norm=norm, transform = projection, s=10, )
    
    
    # if add_colorbar == True:
    
    #     cbar = plt.colorbar(mappable= sc, ax=ax, cax=cbar_ax, label= "XCO2 (ppm)",
    #                         shrink=0.5, extend="both", pad=0.05)
    #     cbar.ax.tick_params(labelsize=10)
        
    
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

    
# function that get file globe
plot_projection = ccrs.Robinson(central_longitude=0, globe=None)
projection = ccrs.PlateCarree()

path_to_plots = "F:\scratch"


vmax=420
vmin=400
levels=20
 
#cmap_levels =np.linspace(data["mean"].min(), data["mean"].max(), levels)
cmap_levels =np.linspace(vmin, vmax, levels)
cmap_adj = plt.cm.get_cmap("jet", len(cmap_levels)-1)
norm_adj = mcolors.BoundaryNorm(boundaries=cmap_levels, ncolors=cmap_adj.N) 

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18,15), subplot_kw={"projection":plot_projection})
plot_glob_datasets(ds, vmax=420, vmin=400, levels=20, plot_projection=None, ax=ax, fig=fig,
                add_colorbar=True,cmap=cmap_adj, path_to_store=path_to_plots, 
                norm=norm_adj, orientation="vertical", use_colorbar_default=False, title="XCO2 level 2 month=")

print(None)

# # Define the output NetCDF file path
# output_nc_file = "output_file.nc"

# # Save the xarray dataset to a NetCDF file
# ds.to_netcdf(output_nc_file)

# print(f"Saved the modified data to {output_nc_file}")
