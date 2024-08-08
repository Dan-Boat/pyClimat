# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 12:53:37 2024

@author: dboateng
"""

import xarray as xr 
import numpy as np
import matplotlib.pyplot as plt 
import os
import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs

from pyClimat.plots import plot_echam_topo
from pyClimat.plot_utils import terrain


def assign_coords(file: str):
    """
    Takes a netcdf file and assigns the coordinates to the file based on the x and y dimensions
    Returns the xarray dataset object with the coordinates assigned
    """
    ds = xr.open_dataset(file)
    # find the variable name from vars and assign it to the variable_name
    variable_name = [i for i in ds.data_vars][0]

    # check if the file has x and y dimensions and the bigger one is going to be used as longitude
    x_size = ds.dims['x']
    y_size = ds.dims['y']

    
    if x_size > y_size:
        # assign coords from 1 to 360 but NOT including 360 to the longitude
        ds = ds.assign_coords(lon = ("x", np.linspace(0, 360, x_size, endpoint=False)))
        interval = 180 / (y_size - 1)
        ds = ds.assign_coords(lat = ("y", np.linspace(-90 + interval/2, 90 - interval/2, y_size)))
        # Swap dimensions 'x' and 'y' with 'lat' and 'lon' for the specific variable
        ds[variable_name] = ds[variable_name].swap_dims({'x': 'lon', 'y': 'lat'})
    else:
        ds = ds.assign_coords(lon = ("x", np.linspace(0, 360, y_size, endpoint=False)))
        interval = 180 / (x_size - 1)
        ds = ds.assign_coords(lat = ("y", np.linspace(-90 + interval/2, 90 - interval/2, x_size)))
        # Swap dimensions 'x' and 'y' with 'lat' and 'lon' for the specific variable
        ds[variable_name] = ds[variable_name].swap_dims({'x': 'lat', 'y': 'lon'})

    # add attributes to the coordinates such as long_name and units
    ds['lon'].attrs['long_name'] = 'longitude'
    ds['lon'].attrs['units'] = 'degrees_east'
    ds['lat'].attrs['long_name'] = 'latitude'
    ds['lat'].attrs['units'] = 'degrees_north'
    # add attributes to the ccordinates such as standard_name and axis
    ds['lon'].attrs['standard_name'] = 'longitude'
    ds['lon'].attrs['axis'] = 'X'
    ds['lat'].attrs['standard_name'] = 'latitude'
    ds['lat'].attrs['axis'] = 'Y'

    return ds


path_to_topo = "D:/Datasets/topo/plasim_top.nc"

data = assign_coords(path_to_topo)

#data = data.assign_coords({"lon": (((data.lon + 180) % 360) - 180)})


maxlat, minlat, maxlon, minlon = 50, 0, 300, 275
minelev = 20000
#convert lon to -180 to 180
#data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})

lat_range = (data.lat >= minlat) & (data.lat <= maxlat)
lon_range = (data.lon >= minlon) & (data.lon <= maxlon)
data_extract = np.ones((data.dims["lat"], data.dims["lon"])) * data.where((lat_range & lon_range) & (data.var129 >= minelev))
# add as mask and extract values

data.coords["mask"] = (("time", "lat", "lon"), data_extract.var129.data)

data["var129_modify"] = xr.where(np.isnan(data["mask"])==False, data["var129"] * 2, data["var129"])

fig, (ax1,ax2) = plt.subplots(nrows = 2, ncols = 1, figsize=(25, 17))

data.var129_modify.isel(time=0).plot(ax=ax1, vmax=80000)

data.var129.isel(time=0).plot(ax=ax2, vmax=80000)



new_data = xr.open_dataset(path_to_topo)

new_data["var129"].values = data.var129_modify.values

new_data.to_netcdf("D:/Datasets/topo/plasim_top_modify.nc")

# cdo -f nc4 outputsrv plasim_top_modify.nc > N032_surf_0129new.sra

#set coordinates to change


