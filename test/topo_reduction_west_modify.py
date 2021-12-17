#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 12:04:36 2021

@author: dboateng

This script is a test for writing a function for topo reduction in the West setup:
"""

# importing modules 
import xarray as xr 
import os 
import math 
import pandas as pd
import numpy as np 

# import from Package
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 
from Package import *

    
# defining paths and filename 
path_to_tile ="/home/dboateng/Datasets/Topo_files/"
tile_name = "W020N90.nc"
path_to_store = "/home/dboateng/Datasets/Topo_files/Modified_topo"

def haversine(lon1, lat1, lon2, lat2):
# convert decimal degrees to radians 
    lon1 = np.deg2rad(lon1)
    lon2 = np.deg2rad(lon2)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6371
    return c * r    #km 



# reading topo 

Dataset = read_Gtopo(path=path_to_tile, tile_name=tile_name,)

maxlat, minlat, maxlon, minlon = 48, 43, 10, 5
minelev = 1000
#convert lon to -180 to 180
#data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})

lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z >= minelev))
# add as mask and extract values

Dataset.coords["west_alps_mask"] = (("lat", "lon"), data_extract.z)

Dataset["z_modified"] = xr.where(np.isnan(Dataset["west_alps_mask"])==False, Dataset["z"] * 3, Dataset["z"])

# resolve part in southern Alps that was uplifted 
maxlat, minlat, maxlon, minlon = 44.8875, 43.4625, 12.0125, 9.0875
lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z_modified >= minelev))
# add as mask and extract values

Dataset.coords["west_alps_mask_fixed"] = (("lat", "lon"), data_extract.z_modified)

Dataset["z_modified"] = xr.where(np.isnan(Dataset["west_alps_mask_fixed"])==False, Dataset["z_modified"]/3, Dataset["z_modified"])

#trying gradient buffer
max_dist = 100 # 2 grids when interpolated

buffer_lat, buffer_lon = 46.8, 9.8
lon_limit = 16
distance = haversine(Dataset.lon[:], Dataset.lat[:], buffer_lon, buffer_lat)
Dataset["distance"] = (("lat","lon"), distance)

Dataset["z_modified"] = xr.where(((Dataset.distance <= max_dist) & (Dataset.lon < lon_limit)) & (Dataset.z_modified <= 4000), Dataset["z_modified"]
                                 + (1500*np.tan(0.015)*(max_dist - Dataset.distance)), Dataset["z_modified"])

Dataset["z_modified"] = xr.where(((Dataset.distance <= max_dist) & ((Dataset.lon > 9.5) & (Dataset.lon < 16))) & (Dataset.z_modified > 1000), Dataset["z_modified"]
                                 - (1500*np.tan(0.015)*(max_dist - Dataset.distance)), Dataset["z_modified"])

Dataset["z_modified"] = xr.where(((Dataset.distance <= max_dist) & ((Dataset.lon > 11) & (Dataset.lon < 16))) & (Dataset.z_modified > 1000), Dataset["z_modified"]
                                 - (500*np.tan(0.015)*(max_dist - Dataset.distance)), Dataset["z_modified"])


# for east 0
maxlat, minlat, maxlon, minlon = 47.7, 45, 17, 11
lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z_modified >= 1200))
# add as mask and extract values

Dataset.coords["east_alps_mask"] = (("lat", "lon"), data_extract.z_modified)

Dataset["z_modified"] = xr.where(np.isnan(Dataset["east_alps_mask"])==False, Dataset["z_modified"]*0.05, Dataset["z_modified"])

# # applying a buffer between 11 and 13
# maxlat, minlat, maxlon, minlon = 47.8125, 45.64583, 13, 11
# lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
# lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
# data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z_modified >= 4000))
# Dataset.coords["west_alps_mask_buffer"] = (("lat", "lon"), data_extract.z_modified)
# Dataset["z_modified"] = xr.where(np.isnan(Dataset["west_alps_mask_buffer"])==False, Dataset["z_modified"]*0.5, Dataset["z_modified"])



# data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z_modified <= 4000))
# Dataset.coords["west_alps_mask_buffer"] = (("lat", "lon"), data_extract.z_modified)
# Dataset["z_modified"] = xr.where(np.isnan(Dataset["west_alps_mask_buffer"])==False, Dataset["z_modified"]*1.5, Dataset["z_modified"])



#drop the original data and masks
Dataset = Dataset.drop_vars(["z", "west_alps_mask", "west_alps_mask_fixed", "distance", ])

#rename z_modified to z because of topo scripts (ncl) and add headers, attributes like original data
Dataset = Dataset.rename({"z_modified":"z"})
Dataset["z"].attrs["long_name"] = "m"
Dataset["z"].attrs["actual_range"] = np.array([Dataset.z.min(), Dataset.z.max()])
Dataset["lon"].attrs["long_name"] = "longitude"
Dataset["lon"].attrs["actual_range"] = np.array([-20., 20.])
Dataset["lon"].attrs["units"] = "degrees_east"
Dataset["lat"].attrs["long_name"] = "latitude"
Dataset["lat"].attrs["actual_range"] = np.array([40., 90.])
Dataset["lat"].attrs["units"] = "degrees_north"
Dataset = Dataset.astype(dtype = "float32", casting= "same_kind")

#saving to path (using netcdf3_classic because of old version of ncl scripts use for upscaling into T159)
Dataset.to_netcdf(os.path.join(path_to_store, "Alps_modified_.nc"), format = "NETCDF3_CLASSIC")

