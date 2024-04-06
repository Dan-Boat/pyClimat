#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 18:53:11 2021

@author: dboateng
This routine contains funtional utilities required in the other modules
"""

#import xarray as xr
#import pandas as pd 
import numpy as np 
# import pickle 
# import os 
# from collections import OrderedDict 
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.colors as colors
# from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
#                                 LatitudeLocator)


# function to estimate distance 
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
        

def extract_indices_around(dataset, lat, lon, radius):
    close_grids = lambda lat_, lon_: haversine(lat, lon, lat_, lon_) <= radius
    
    
    if hasattr(dataset, "longitude"):
        LON, LAT = np.meshgrid(dataset.longitude, dataset.latitude)
    else:
        LON, LAT = np.meshgrid(dataset.lon, dataset.lat)
        
    grids_index = np.where(close_grids(LAT, LON))
    
    return grids_index




def vert_coord_convertion(data, units):
    """
    

    Parameters
    ----------
    data : TYPE: Dataarray (MxNxV)
        DESCRIPTION. Dataarray on vertical coordinates with Pa as default unit
    units : TYPE: str
        DESCRIPTION. The required unit for conversion eg. hPa

    Raises
    ------
    ValueError
        DESCRIPTION. When other units is defined aside hPa 

    Returns
    -------
    data : TYPE: Datarray 
        DESCRIPTION. Data vertical units converted

    """
    if hasattr(data, "mlev"):
        data = data.rename({"mlev":"lev"})
    elif hasattr(data, "plev"):
        data = data.rename({"plev":"lev"})
        
    if units =="hPa":
        data["lev"] = data.lev/100  #Pa -->hPa
        data["lev"].attrs["units"] = "hPa"
    else:
        raise ValueError("The units may be incorrect or not implemented yet!")
    
    return data


def extract_region(data, maxlon=60, minlon=-80, maxlat=80, minlat=10, time="season", 
                              season="DJF", month=None, regional_mean=False, select_time=False, 
                              time_range=None, use_slice=False, slice_start=None, slice_end=None):
    if hasattr(data, "longitude"):
        data = data.rename({"longitude":"lon", "latitude":"lat"})
        
    
    data = data.assign_coords({"lon": (((data.lon + 180) % 360) - 180)})

    data = data.where((data.lat >= minlat) & (data.lat <= maxlat), drop=True)
    data = data.where((data.lon >= minlon) & (data.lon <= maxlon), drop=True)
    
    if select_time:
        if use_slice:
            data = data.sel(time= slice(slice_start, slice_end))
        else:
            data = data.sel(time= time_range)
            
        
        
    if time =="season":
        data_group = data.groupby("time.season")
        
        if season is not None:
            data = data_group[season]
            
        else:
            data = data_group["DJF"]  # set to NH Winter as default
            
    elif time == "month":
        
        data_group = data.groupby("time.month")
        
        if isinstance(month, int):
            data = data_group.get(month)
        
        elif month == "ONDJFM":
            data = data.sel(time=data.time.dt.month.isin([1, 2, 3, 10, 11, 12]))
            
        elif month == "JJAS":
            data = data.sel(time=data.time.dt.month.isin([6, 7, 8, 9]))
            
        elif month == "AMJJAS":
            data = data.sel(time=data.time.dt.month.isin([4, 5, 6, 7, 8, 9]))
            
        elif month == "JA": # July and August for the SNAO
            data = data.sel(time=data.time.dt.month.isin([7, 8]))
            
        else:
            raise ValueError("The define month parameter is not recognized")
            
    if regional_mean:
        data = data.mean(dim=("lat", "lon"))
            
    return data   