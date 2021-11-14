#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 19:27:52 2021

@author: dboateng
"""

#importing modules 
import os
import math 
import xarray as xr 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

try:
    from .Climat_data import read_Gtopo
except:
    from Climat_data import read_Gtopo


def extract_alps(path, tile_name, minelev, extract_var=None, east_alps=None, west_alps=None, west_buffer=None, central_buffer=None,
                 east_buffer=None):
    """
    

    Parameters
    ----------
    path : TYPE: str
        DESCRIPTION. The directrory to all the tile files 
    tile_name : TYPE: str
        DESCRIPTION. : hich tile to use for modification (check the image in the original files folder)
    minelev : TYPE: float
        DESCRIPTION. Minimum elevation to extaract the Alps
    extract_var : TYPE, optional or yes
        DESCRIPTION. The default is None or To extract only the values to datarray
    east_alps : TYPE, optional or yes
        DESCRIPTION. The default is None. to mask out the eastern part
    west_alps : TYPE, optional or yes 
        DESCRIPTION. The default is None. Or mask ou the west part of the Alps 
    west_buffer : TYPE, optional
        DESCRIPTION. The default is None. Or yes to apply buffer gradient 
    central_buffer : TYPE, optional
        DESCRIPTION. The default is None. or yes to apply buffer gradient at the central Alps
    east_buffer : TYPE, optional
        DESCRIPTION. The default is None. or to apply buffer gradient at the eastern side after modification

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    # read dataset and var
    Dataset = read_Gtopo(path=path, tile_name=tile_name, extract_var=None) 
    
    if east_alps is not None:
        if east_alps in ["Yes", "yes", "yep"]:
            maxlat, minlat, maxlon, minlon = 48, 44, 16, 12
            #convert lon to -180 to 180
            #data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})
    
            lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
            lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
            data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z >= minelev))
            # add as mask and extract values
            
            Dataset.coords["east_alps_mask"] = (("lat", "lon"), data_extract.z)
            
            return Dataset
        else:
            print("define east_alps as yes")
    elif west_alps is not None:
        if west_alps in ["Yes", "yes", "yep"]:
            maxlat, minlat, maxlon, minlon = 48, 43, 10, 5
            #convert lon to -180 to 180
            #data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})
    
            lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
            lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
            data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z >= minelev))
            # add as mask and extract values
            
            Dataset.coords["west_alps_mask"] = (("lat", "lon"), data_extract.z)
            
            return Dataset
            
        else:
            print("define the west_alps as yes")
            
    elif east_buffer is not None:
        if east_buffer in ["Yes", "yes", "yep"]:
            maxlat, minlat, maxlon, minlon = 48, 44, 20, 15.5
            #convert lon to -180 to 180
            #data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})
    
            lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
            lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
            data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z >= minelev))
            # add as mask and extract values
            
            Dataset.coords["east_buffer_mask"] = (("lat", "lon"), data_extract.z)
            
            return Dataset
        else:
            None
            
    elif west_buffer is not None:
        if west_buffer in ["Yes", "yes", "yep"]:
            maxlat, minlat, maxlon, minlon = 48, 43, 5.5, 2
            #convert lon to -180 to 180
            #data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})
    
            lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
            lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
            data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z >= minelev))
            # add as mask and extract values
            
            Dataset.coords["west_buffer_mask"] = (("lat", "lon"), data_extract.z)
            
            return Dataset
        else:
            None
            
    elif central_buffer is not None:
        if west_buffer in ["Yes", "yes", "yep"]:
            maxlat, minlat, maxlon, minlon = 48, 43, 12.5, 9.5
            #convert lon to -180 to 180
            #data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})
    
            lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
            lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
            data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z >= minelev))
            # add as mask and extract values
            
            Dataset.coords["central_buffer_mask"] = (("lat", "lon"), data_extract.z)
            
            return Dataset
    else:
        maxlat, minlat, maxlon, minlon = 48, 43, 16, 2
        #convert lon to -180 to 180
        #data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})
    
        lat_range = (Dataset.lat >= minlat) & (Dataset.lat <= maxlat)
        lon_range = (Dataset.lon >= minlon) & (Dataset.lon <= maxlon)
        data_extract = np.ones((Dataset.dims["lat"], Dataset.dims["lon"])) * Dataset.where((lat_range & lon_range) & (Dataset.z >= minelev))
        # add as mask and extract values
            
        Dataset.coords["alps_mask"] = (("lat", "lon"), data_extract.z)
            
        return Dataset
   

def modify_Gtopo(factor, part_of_alps, path_to_store, path, tile_name, minelev, extract_var=None, east_alps=None, west_alps=None, west_buffer=None, central_buffer=None,
                 east_buffer=None):
    """
    

    Parameters
    ----------
    factor : TYPE: float
        DESCRIPTION. The factor to multiply for changing the original topography
    part_of_alps : TYPE: str
        DESCRIPTION. Eg. full or whole, East_alps, West_alps 
    path_to_store : TYPE:str
        DESCRIPTION. The directory to save the output
    path : TYPE: str
        DESCRIPTION. The directory to load the original paths
   tile_name : TYPE: str
        DESCRIPTION. : hich tile to use for modification (check the image in the original files folder)
    minelev : TYPE: float
        DESCRIPTION. Minimum elevation to extaract the Alps
    extract_var : TYPE, optional or yes
        DESCRIPTION. The default is None or To extract only the values to datarray
    east_alps : TYPE, optional or yes
        DESCRIPTION. The default is None. to mask out the eastern part
    west_alps : TYPE, optional or yes 
        DESCRIPTION. The default is None. Or mask ou the west part of the Alps 
    west_buffer : TYPE, optional
        DESCRIPTION. The default is None. Or yes to apply buffer gradient 
    central_buffer : TYPE, optional
        DESCRIPTION. The default is None. or yes to apply buffer gradient at the central Alps
    east_buffer : TYPE, optional
        DESCRIPTION. The default is None. or to apply buffer gradient at the eastern side after modification

    Returns
    -------
    Dataset : TYPE: Dataset
        DESCRIPTION. Returns Dataset with similar headers like the original file from Gtopo

    """
    
    if part_of_alps in ["full", "whole"]:
        Dataset = extract_alps(path=path, tile_name=tile_name, minelev=800)
        if factor == 0:
            Dataset["z_modified"] = xr.where(np.isnan(Dataset["alps_mask"])==False, 250, Dataset["z"])
        else:
            Dataset["z_modified"] = xr.where(np.isnan(Dataset["alps_mask"])==False, Dataset["z"] * factor, Dataset["z"])
        
        #drop the original data and masks
        Dataset = Dataset.drop_vars(["z", "alps_mask"])
        
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
        Dataset.to_netcdf(os.path.join(path_to_store, "Alps_modified_"+ part_of_alps + "_" + str(factor) + ".nc"), format = "NETCDF3_CLASSIC")
        
        return Dataset
    elif part_of_alps in ["East_alps", "east_alps", "East_Alps"]:
        Dataset = extract_alps(path=path, tile_name=tile_name, minelev=800, east_alps="yes")
        
        print("Not yet implemented")
    elif part_of_alps in ["West_alps", "west_alps", "West_Alps"]:
        print("No yet implemented")
    else:
        None 
            
def visualise_topo(data, path_to_store, name_of_plot):
    """
    The function to visualise the modified topo file 

    Parameters
    ----------
    data : TYPE: Dataset
        DESCRIPTION. The Dataset containing the modified topography
    path_to_store : TYPE: str
        DESCRIPTION. The Dataset to store the figure
    name_of_plot : TYPE:str
        DESCRIPTION. The name of the file

    Returns
    -------
    None.

    """
    projection = ccrs.PlateCarree()
    fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":projection})
    p = data["z"].plot.imshow(ax =ax, cmap=plt.cm.terrain, vmin=0, vmax=10000, levels=25, transform = ccrs.PlateCarree(), 
                                      norm=None, cbar_kwargs= {"drawedges": True, "orientation": "vertical", "shrink": 0.40},
                                      extend= "neither")
    p.colorbar.set_label(label="Elevation" + " [m]", size = 20)  
    p.colorbar.ax.tick_params(labelsize = 20)
    p.axes.set_global()                    # setting global axis 
    p.axes.coastlines(resolution = "50m")  # add coastlines outlines to the current axis
    p.axes.add_feature(cfeature.BORDERS, color="black", linewidth = 0.5)
    p.axes.set_extent([-15, 20, 40, 65], ccrs.PlateCarree())
    gl= p.axes.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 0.5,
                          color = "gray", linestyle = "--")
        
    gl.top_labels = True                  # labesl at top
    gl.right_labels = True
    gl.xformatter = LONGITUDE_FORMATTER     # axis formatter
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": 15, "color": "black", "weight": "bold"}   #axis style 
    gl.ylabel_style = {"color": "black", "size": 15, "weight": "bold"}
    
    plt.tight_layout()
    plt.savefig(os.path.join(path_to_store, name_of_plot), format="tiff", bbox_inches="tight")
    
        
    

# #### control script ####
# path = "/home/dboateng/Gtopo30/001_original_nc/"
# path_to_store = "/home/dboateng/Gtopo30/Modified_topo/"
# tile_name = "W020N90.nc"

# data = extract_alps(path=path, tile_name=tile_name, minelev=800, east_alps=None, west_alps=None)

# data_modi = modify_Gtopo(factor=0, part_of_alps="full", path_to_store=path_to_store, path=path, tile_name=tile_name, minelev=800, extract_var=None, east_alps=None, west_alps=None, west_buffer=None, central_buffer=None,
#                  east_buffer=None)

# visualise_topo(data_modi, path_to_store,"Alps_0_full.tiff")