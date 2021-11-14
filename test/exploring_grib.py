#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 14:19:07 2021

@author: dboateng
"""

import xarray as xr
path_to_data = "/home/dboateng/flexpart/input/era.grib"
path_to_era5_surf = "/home/dboateng/Desktop/GribtoArl/ERA5_2017_Aug22_surface.grib"
path_to_era5_sig = "/home/dboateng/Desktop/GribtoArl/ERA5_2017_Aug22_sigma_part.grib"

#opeaning data
#data = xr.open_dataset(path_to_data, engine= "cfgrib")

data_surf = xr.open_dataset(filename_or_obj = path_to_era5_surf, engine="cfgrib")
data_sig = xr.open_dataset(filename_or_obj= path_to_era5_sig, engine="cfgrib")


print(data_surf)