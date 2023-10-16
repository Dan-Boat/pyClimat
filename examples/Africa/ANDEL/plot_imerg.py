# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:22:21 2023

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

from pyClimat.plots import plot_annual_mean
from pyClimat.plot_utils import *


path_to_data = "D:/Datasets/IMERG/"
path_to_store = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/ANDEL/"


datasets = xr.open_mfdataset(path_to_data+"*.nc4")

# select the prec
prec = datasets.precipitationCal.compute()
prec = prec.transpose("time", "lat", "lon")

for i in range(len(prec.time)):
    data = prec[i]
    
    projection = ccrs.PlateCarree()
    plot_annual_mean(variable="PrecCal", data_alt=data, cmap=PrecipitationCal, units="mm/hr",
                      vmax=40, vmin=1, levels=20, level_ticks=9, orientation="vertical", 
                      domain="West Africa ANDEL", center=False, title = str(data.time.values)[:19],
                      plot_borders=True, coast_resolution="50m", output_name=str(i),
                      output_format="png", path_to_store=path_to_store)


from PIL import Image
import os
import glob

# Directory containing PNG files
output_gif = "output.gif"

# Get a list of all PNG files in the directory
png_files = glob.glob(os.path.join(path_to_store, "*.png"))

# Sort the list of PNG files by filename, if necessary
png_files.sort()

# Create a list to store the image objects
images = []

# Loop through the PNG files and append them to the images list
for png_file in png_files:
    img = Image.open(png_file)
    images.append(img)

# Save the list of images as a GIF
images[0].save(os.path.join(path_to_store,output_gif), save_all=True, append_images=images[1:], loop=0, duration=100)

#print(f"GIF saved as {output_gif}")
