# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:02:11 2023

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

from pyClimat.plot_utils import cut_terrain_map
from pyClimat.plots import plot_echam_topo



path_to_ctl = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/orography.nc"
path_to_store = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/"
gtopo = "D:/Datasets/ECHAM5/Inputs/CTL_gtopo/global_gtopo10.nc"


# read the modern topo
gtopo_data = xr.open_dataset(gtopo)
mio_data = xr.open_dataset(path_to_ctl)

gtopo_data_interp = gtopo_data.interp(lat=mio_data.lat).interp(lon=mio_data.lon)





# Define the subsection with minimum and maximum longitude and latitude values
min_lon = 26
max_lon = 60
min_lat = 0
max_lat = 20


data = mio_data

lat_range = (data.lat >= min_lat) & (data.lat <= max_lat)
lon_range = (data.lon >= min_lon) & (data.lon <= max_lon)
# mask the whole africa


data_east_africa_grid = np.ones((data.dims["lat"], data.dims["lon"])) * data.where((lat_range & lon_range)
                                                                              & (data.topo >= 400))

data.coords["east_africa"] = (("lat", "lon"), data_east_africa_grid.topo.data)

data["topo_modify"] = xr.where(np.isnan(data["east_africa"])==False, 
                                       data["topo"] * 3 , data["topo"])



# Define the subsection with minimum and maximum longitude and latitude values
min_lon = 2
max_lon = 60
min_lat = -35
max_lat = 2


lat_range = (data.lat >= min_lat) & (data.lat <= max_lat)
lon_range = (data.lon >= min_lon) & (data.lon <= max_lon)
# mask the whole africa


data_south_africa_grid = np.ones((data.dims["lat"], data.dims["lon"])) * data.where((lat_range & lon_range)
                                                                              & (data.topo >= 450))

data.coords["south_africa"] = (("lat", "lon"), data_south_africa_grid.topo.data)

data["topo_modify"] = xr.where(np.isnan(data["south_africa"])==False, 
                                       data["topo"] * 2.5, data["topo_modify"])



min_lon = 5
max_lon = 42
min_lat = -35
max_lat = 20


lat_range = (data.lat >= min_lat) & (data.lat <= max_lat)
lon_range = (data.lon >= min_lon) & (data.lon <= max_lon)
# mask the whole africa


data_africa_grid = np.ones((data.dims["lat"], data.dims["lon"])) * data.where((lat_range & lon_range)
                                                                              & (data.topo_modify >= 1400))

data.coords["africa"] = (("lat", "lon"), data_africa_grid.topo.data)

data["topo_modify"] = xr.where(np.isnan(data["africa"])==False, 
                                       data["topo_modify"] + 600, data["topo_modify"])


data = data.drop_vars(["topo", "east_africa", "south_africa", "africa"])
data = data.rename({"topo_modify": "topo"})

data.to_netcdf(os.path.join(path_to_store, "East_africa_high"), format = "NETCDF3_CLASSIC")



levels = [i for i in range(-100, 4000, 100)]
terrain_new = mpl.cm.get_cmap("terrain", 256)
terrain_adjust = ListedColormap(terrain_new(np.linspace(0.23, 1, 256)))
new_colors = terrain_adjust(np.linspace(0,1,256))

blue = np.array([135/256, 206/256, 250/256, 1])
new_colors[:1, :] = blue
terrain_shift = ListedColormap(new_colors)

norm_new = col.BoundaryNorm(levels, ncolors=terrain_shift.N, clip=True)

projection = ccrs.PlateCarree()


topo_data = data.topo
topo_data["lon"] = gtopo_data.lon.values
plot_echam_topo(variable="Elevation", data=topo_data, cmap=terrain_shift, units="m", 
                vmax=4000, vmin=-100, levels=31, level_ticks=6,
                domain="Africa", cbar=True, cbar_position= [0.35, 0.05, 0.25, 0.02], cbar_orientation="horizontal",
                projection=projection, norm=norm_new, plot_coastlines=True, plot_borders=False, 
                title="[A] PI-W1E1", bottom_labels=True)