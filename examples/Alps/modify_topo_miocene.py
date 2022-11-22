# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 15:07:49 2022

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


# path to the Herold original files

path_to_ctl = "/home/dboateng/Model_output_pst/Miotopofiles/CTL_Mio_Herold/orography.nc"
path_to_ctl_t159 = "/home/dboateng/Model_output_pst/Miotopofiles/CTL_Mio_Herold/T159_MIO_W2E1_278ppm_jan_surf_Herold.nc"
path_to_store = "/home/dboateng/Model_output_pst/Miotopofiles/CTL_Mio_Herold/"


# load data 
def modify_topo_mio(filename):
    data = xr.open_dataset(path_to_ctl)
    
    
    
    #extract the Alps and its sections
    
    # add western Alps coords
    
    maxlat, minlat, maxlon, minlon = 48, 43, 8, 1
    minelev = 800
    #convert lon to -180 to 180
    #data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})
    
    lat_range = (data.lat >= minlat) & (data.lat <= maxlat)
    lon_range = (data.lon >= minlon) & (data.lon <= maxlon)
    data_extract = np.ones((data.dims["lat"], data.dims["lon"])) * data.where((lat_range & lon_range) & (data.topo >= minelev))
    # add as mask and extract values
    
    data.coords["west_alps_mask"] = (("lat", "lon"), data_extract.topo.data)
    
    data["topo_modify"] = xr.where(np.isnan(data["west_alps_mask"])==False, data["topo"] * 3, data["topo"])
    
    # inrease the eastern Alps a litle to avoid steps
    
    maxlat, minlat, maxlon, minlon = 48, 44, 15, 8
    minelev = 600
    #convert lon to -180 to 180
    #data_z = data_z.assign_coords({"lon": (((data_z.lon + 180) % 360) - 180)})
    
    lat_range = (data.lat >= minlat) & (data.lat <= maxlat)
    lon_range = (data.lon >= minlon) & (data.lon <= maxlon)
    data_extract = np.ones((data.dims["lat"], data.dims["lon"])) * data.where((lat_range & lon_range) & (data.topo >= minelev))
    # add as mask and extract values
    
    data.coords["east_alps_mask"] = (("lat", "lon"), data_extract.topo.data)
    
    if filename == "W2E1_orography.nc":
        
        data["topo_modify"] = xr.where(np.isnan(data["east_alps_mask"])==False, 
                                        data["topo_modify"] + data["topo_modify"] * 0.3 , data["topo_modify"])
        
    elif filename == "W2E0_orography.nc":
        data["topo_modify"] = xr.where(np.isnan(data["east_alps_mask"])==False, 
                                        250 , data["topo_modify"])
        
    elif filename == "W2E1.5_orography.nc":
        data["topo_modify"] = xr.where(np.isnan(data["east_alps_mask"])==False, 
                                        data["topo_modify"] + data["topo_modify"] * 0.9 , data["topo_modify"])
        
    elif filename == "W2E2_orography.nc":
        data["topo_modify"] = xr.where(np.isnan(data["east_alps_mask"])==False, 
                                        data["topo_modify"] * 3, data["topo_modify"])
    
    #save files 
    
    data = data.drop_vars(["topo", "west_alps_mask", "east_alps_mask"])
    data = data.rename({"topo_modify": "topo"})
    data.to_netcdf(os.path.join(path_to_store, filename), format = "NETCDF3_CLASSIC")

filenames = ["W2E1_orography.nc", "W2E1.5_orography.nc", "W2E0_orography.nc", "W2E2_orography.nc"]

for filename in filenames:
    print ("running for", filename)
    modify_topo_mio(filename)
    

def plot_jan_surf(path_to_ctl_t159):
    # plot jan_surf 
    ctl_jan_surf = xr.open_dataset(path_to_ctl_t159)   
    oromea = ctl_jan_surf.OROMEA
    oromea = xr.where(ctl_jan_surf.SLM == 0, -100, oromea)
       
    # # # add cyclic on the longitude coordinates
    # # lon = mmco_topo.coords["lon"]
    # # lon_index = mmco_topo.dims.index("lon")
    
    # # wrap_data, wrap_lon = add_cyclic_point(mmco_topo.values, coord=lon, axis=lon_index)
    
    levels = [i for i in range(-100, 3000, 100)]
    terrain_new = mpl.cm.get_cmap("terrain", 256)
    terrain_adjust = ListedColormap(terrain_new(np.linspace(0.23, 1, 256)))
    new_colors = terrain_adjust(np.linspace(0,1,256))
    
    blue = np.array([135/256, 206/256, 250/256, 1])
    new_colors[:1, :] = blue
    terrain_shift = ListedColormap(new_colors)
    
    norm_new = col.BoundaryNorm(levels, ncolors=terrain_shift.N, clip=True)
    projection = ccrs.EuroPP()
    # plt.subplots(nrows = 1, ncols = 1, figsize=(20, 15), subplot_kw={"projection":projection})
    
    plot_echam_topo(variable="Elevation", data=oromea, cmap=terrain_shift, units="[m]", 
                    vmax=3000, vmin=-100, levels=31, level_ticks=5,
                    domain="Europe", cbar=True, cbar_position= [0.90, 0.30, 0.02, 0.40], cbar_orientation="vertical",
                    projection=projection, norm=norm_new, plot_coastlines=True, plot_borders=False)
