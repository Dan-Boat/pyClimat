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

from pyClimat import *


# import paths
path_to_mmg = "/home/dboateng/Model_output_pst/Miotopofiles/mmg_topo_bathy_2160x1080.nc"
path_to_mmco = "/home/dboateng/Model_output_pst/Miotopofiles/CTL_Mio450/T159_MIO_450ppm_jan_surf_Frigola.nc"


# read data before editing (plotting)
mmg_data = xr.open_dataset(path_to_mmg)
mmco_data = xr.open_dataset(path_to_mmco)

mmco_topo = mmco_data.OROMEA
mmg_topo = mmg_data.topo

# try to plot with pyClimat

levels = [i for i in range(-100, 3000, 100)]
terrain_new = mpl.cm.get_cmap("terrain", 256)
terrain_adjust = ListedColormap(terrain_new(np.linspace(0.23, 1, 256)))
new_colors = terrain_adjust(np.linspace(0,1,256))

blue = np.array([135/256, 206/256, 250/256, 1])
new_colors[:1, :] = blue
terrain_shift = ListedColormap(new_colors)

norm_new = col.BoundaryNorm(levels, ncolors=terrain_shift.N, clip=True)
projection = ccrs.EuroPP()
plt.subplots(nrows = 1, ncols = 1, figsize=(20, 15), subplot_kw={"projection":projection})

plot_echam_topo(variable="Elevation", data=mmco_topo, cmap=terrain_shift, units="[m]", 
                vmax=2800, vmin=-100, levels=31, level_ticks=6,
                domain="Europe", cbar=True, cbar_position= [0.90, 0.30, 0.02, 0.40], cbar_orientation="vertical",
                projection=projection, norm=norm_new)
