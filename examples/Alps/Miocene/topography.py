# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:22:54 2023

@author: dboateng

This script plots the topography configuration used for the Miocene experiments 
(no need to plot the first, just reference)
Note that, the coastlines are plotted with confour of the land-sea mask!
"""

# import models
import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as col
import matplotlib as mpl 
from cartopy.util import add_cyclic_point
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_echam_topo


# set paths to data and store 
path_to_data = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/jan_surf_files"
W1E1_filename = "T159_MIO_W1E1_jan_surf_Herold.nc"
W2E1_filename = "T159_MIO_W2E1_jan_surf_Herold.nc"
W2E0_filename = "T159_MIO_W2E0_jan_surf_Herold.nc"
W2E15_filename = "T159_MIO_W2E1.5_jan_surf_Herold.nc"
W2E2_filename = "T159_MIO_W2E2_jan_surf_Herold.nc"
CTL_filename = "T159_PI_W1E1_jan_surf.nc"


def read_jan_surf_oromea(path, filename, return_slm=False):
    """
    Generate the orography with -100 for the ocean surface. This help is creating +
    colormap for the topography

    Parameters
    ----------
    path : TYPE: STR
        DESCRIPTION. Dir containing the data
    filename : TYPE: STR
        DESCRIPTION. The name of the file

    Returns
    -------
    oromea : TYPE: DATARRAY
        DESCRIPTION.

    """
    # plot jan_surf 
    jan_surf = xr.open_dataset(os.path.join(path, filename))   
    oromea = jan_surf.OROMEA
    
    oromea = xr.where(jan_surf.SLM == 0, -100, oromea)
    
    if return_slm == True:
        slm = jan_surf.SLM
        
        return oromea, slm
    else:
        return oromea

def creat_norm():
    levels = [i for i in range(-100, 4000, 100)]
    terrain_new = mpl.cm.get_cmap("terrain", 256)
    terrain_adjust = ListedColormap(terrain_new(np.linspace(0.23, 1, 256)))
    new_colors = terrain_adjust(np.linspace(0,1,256))
    
    blue = np.array([135/256, 206/256, 250/256, 1])
    new_colors[:1, :] = blue
    terrain_shift = ListedColormap(new_colors)

    norm_new = col.BoundaryNorm(levels, ncolors=terrain_shift.N, clip=True)
    
    return norm_new, terrain_shift



CTL_topo, CTL_slm = read_jan_surf_oromea(path=path_to_data, filename=CTL_filename, return_slm=True)
W1E1_topo, Mio_slm = read_jan_surf_oromea(path=path_to_data, filename=W1E1_filename, return_slm=True)


projection = ccrs.EuroPP()
fig, ((ax1,ax2,ax3), (ax4, ax5,ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(26, 20), 
                                                    subplot_kw={"projection": projection})

norm, terrain = creat_norm()

plot_echam_topo(variable="Elevation", data=CTL_topo, cmap=terrain, units="m", 
                vmax=4000, vmin=-100, levels=31, level_ticks=6,
                domain="Europe", cbar=True, cbar_position= [0.35, 0.05, 0.25, 0.02], cbar_orientation="horizontal",
                projection=projection, norm=norm, plot_coastlines=True, plot_borders=False, ax=ax1, 
                title="CTL (PI)", bottom_labels=False, fig=fig)

plot_echam_topo(variable="Elevation", data=W1E1_topo, cmap=terrain, units="m", 
                vmax=4000, vmin=-100, levels=31, level_ticks=6,
                domain="Europe", cbar=False, projection=projection, norm=norm, plot_coastlines=False, plot_borders=False, 
                ax=ax2, title="W1E1 (MIO)", bottom_labels=False, left_labels=False)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.01, hspace=0.06)
plt.show()