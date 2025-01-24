# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 13:57:33 2025

@author: dboateng
"""
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
import matplotlib.patches as patches


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_echam_topo
from pyClimat.analysis import extract_profile

path = "D:/Datasets/topo/modified_topo_PD/global_gtopo10_smoothW1E2.nc"
path_to_plots = "C:/Users/dboateng/Desktop"


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


data = xr.open_dataset(path,)
data = data.topo10_smth


apply_style(fontsize=28, style=None, linewidth=2.5) 
            
projection = ccrs.Robinson(central_longitude=0, globe=None)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(13,12), subplot_kw={"projection":projection})

norm, terrain = creat_norm()

plot_echam_topo(variable="Elevation", data=data, ax=ax, cmap=terrain, units="m", vmax=4000, vmin=-100, 
                levels=31, level_ticks=6, cbar=True, cbar_position= [0.35, 0.05, 0.45, 0.02], 
                cbar_orientation="horizontal", norm=norm, plot_coastlines=True, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=True, domain="Alps")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95)
plt.savefig(os.path.join(path_to_plots, "topoW1E2.png"), format= "png", bbox_inches="tight", dpi=600)

plt.show()