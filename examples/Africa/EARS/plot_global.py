# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 20:40:07 2023

@author: dboateng
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
import matplotlib.patches as patches


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_echam_topo
from pyClimat.analysis import extract_profile


# set paths to data and store 
path_to_data = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/jan_surf_files"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/EARS"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"




W1E1_filename = "T159_EA_high_MIO_278ppm_jan_surf_Herold.nc"

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


W1E1_topo = read_jan_surf_oromea(path=path_to_data, filename=W1E1_filename, return_slm=False)

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM

def plot_global_topo(): 
    apply_style(fontsize=28, style=None, linewidth=2.5) 
            
    projection = ccrs.Robinson(central_longitude=0, globe=None)
    
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(13,12), subplot_kw={"projection":projection})
    
    norm, terrain = creat_norm()
    
    plot_echam_topo(variable="Elevation", data=W1E1_topo, ax=ax, cmap=terrain, units="m", vmax=4000, vmin=-100, 
                    levels=31, level_ticks=6, cbar=True, cbar_position= [0.35, 0.05, 0.45, 0.02], 
                    cbar_orientation="horizontal", norm=norm, plot_coastlines=False, bottom_labels=False,
                    left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                    projection=projection)
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95)
    plt.savefig(os.path.join(path_to_plots, "global_topo_with_high_EA.png"), format= "png", bbox_inches="tight", dpi=600)

    plt.show()
    
    
plot_global_topo()