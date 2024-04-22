# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 12:31:39 2024

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl 
import cartopy.crs as ccrs
import matplotlib.colors as col
from matplotlib.colors import ListedColormap


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_echam_topo
from pyClimat.data import read_from_path, read_ERA_processed
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"


#path to data
path_to_data = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/jan_surf_files"
mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"


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

control_topo = read_jan_surf_oromea(path=path_to_data, filename="T159_MIO_W1E1_jan_surf_Herold.nc", 
                                    return_slm=False)

east_africa_high_topo = read_jan_surf_oromea(path=path_to_data, filename="T159_EA_high_MIO_278ppm_jan_surf_Herold.nc", 
                                    return_slm=False)

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# plotting
apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
            
projection = ccrs.Robinson(central_longitude=0, globe=None)

pprojection_trans = ccrs.PlateCarree()

fig,(ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16,13), subplot_kw={"projection":projection})


norm, terrain = creat_norm()
    
plot_echam_topo(variable="Elevation", data=control_topo, ax=ax1, cmap=terrain, units="m", vmax=4000, vmin=-100, 
                levels=31, level_ticks=6, cbar=True, cbar_position= [0.35, 0.05, 0.45, 0.02], 
                cbar_orientation="horizontal", norm=norm, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                projection=projection, title="(a) Miocene paleotopography (CTL)")


plot_echam_topo(variable="Elevation", data=east_africa_high_topo, ax=ax2, cmap=terrain, units="m", vmax=4000, vmin=-100, 
                levels=31, level_ticks=6, cbar=True, cbar_position= [0.35, 0.05, 0.45, 0.02], 
                cbar_orientation="horizontal", norm=norm, plot_coastlines=False, bottom_labels=False,
                left_labels=False, fig=fig, plot_borders=False, sea_land_mask=mio_slm,
                projection=projection, title="(b) Miocene paleotopography (High EARS)")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "poster_fig3.pdf"), format= "pdf", bbox_inches="tight", dpi=600)