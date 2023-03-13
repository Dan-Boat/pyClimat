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
import matplotlib.patches as patches


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_echam_topo
from pyClimat.analysis import extract_profile


# set paths to data and store 
path_to_data = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/jan_surf_files"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"
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
W2E1_topo  = read_jan_surf_oromea(path=path_to_data, filename=W2E1_filename,)
W2E15_topo  = read_jan_surf_oromea(path=path_to_data, filename=W2E15_filename,)
W2E0_topo  = read_jan_surf_oromea(path=path_to_data, filename=W2E0_filename,)
W2E2_topo = read_jan_surf_oromea(path=path_to_data, filename=W2E2_filename,)



def extract_topo_profile(maxlon=20, minlon=-5, maxlat=47, minlat=46, 
                          dim="lon"):
    
    w1e1_mio = extract_profile(data = W1E1_topo, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, 
                              dim=dim, to_pandas=True)
    w2e1_mio = extract_profile(data = W2E1_topo, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, 
                              dim=dim, to_pandas=True)
    
    w2e0_mio = extract_profile(data = W2E0_topo, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, 
                              dim=dim, to_pandas=True)
    
    w2e2_mio = extract_profile(data = W2E2_topo, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, 
                              dim=dim, to_pandas=True)
    w2e15_mio = extract_profile(data = W2E15_topo, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, 
                              dim=dim, to_pandas=True)
    
    
    topo_data = [w1e1_mio, w2e1_mio, w2e0_mio, w2e15_mio, w2e2_mio]
    topo_names = ["W1E1 (MIO)", "W2E1 (MIO)", "W2E0 (MIO)", "W2E1.5 (MIO)" ,"W2E2 (MIO)"]
    
    for i,topo in enumerate(topo_names):
        if i ==0:
            
            df = pd.DataFrame(index=topo_data[i].index.values, columns=topo_names)
        df[topo] = topo_data[i]
        
    return df

def plot_topo_profiles(varname, units, data_topo, ax=None, path_to_store=None, filename=None,
                       colors=None, xlabel=True, ylabel=True, title=None, ax_legend=True,
                       ymin=None, ymax=None, dim="lon"):
    
    topo_names = ["W1E1 (MIO)", "W2E1 (MIO)", "W2E0 (MIO)", "W2E1.5 (MIO)" ,"W2E2 (MIO)"]
    
    if ax is None:
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
        
    if colors is None:
        colors = ['#21130d', '#2ca02c', '#1f77b4', '#fc03db', '#d62728']
    data_topo.plot(ax=ax, linestyle="-", color=colors, linewidth=3)    #marker="o", markersize=7
            
            
    if ylabel:
        ax.set_ylabel(varname + " [" + units + "]", fontweight="bold", fontsize=20)
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_yticklabels([])
    
    if xlabel is not None:
        if dim == "lon":
             ax.set_xlabel("Longitude [E°]", fontsize=22, fontweight="bold")
        elif dim == "lat":
             ax.set_xlabel("Latitude [N°]", fontsize=22, fontweight="bold")
        else:
            raise ValueError("Define dim as lat or lon")
        
        
        ax.grid(True, linestyle="--", color="grey")
    else:
        ax.grid(True, linestyle="--", color="grey")
        ax.set_xticklabels([])
        
    if all(parameter is not None for parameter in [ymax, ymin]):
        ax.set_ylim(ymin, ymax)
        
    if ax_legend:
        ax.legend(frameon=True, fontsize=22, loc="upper right")
    else:
        ax.legend([],[], frameon=False)
        
       
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="center")
        
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
    
    if path_to_store is not None:
        plt.savefig(os.path.join(path_to_store, filename), bbox_inches="tight", format= "svg")


df_lon = extract_topo_profile(maxlon=20, minlon=-5, maxlat=47, minlat=46, 
                          dim="lon")

df_lat = extract_topo_profile(maxlon=10, minlon=9, maxlat=54, minlat=40, 
                          dim="lat")


apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(18, 14), sharey=True)
plot_topo_profiles(varname="Elevation", units="m", data_topo=df_lon, ax=ax1,
                   ax_legend=True,ymax=3500, ymin=0, dim="lon")

plot_topo_profiles(varname="Elevation", units="m", data_topo=df_lat, ax=ax2,
                   ax_legend=False,ymax=3500, ymin=0, dim="lat")

plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.05)
plt.savefig(os.path.join(path_to_plots, "topo_profiles.svg"), format= "svg", bbox_inches="tight", dpi=600)
plt.show()

def plot_line(ax):
    
    projection = ccrs.PlateCarree()
    ax.plot([-5, 20],[46.5, 46.5],transform=projection, color="black", linestyle="--",
            )
    ax.plot([10, 10],[40, 54],transform=projection, color="red", linestyle="--",
            )

def add_patches(ax):
    projection = ccrs.PlateCarree()
    #west --> 44, 47 N 1, 8 E
    lat_w, h_w = 44 , 3  # lat and height
    lon_w , w_w = 1, 7  # long and width 
    
    #Europe --> 43, 51 N 2E, 16 E
    lat_e, h_e = 42, 9
    lon_e, w_e = 0, 19
    
    #north --> 46.5, 50 N 5, 16 E
    lat_n, h_n = 46.5, 4
    lon_n, w_n = 5, 11
    
    # south--> 43, 47 N 7.5, 15 E
    lat_s, h_s = 43, 4
    lon_s, w_s = 8, 7.5
    
    ax.add_patch(patches.Rectangle(xy =(lon_w, lat_w), width= w_w, height=h_w, ls= "-", color= red, transform = projection, 
                                fc="None", lw=2.5,))

    ax.add_patch(patches.Rectangle(xy =(lon_n, lat_n), width= w_n, height=h_n, ls= "-", color= golden, transform = projection, 
                                    fc="None", lw=2.5))
    
    ax.add_patch(patches.Rectangle(xy =(lon_s, lat_s), width= w_s, height=h_s, ls= "-", color= magenta, transform = projection, 
                                    fc="None", lw=2.5))
    
    ax.add_patch(patches.Rectangle(xy =(lon_e, lat_e), width= w_e, height=h_e, ls= "-", color= black, transform = projection, 
                                    fc="None", lw=2.5))
    
    
    
def plot_topo():    
    projection = ccrs.PlateCarree()
    fig, ((ax1,ax2,ax3), (ax4, ax5,ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(26, 17), 
                                                        subplot_kw={"projection": projection})
    
    norm, terrain = creat_norm()
    
    plot_echam_topo(variable="Elevation", data=CTL_topo, cmap=terrain, units="m", 
                    vmax=4000, vmin=-100, levels=31, level_ticks=6,
                    domain="Europe", cbar=True, cbar_position= [0.35, 0.05, 0.25, 0.02], cbar_orientation="horizontal",
                    projection=projection, norm=norm, plot_coastlines=True, plot_borders=False, ax=ax1, 
                    title="CTL (PI)", bottom_labels=False, fig=fig)
    
    
    add_patches(ax1)
    plot_echam_topo(variable="Elevation", data=W1E1_topo, cmap=terrain, units="m", 
                    vmax=4000, vmin=-100, levels=31, level_ticks=6,
                    domain="Europe", cbar=False, projection=projection, norm=norm, plot_coastlines=False, plot_borders=False, 
                    ax=ax2, title="W1E1 (MIO)", bottom_labels=False, left_labels=False, sea_land_mask= Mio_slm)
    
    plot_line(ax=ax2)
    
    plot_echam_topo(variable="Elevation", data=W2E2_topo, cmap=terrain, units="m", 
                    vmax=4000, vmin=-100, levels=31, level_ticks=6,
                    domain="Europe", cbar=False, projection=projection, norm=norm, plot_coastlines=False, plot_borders=False, 
                    ax=ax3, title="W2E2 (MIO)", bottom_labels=False, left_labels=False, sea_land_mask= Mio_slm)
    
    plot_echam_topo(variable="Elevation", data=W2E0_topo, cmap=terrain, units="m", 
                    vmax=4000, vmin=-100, levels=31, level_ticks=6,
                    domain="Europe", cbar=False, projection=projection, norm=norm, plot_coastlines=False, plot_borders=False, 
                    ax=ax4, title="W2E0 (MIO)", bottom_labels=True, left_labels=True, sea_land_mask= Mio_slm)
    
    plot_echam_topo(variable="Elevation", data=W2E1_topo, cmap=terrain, units="m", 
                    vmax=4000, vmin=-100, levels=31, level_ticks=6,
                    domain="Europe", cbar=False, projection=projection, norm=norm, plot_coastlines=False, plot_borders=False, 
                    ax=ax5, title="W2E1 (MIO)", bottom_labels=True, left_labels=False, sea_land_mask= Mio_slm)
    
    plot_echam_topo(variable="Elevation", data=W2E15_topo, cmap=terrain, units="m", 
                    vmax=4000, vmin=-100, levels=31, level_ticks=6,
                    domain="Europe", cbar=False, projection=projection, norm=norm, plot_coastlines=False, plot_borders=False, 
                    ax=ax6, title="W2E1.5 (MIO)", bottom_labels=True, left_labels=False, sea_land_mask= Mio_slm)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
    plt.savefig(os.path.join(path_to_plots, "topography_mio.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.show()
    
    
plot_topo()