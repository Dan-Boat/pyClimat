# -*- coding: utf-8 -*-
"""
Created on Sat May 20 17:27:36 2023

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


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_ECHAM_processed
from pyClimat.analysis import extract_var, compute_lterm_mean, compute_lterm_diff


path_to_data = "D:/Datasets/Model_output_pst"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM


# file names 

W1E1_PI_filename = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
W1E1_278_filename = "a015_hpc-bw_e5w2.3_t159_MIO_W1E1_278ppm_t159l31.6h"
W1E1_450_filename = "a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h"

W2E1_PI_filename = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
W2E1_Mio278_filename = "a017_hpc-bw_e5w2.3_t159_MIO_W2E1_278ppm_t159l31.6h"
W2E1_Mio450_filename = "a016_hpc-bw_e5w2.3_t159_MIO_W2E1_450ppm_t159l31.6h"



# read data (long-term means)
years = "1003_1017"
period = "1m"


W1E1_PI_data, W1E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_278_data, W1E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_278_filename, years=years, 
                                          period=period, read_wiso=True)

W1E1_450_data, W1E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W1E1_450_filename, years=years, 
                                          period=period, read_wiso=True)


W2E1_PI_data, W2E1_PI_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_PI_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_278_data, W2E1_278_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio278_filename, years=years, 
                                          period=period, read_wiso=True)

W2E1_450_data, W2E1_450_wiso = read_ECHAM_processed(main_path=path_to_data, exp_name=W2E1_Mio450_filename, years=years, 
                                          period=period, read_wiso=True)


def extract_winds_and_analysis(exp_data, W1E1_data, return_ctl=False):
    
    # extract prec
    v10_data_ctl = extract_var(Dataset=W1E1_data, varname="v10")
    u10_data_ctl = extract_var(Dataset=W1E1_data, varname="u10")
    wind10_data_ctl = extract_var(Dataset=W1E1_data, varname="wind10")
    elev_data_ctl = extract_var(Dataset=W1E1_data, varname="elev")
    
    v10_data = extract_var(Dataset=exp_data, varname="v10")
    u10_data = extract_var(Dataset=exp_data, varname="u10")
    wind10_data = extract_var(Dataset=exp_data, varname="wind10")
    elev_data = extract_var(Dataset=exp_data, varname="elev")
    
    
    exp_v10_alt = compute_lterm_mean(data=v10_data, time="annual")
    exp_u10_alt = compute_lterm_mean(data=u10_data, time="annual")
    exp_wind10_alt = compute_lterm_mean(data=wind10_data, time="annual")
    
    # v10_diff_alt = compute_lterm_diff(data_control=v10_data_ctl, data_main=v10_data,
    #                                       time="annual")
    # u10_diff_alt = compute_lterm_diff(data_control=u10_data_ctl, data_main=u10_data,
    #                                       time="annual")
    
    # wind10_diff_alt = compute_lterm_diff(data_control=wind10_data_ctl, data_main=wind10_data,
    #                                       time="annual")
    
    exp_elev_alt = compute_lterm_mean(data=elev_data, time="annual")
    
    
    if return_ctl:
        v10_alt = compute_lterm_mean(data=v10_data_ctl, time="annual")
        u10_alt = compute_lterm_mean(data=u10_data_ctl, time="annual")
        wind10_alt = compute_lterm_mean(data=wind10_data_ctl, time="annual")
        elev_alt = compute_lterm_mean(data=elev_data_ctl, time="annual")
        
        return {"ctl_v10":v10_alt, "ctl_u10":u10_alt, "ctl_wind10":wind10_alt, "ctl_elev":elev_alt,
                "exp_v10":exp_v10_alt, "exp_u10":exp_u10_alt,
                "exp_wind10": exp_wind10_alt, "exp_elev":exp_elev_alt}
    else:
        return {"exp_v10":exp_v10_alt, "exp_u10":exp_u10_alt, 
                "exp_wind10": exp_wind10_alt, "exp_elev":exp_elev_alt}
    

W2E1_PI_diff = extract_winds_and_analysis(exp_data=W2E1_PI_data, W1E1_data=W1E1_PI_data, return_ctl=True)
W2E1_278_diff = extract_winds_and_analysis(exp_data=W2E1_278_data, W1E1_data=W1E1_278_data, return_ctl=True)
W2E1_450_diff = extract_winds_and_analysis(exp_data=W2E1_450_data, W1E1_data=W1E1_450_data, return_ctl=True)




def plot_winds(data_u, data_v, data, cmap=None, wind_scale=40, speed=None, ax=None, fig=None, plot_projection=None, add_colorbar=False,
               domain="Europe", left_labels=False, right_labels=False, bottom_labels=False, cbar_pos=None,
               plot_coastlines=False, coast_resolution=None, sea_land_mask=None, plot_borders=False,
               use_quiver=True, use_streamplot=False, show_arrow_scale=True, vmax=None, vmin=None, level_ticks=None,
               levels=None, cbar_orientation="horizontal", use_colorbar_default=True, units=None, variable=None,
               title=None):
    
    projection = ccrs.PlateCarree()
    
    if plot_projection is None:
        plot_projection = projection
        
    #generating plot using geoaxis predefined or from default
    if ax is None:
        fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":plot_projection})
        
    if add_colorbar == True:    
        if cbar_pos is None:
            cbar_pos = [0.90, 0.30, 0.03, 0.45]
        
        cbar_ax = fig.add_axes(cbar_pos)   # axis for subplot colorbar # left, bottom, width, height
        
        if cbar_orientation == "vertical":
            cbar_ax.get_xaxis().set_visible(False)
            cbar_ax.yaxis.set_ticks_position('right')
            cbar_ax.set_yticklabels([])
            cbar_ax.tick_params(size=0)
        else:
            cbar_ax.get_yaxis().set_visible(False)
            cbar_ax.xaxis.set_ticks_position('bottom')
            cbar_ax.set_xticklabels([])
            cbar_ax.tick_params(size=0)
     
                
    ticks = np.linspace(0, vmax, level_ticks)
    
    if add_colorbar:
        p = data.plot.imshow(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, transform = projection,
                                    cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": cbar_orientation, 
                                                  "shrink": 0.5, "format": "%.0f", "ticks":ticks}, extend= "neither", add_colorbar=add_colorbar, 
                                    cbar_ax=cbar_ax,
                                    add_labels=False)
                
                
        p.colorbar.set_label(label=variable + " [" + units + "]", size= 20, fontweight="bold")
        p.colorbar.ax.tick_params(labelsize=20, size=0,)
        
    else:
        p = data.plot.imshow(ax =ax, cmap=cmap, transform = projection, vmin=vmin, vmax=vmax, levels=levels,
                                 add_colorbar=False, add_labels=False)
        
            
    # ploting background extent
    plot_background(p, domain= domain, left_labels=left_labels, bottom_labels=bottom_labels,
                    plot_coastlines=plot_coastlines, coast_resolution=coast_resolution, plot_borders=plot_borders,
                    coast_color="lightgray")


        
    x = data_v.coords["lon"].data
    y = data_v.coords["lat"].data

    u = data_u.data
    v = data_v.data
    
    
    if use_quiver:
        X,Y = np.meshgrid(x,y)
        skip = (slice(None, None, 3), slice(None, None, 3))  #for extracting the data on interval or use data[::3, ::3]
        
        # ploting winds using quiver 
        q = ax.quiver(X[skip], Y[skip], u[skip], v[skip], transform=projection,  pivot= "mid", scale=wind_scale,
                      headwidth=6, headlength=10, headaxislength=6,color="black")
        
        
        if show_arrow_scale:
            qk = ax.quiverkey(q, 0.90, -0.1, 2, r'$4 \frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties=
                              {"size": 20, "weight":"bold"})
            
            
    elif use_streamplot:
        strm = ax.streamplot(x, y, u, v, transform=projection, density=4, )
        
        
        
    if plot_coastlines==False:
        if sea_land_mask is not None:
            sea_land_mask.plot.contour(colors="lightgray", linestyles="-", ax=ax, transform=projection, levels=[0], linewidths=2.0,
                                       add_labels=False)
            
    if title is not None:
        ax.set_title(title, fontsize=20, weight="bold", loc="left")
    


projection = ccrs.Robinson(central_longitude=0, globe=None)
apply_style2(fontsize=22, style=None, linewidth=2.5) 

norm, terrain = creat_norm()

fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex=False, figsize= (26, 15), subplot_kw= {"projection":projection})

plot_winds(data_u=W2E1_PI_diff.get("ctl_u10"), data_v=W2E1_PI_diff.get("ctl_v10"), wind_scale=40,
           speed=W2E1_PI_diff.get("ctl_wind10"), cmap=terrain, data=W2E1_PI_diff.get("ctl_elev"),
           plot_coastlines=True, bottom_labels=False, left_labels=True, fig=fig, plot_borders=False,
           plot_projection=projection, domain="Europe", use_quiver=True,
           units="m", vmax=4000, vmin=-100, levels=31, level_ticks=6, add_colorbar=True, cbar_pos= [0.35, 0.05, 0.45, 0.02], 
                           cbar_orientation="horizontal", variable="Elevation", ax=ax1, title="(a) W1E1 (PI)") 


plot_winds(data_u=W2E1_278_diff.get("ctl_u10"), data_v=W2E1_278_diff.get("ctl_v10"), wind_scale=40,
           speed=W2E1_278_diff.get("ctl_wind10"), cmap=terrain, data=W2E1_278_diff.get("ctl_elev"),
           plot_coastlines=False, bottom_labels=False, left_labels=False, fig=fig, plot_borders=False,
           sea_land_mask=mio_slm, plot_projection=projection, domain="Europe", use_quiver=True,
           units="m", vmax=4000, vmin=-100, levels=31, level_ticks=6, add_colorbar=False, variable="Elevation", ax=ax2,
           title="(b) W1E1 (Mio278)") 

plot_winds(data_u=W2E1_450_diff.get("ctl_u10"), data_v=W2E1_450_diff.get("ctl_v10"), wind_scale=40,
           speed=W2E1_450_diff.get("ctl_wind10"), cmap=terrain, data=W2E1_450_diff.get("ctl_elev"),
           plot_coastlines=False, bottom_labels=False, left_labels=False, fig=fig, plot_borders=False,
           sea_land_mask=mio_slm, plot_projection=projection, domain="Europe", use_quiver=True,
           units="m", vmax=4000, vmin=-100, levels=31, level_ticks=6, add_colorbar=False, variable="Elevation", ax=ax3, 
           title="(c) W1E1 (Mio450)")      


plot_winds(data_u=W2E1_PI_diff.get("exp_u10"), data_v=W2E1_PI_diff.get("exp_v10"), wind_scale=40,
           speed=W2E1_PI_diff.get("exp_wind10"), cmap=terrain, data=W2E1_PI_diff.get("exp_elev"),
           plot_coastlines=True, bottom_labels=True, left_labels=True, fig=fig, plot_borders=False,
           plot_projection=projection, domain="Europe", use_quiver=True,
           units="m", vmax=4000, vmin=-100, levels=31, level_ticks=6, add_colorbar=False, 
           variable="Elevation", ax=ax4, title="(d) W2E1 (PI)") 


plot_winds(data_u=W2E1_278_diff.get("exp_u10"), data_v=W2E1_278_diff.get("exp_v10"), wind_scale=40,
           speed=W2E1_278_diff.get("exp_wind10"), cmap=terrain, data=W2E1_278_diff.get("exp_elev"),
           plot_coastlines=False, bottom_labels=True, left_labels=False, fig=fig, plot_borders=False,
           sea_land_mask=mio_slm, plot_projection=projection, domain="Europe", use_quiver=True,
           units="m", vmax=4000, vmin=-100, levels=31, level_ticks=6, add_colorbar=False, variable="Elevation", ax=ax5,
           title="(e) W2E1 (Mio278)") 

plot_winds(data_u=W2E1_450_diff.get("exp_u10"), data_v=W2E1_450_diff.get("exp_v10"), wind_scale=40,
           speed=W2E1_450_diff.get("exp_wind10"), cmap=terrain, data=W2E1_450_diff.get("exp_elev"),
           plot_coastlines=False, bottom_labels=True, left_labels=False, fig=fig, plot_borders=False,
           sea_land_mask=mio_slm, plot_projection=projection, domain="Europe", use_quiver=True,
           units="m", vmax=4000, vmin=-100, levels=31, level_ticks=6, add_colorbar=False, variable="Elevation", ax=ax6,
           title="(f) W2E1 (Mio450)")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.10, wspace=0.05)
plt.savefig(os.path.join(path_to_plots, "winds_on_topo.pdf"), format= "pdf", bbox_inches="tight", dpi=300)                                        