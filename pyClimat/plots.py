#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 18:49:46 2021

@author: dboateng
This module contains all the functions required for generating annual, seasonal and monthly plots. It also contains all the analysis plots 
like isotopic profile plots, lapse rate scatter plots
"""
# Import modules
import xarray as xr
import os
import  pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.mpl.ticker import (LatitudeLocator, LongitudeLocator) 
from cartopy.util import add_cyclic_point
import calendar


#Import related package modules
try:
    from .plot_utils import *
    from .analysis import *
except:
    from plot_utils import *
    from analysis import *

# annual plots 

def plot_annual_mean(variable, data_alt, cmap, units, ax=None, vmax=None, vmin=None, levels=None, domain=None, center= True, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None, data_v=None, data_u=None, GNIP_data=None,
                     left_labels= True, bottom_labels=True, add_colorbar=True, plot_stats= False, compare_data1=None, compare_data2=None, max_pvalue=None,
                     hatches=None, fig=None, cbar_pos=None, use_colorbar_default=False, plot_winds=False,
                     orientation = "horizontal", time=None, plot_projection=None):
    """
    

    Parameters
    ----------
    variable : TYPE: str
        DESCRIPTION. The variable to be plotted. Note, it will be display as colorbar name
    data_alt : TYPE: datarray
        DESCRIPTION. The processed data to be visualized
    cmap : TYPE: plt.cmap 
        DESCRIPTION. Color map handle from matplotlib
    units : TYPE: str
        DESCRIPTION. The unit of the dataset to be visualized 
    ax : TYPE: GeoAxis using Matplotlib, optional or defined in control script if subplots are required for different variables
        DESCRIPTION. The default is None. Figure handle to contain plot 
    vmax : TYPE: float, optional
        DESCRIPTION. The default is None. maximum value limit of the variable to be ploted 
    vmin : TYPE: float, optional
        DESCRIPTION. The default is None. minimum value limit of the variable to be ploted 
    levels : TYPE: float, optional
        DESCRIPTION. The default is None. the number of levels for colorbar scale
    domain : TYPE: str, optional
        DESCRIPTION. The default is None. eg. Africa, Asia, Europe
    output_name : TYPE: str, optional
        DESCRIPTION. The default is None. Filename of generated figures
    output_format : TYPE: str, optional
        DESCRIPTION. The default is None. Format to save figure eg. pdf, svg, tiff
    level_ticks : TYPE: float, optional
        DESCRIPTION. The default is None. Interval of ticks for colorbar
    title : TYPE: str, optional
        DESCRIPTION. The default is None. Title of plots
    path_to_store : TYPE: str, optional
        DESCRIPTION. The default is None. Directory to store data
        
    data_v10 = datarray (required for ploting winds)
    data_u10 = datarray (required for ploting winds)
    
    GNIP_data = DataFrame with lon, lat and d18Op for plotting a scatter circles with filled colormap 
    left_labels: TYPE: Boolean, Default is True
        DESCRIPTION. To add lat coordinates on the left of the plots, optioanl 
    bottom_labels: TYPE: Boolean, Default is True
        DESCRIPTION. To add lon coordinates on the bottom of the plots, optioanl 
    add_colorbar: TYPE: Boolean, Default is True
        DESCRIPTION. To add colormap to the plot

    plot_stats: TYPE: Boolean, Default
        DESCRIPTION: plot the statiscal difference between two varied datasets

    compare_data1: TYPE: datarray
        DESCRIPTION: dataset 1 if plot_stats == true
    
    compare_data2: TYPE: datarray
        DESCRIPTION: dataset 2 if plot_stats == true
    center: TYPE: Boolean, True to apply norm for centering zero

    max_pvalue: TYPE: float, optional
        DESCRIPTION: pvalue for the student t-test significance testing

    hatches: TYPE: str, optional:
        DESCRIPTION: hatches from matplotlib 

    Returns
    -------
    None.

    """
    norm = MidpointNormalize(midpoint = 0)
    projection = ccrs.PlateCarree()
    
    if plot_projection is None:
        plot_projection = ccrs.PlateCarree()
        
    #generating plot using geoaxis predefined or from default
    if ax is None:
        fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":plot_projection})
        
    if add_colorbar == True:    
        if cbar_pos is None:
            cbar_pos = [0.90, 0.30, 0.03, 0.45]
        
        if use_colorbar_default == False:
            
            
            cbar_ax = fig.add_axes(cbar_pos)   # axis for subplot colorbar # left, bottom, width, height
            
            if orientation == "vertical":
                cbar_ax.get_xaxis().set_visible(False)
                cbar_ax.yaxis.set_ticks_position('right')
                cbar_ax.set_yticklabels([])
                cbar_ax.tick_params(size=0)
            else:
                cbar_ax.get_yaxis().set_visible(False)
                cbar_ax.xaxis.set_ticks_position('bottom')
                cbar_ax.set_xticklabels([])
                cbar_ax.tick_params(size=0)
        
    
    
    if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
        ticks = np.linspace(vmin, vmax, level_ticks)
        if vmin < 0:
            
            if center==True:
                if add_colorbar ==True:
                    
                    if use_colorbar_default == True:
                        
                        p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, center=0, 
                                        levels=levels, transform = projection, norm=norm, 
                                        cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": orientation, 
                                                      "shrink": 0.70, "format": "%.0f", "ticks":ticks}, extend= "neither",
                                        add_colorbar=True, add_labels=False)
                    else:
                        
            
                         p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, center=0, 
                                         levels=levels, transform = projection, norm=norm, 
                                         cbar_kwargs= {"pad":0.05, "drawedges": True, "orientation": orientation, 
                                                       "shrink": 0.30, "format": "%.0f", "ticks":ticks}, extend= "neither",
                                         add_colorbar=True, cbar_ax = cbar_ax, add_labels=False)
                else:
                    p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, center=0, 
                                    levels=levels, transform = projection, norm=norm, add_colorbar=False, add_labels=False) 
                                   
                    
            else:
                if add_colorbar == True:
                    if use_colorbar_default == True:
                        
                        p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                        levels=levels, transform = projection, 
                                        cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": orientation, 
                                                      "shrink": 0.70, "format": "%.0f", "ticks":ticks}, extend= "neither", 
                                        add_colorbar=True, add_labels=False)
                    else:
                        
                        p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                        levels=levels, transform = projection, 
                                        cbar_kwargs= {"pad":0.05, "drawedges": True, "orientation": orientation, 
                                                      "shrink": 0.30, "format": "%.0f", "ticks":ticks}, extend= "neither", 
                                        add_colorbar=True, cbar_ax = cbar_ax, add_labels=False)
                else:
                    p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                    levels=levels, transform = projection, add_colorbar=False, add_labels=False,)
                                    
        else:
            if add_colorbar == True:
                
                if use_colorbar_default == True:
                    p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                         levels=levels, transform = projection, 
                                         cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": orientation, 
                                                       "shrink": 0.70, "format": "%.0f", "ticks":ticks}, extend= "both",
                                         add_colorbar=True, add_labels=False)
                else:
                    
                    
                
                    p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                         levels=levels, transform = projection, 
                                         cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": orientation, 
                                                       "shrink": 0.70, "format": "%.0f", "ticks":ticks}, extend= "both",
                                         add_colorbar=True, cbar_ax = cbar_ax, add_labels=False)
            else:
                p = data_alt.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                                     levels=levels, transform = projection, add_colorbar=False, add_labels=False)
                                     
    # when limits are not defined for the plot            
    else:
        p = data_alt.plot.imshow(ax =ax, cmap=cmap, transform = projection, 
                                 cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": orientation, 
                                               "shrink": 0.70, "format": "%.0f", "ticks":ticks}, extend= "neither", 
                                 add_colorbar=True, cbar_ax = cbar_ax, add_labels=False)
    
    
    if add_colorbar == True:
        
        p.colorbar.set_label(label=variable + " [" + units + "]", size= 20, fontweight="bold")
        p.colorbar.ax.tick_params(labelsize=20, size=0,)
    
    # ploting background extent
    plot_background(p, domain= domain, left_labels=left_labels, bottom_labels=bottom_labels)
    
    
    if plot_winds == True: 
        
        if all(data is not None for data in [data_v, data_u]):
            # extracting variables for quiver 
            x = data_v.coords["lon"].data
            y = data_v.coords["lat"].data
        
            u = data_u.data
            v = data_v.data
        
            X,Y = np.meshgrid(x,y)
            skip = (slice(None, None, 3), slice(None, None, 3))  #for extracting the data on interval or use data[::3, ::3]
            
            # ploting winds using quiver 
            q = ax.quiver(X[skip], Y[skip], u[skip], v[skip], transform=projection,  pivot= "mid", scale= 100,
                          headwidth=3, headlength=5, headaxislength=4.5)
            
            qk = ax.quiverkey(q, 0.90, -0.1, 5, r'$5 \frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties=
                              {"size": 20, "weight":"bold"})
        
        
    if plot_stats == True:
        data1 = compare_data1
        data2 = compare_data2
        
        if domain == "Europe":
            minlat, maxlat, minlon, maxlon = 35, 65, -15, 40
            stats_results = student_t_test_btn_datasets(dataA=data1, dataB=data2, return_pvalue=True, 
                                                        minlat=minlat, minlon=minlon, maxlon=maxlon, maxlat=maxlat,
                                                        max_pvalue=max_pvalue, time=time)
        elif domain == "West Africa":
            minlat, maxlat, minlon, maxlon = -5, 35, -25, 40
            stats_results = student_t_test_btn_datasets(dataA=data1, dataB=data2, return_pvalue=True, 
                                                        minlat=minlat, minlon=minlon, maxlon=maxlon, maxlat=maxlat,
                                                        max_pvalue=max_pvalue, time=time)
        else:
            stats_results = student_t_test_btn_datasets(dataA=data1, dataB=data2, return_pvalue=True, max_pvalue=max_pvalue, time=time)
        
        if hatches is not None:
            ax.contourf(stats_results.lon.values, stats_results.lat.values, stats_results.t_statistic.values, colors="none", hatches=[hatches])
        else:
            ax.contourf(stats_results.lon.values, stats_results.lat.values, stats_results.t_statistic.values, colors="none", hatches=["//"])
            
            
            
        
    if GNIP_data is not None:
        
        
        if center == True:
            ax.scatter(x=GNIP_data["lon"], y=GNIP_data["lat"], c=GNIP_data["d18op"], cmap=cmap, vmax=vmax, vmin=vmin, norm=norm, edgecolor="k", s= 140)
        else:
            ax.scatter(x=GNIP_data["lon"], y=GNIP_data["lat"], c=GNIP_data["d18op"], cmap=cmap, vmax=vmax, vmin=vmin, edgecolor="k", s= 140)
    
    if title is not None:
        ax.set_title(title, fontsize=20, weight="bold", loc="left")
        
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
    else:
        print("The output would not be save on directory")
        
        
        
        


def plot_seasonal_mean(variable, data_slt, cmap, units, seasons, axes=None, fig=None, vmax=None, vmin=None, levels=None, domain=None, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None, data_v=None,
                     plot_winds_pattern=False, plot_winds_streamline=False,
                     data_u=None, cbar_pos=None, fig_title=None, season_label=None, plot_stats= False, compare_data1=None, compare_data2=None, max_pvalue=None,
                     hatches=None, add_colorbar = True, left_labels= True, bottom_labels=True, show_arrow_scale=True, center=True, 
                     orientation = "vertical"):
    """
    

    Parameters
    ----------
    variable : TYPE: str
        DESCRIPTION. The variable to be plotted. Note, it will be display as colorbar name
    data_slt : TYPE: datarray
        DESCRIPTION. The processed data to be visualized (must contain the season time coordinate)
    cmap : TYPE: plt.cmap 
        DESCRIPTION. Color map handle from matplotlib
    units : TYPE: str
        DESCRIPTION. The unit of the dataset to be visualized 
    seasons : TYPE: List containing str
        DESCRIPTION.List of seasons to be plotted eg. ["JJA", "DJF"] or ["JJA] or list of all seasons
    axes : TYPE, optional
        DESCRIPTION. The default is None.
    fig : TYPE, optional
        DESCRIPTION. The default is None.
    vmax : TYPE: float, optional
        DESCRIPTION. The default is None. maximum value limit of the variable to be ploted 
    vmin : TYPE: float, optional
        DESCRIPTION. The default is None. minimum value limit of the variable to be ploted 
    levels : TYPE: float, optional
        DESCRIPTION. The default is None. the number of levels for colorbar scale
    domain : TYPE: str, optional
        DESCRIPTION. The default is None. eg. Africa, Asia, Europe
    output_name : TYPE: str, optional
        DESCRIPTION. The default is None. Filename of generated figures
    output_format : TYPE: str, optional
        DESCRIPTION. The default is Notime="season", season_calendar="standard"ne. Format to save figure eg. pdf, svg, tiff
    level_ticks : TYPE: float, optional
        DESCRIPTION. The default is None. Interval of ticks for colorbar
    title : TYPE: Bolean, optional
        DESCRIPTION. The default is None. Title of plots
    path_to_store : TYPE: str, optional
        DESCRIPTION. The default is None. Directory to store data
    cbar_pos : TYPE: list, optional
        DESCRIPTION. The default is None. the list defing the position of the color bar eg. [0.90, 0.30, 0.02, 0.40]
    fig_title = None
    seasonal_label: str (fro the label of which season)
    plot_stats: TYPE: Boolean, optional 
        DESCRIPTION. The default is False. True for ploting hatching for signifacne difference using student t-test or correlation with spearmanr cor
    Compare_data1, compare_data2: TYPE: datarray (not optional if plot_stats is set True)
        DESCRIPTION. the datasets required for statistic computation
    hatches: TYPE: str
        DESCRIPTION. the hatche style require for plotting..must be list in matplotlib hatch handle
    max_pvalue: TYPE: float
        DESCRIPTION. The confidence interval range for statistics significance (eg. 0.05 for 95% CI)
    plot_winds_pattern: TYPE: Boolean, optional 
        DESCRIPTION: It plots the winds pattern using arrows on the plot background
    plot_winds_streamline: TYPE: Boolean, optional 
        DESCRIPTION: It plots the wind streamlines on the plot
        
    data_v = datarray (required for ploting winds)
    data_u = datarray (required for ploting winds)
    
    plot_stats: TYPE: Boolean, Default
        DESCRIPTION: plot the statiscal difference between two varied datasets

    compare_data1: TYPE: datarray
        DESCRIPTION: dataset 1 if plot_stats == true
    
    compare_data2: TYPE: datarray
        DESCRIPTION: dataset 2 if plot_stats == true
    center: TYPE: Boolean, True to apply norm for centering zero

    max_pvalue: TYPE: float, optional
        DESCRIPTION: pvalue for the student t-test significance testing

    hatches: TYPE: str, optional:
        DESCRIPTION: hatches from matplotlib 


    Returns
    -------
    None.

    """
    
    norm = MidpointNormalize(midpoint = 0)
    projection = ccrs.PlateCarree()
    
    if axes is None:
        if len(seasons) == 1:
            fig, ax1 = plt.subplots(nrows=1, ncols=1, sharex= True, sharey= True, figsize=(8, 7),
                               subplot_kw={"projection": projection})
            axes = [ax1]
        elif len(seasons) == 2:
            fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex= True, sharey= True, figsize=(8, 10), 
                               subplot_kw={"projection": projection})
            axes = [ax1, ax2]
        elif len(seasons) == 4:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex= True, sharey= True, figsize=(11, 8), 
                               subplot_kw={"projection": projection})
            axes = [ax1, ax2, ax3, ax4]
    cbar_axis = axes[-1]
    
    for i,season in enumerate(seasons):
        if add_colorbar == True:
            
            if cbar_pos is None:
                cbar_pos = [0.90, 0.30, 0.03, 0.45]
    
            cbar_ax = fig.add_axes(cbar_pos)   # axis for subplot colorbar # left, bottom, width, height
            cbar_ax.get_xaxis().set_visible(False)
            cbar_ax.yaxis.set_ticks_position('right')
            cbar_ax.set_yticklabels([])
            cbar_ax.tick_params(size=0)
        
        if axes[i]==cbar_axis:
            if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
                ticks = np.linspace(vmin, vmax, level_ticks)
                if vmin < 0 & center==True:
                    print("---using customized norm for the colormap ------")
                    if add_colorbar ==True:
                        p = data_slt.sel(season=season).plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, center=0, 
                                 levels=levels, transform = projection, norm=norm, 
                                 cbar_kwargs= {"pad":0.05, "drawedges": True, "orientation": orientation, 
                                               "shrink": 0.30, "format": "%.0f", "ticks":ticks}, extend= "neither",
                                 add_colorbar=True, cbar_ax = cbar_ax, add_labels=False)
                    else:
                        p = data_slt.sel(season=season).plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, center=0, 
                                 levels=levels, transform = projection, norm=norm, 
                                 add_colorbar=False, add_labels=False)
                else:
                    print("-----skipping the use of norm for the cmap -------")
                    if add_colorbar == True:
                        
                        p = data_slt.sel(season=season).plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, 
                                     levels=levels, transform = projection, 
                                     cbar_kwargs= {"pad":0.05, "drawedges": True, "orientation": orientation, 
                                                   "shrink": 0.30, "format": "%.0f", "ticks":ticks}, extend= "neither",
                                     add_colorbar=True, cbar_ax=cbar_ax, add_labels=False) 
                    else: 
                        p = data_slt.sel(season=season).plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, 
                                     levels=levels, transform = projection, add_colorbar=False, add_labels=False)
                                     
            else:
                p = data_slt.sel(season=season).plot.imshow(ax =axes[i], cmap=cmap, transform = projection, 
                                 cbar_kwargs= {"pad":0.05, "drawedges": True, "orientation": orientation, 
                                               "shrink": 0.30, "format": "%.0f", "ticks":ticks}, extend= "neither", add_labels=False)
            if add_colorbar == True:
                p.colorbar.set_label(label=variable + " [" + units + "]", size= 22, fontweight= "bold")
                
                p.colorbar.ax.tick_params(labelsize=20, size=0)
    
            # ploting background extent
            plot_background(p, domain= domain, left_labels=left_labels, bottom_labels=bottom_labels)
            
            if plot_winds_pattern == True:
            
                if all(data is not None for data in [data_v, data_u]):
                    # extracting variables for quiver 
                    x = data_v.coords["lon"].data
                    y = data_v.coords["lat"].data
                
                    u = data_u.sel(season=season).data
                    v = data_v.sel(season=season).data
                
                    X,Y = np.meshgrid(x,y)
                    skip = (slice(None, None, 3), slice(None, None, 3))  #for extracting the data on interval or use data[::3, ::3]
                    
                    # ploting winds using quiver 
                    q = axes[i].quiver(X[skip], Y[skip], u[skip], v[skip], transform=projection,  pivot= "mid", scale= 80,
                                  headwidth=3, headlength=5, headaxislength=4.5)
                    
                    if show_arrow_scale==True:
                        qk = axes[i].quiverkey(q, 1.0, -0.02, 5, r'$5 \frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties=
                                               {"size": 22, "weight":"bold"})
            
        else:
            if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
                ticks = np.linspace(vmin, vmax, level_ticks)
                if vmin < 0:
                    p = data_slt.sel(season=season).plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, center=0, 
                                 levels=levels, transform = projection, norm=norm, extend= "neither", add_colorbar=False, add_labels=False)
                else:
                    p = data_slt.sel(season=season).plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, 
                                 levels=levels, transform = projection, extend= "neither", add_colorbar=False, add_labels=False) 
            else:
                p = data_slt.sel(season=season).plot.imshow(ax =axes[i], cmap=cmap, transform = projection, extend= "neither", add_labels=False)
                
            # ploting background extent
            plot_background(p, domain= domain, left_labels=left_labels, bottom_labels=bottom_labels)
            
            if plot_winds_pattern == True:
                
            
                if all(data is not None for data in [data_v, data_u]):
                    # extracting variables for quiver 
                    x = data_v.coords["lon"].data
                    y = data_v.coords["lat"].data
                
                    u = data_u.sel(season=season).data
                    v = data_v.sel(season=season).data
                
                    X,Y = np.meshgrid(x,y)
                    skip = (slice(None, None, 3), slice(None, None, 3))  #for extracting the data on interval or use data[::3, ::3]
                    
                    # ploting winds using quiver 
                    q = axes[i].quiver(X[skip], Y[skip], u[skip], v[skip], transform=projection,  pivot= "mid", scale= 80,
                                  headwidth=3, headlength=5, headaxislength=4.5)
                    
                    if show_arrow_scale==True:
                        qk = axes[i].quiverkey(q, 1.0, -0.02, 5, r'$5 \frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties=
                                               {"size": 22, "weight":"bold"})
                
                
        if plot_stats == True:
            data1 = compare_data1.groupby("time.season")[season]
            data2 = compare_data2.groupby("time.season")[season]
            
            if domain == "Europe":
                minlat, maxlat, minlon, maxlon = 35, 65, -15, 40
                stats_results = student_t_test_btn_datasets(dataA=data1, dataB=data2, return_pvalue=True, 
                                                            minlat=minlat, minlon=minlon, maxlon=maxlon, maxlat=maxlat,
                                                            max_pvalue=max_pvalue)
            else:
                stats_results = student_t_test_btn_datasets(dataA=data1, dataB=data2, return_pvalue=True, max_pvalue=max_pvalue)
            
            if hatches is not None:
                axes[i].contourf(stats_results.lon.values, stats_results.lat.values, stats_results.t_statistic.values, colors="none", hatches=[hatches])
            else:
                axes[i].contourf(stats_results.lon.values, stats_results.lat.values, stats_results.t_statistic.values, colors="none", hatches=["//"])
        
        if plot_winds_streamline == True:
            
            #convert lon to -180 to 180
            data_v = data_v.assign_coords({"lon": (((data_v.lon + 180) % 360) - 180)})
            data_u = data_u.assign_coords({"lon": (((data_u.lon + 180) % 360) - 180)})
            
            
            # extracting variables for quiver 
            x = data_v.coords["lon"].data
            y = data_v.coords["lat"].data
            
            u = data_u.sel(season=season).data
            v = data_v.sel(season=season).data
            
            X,Y = np.meshgrid(x,y)
            skip = (slice(None, None, 3), slice(None, None, 3))  #for extracting the data on interval or use data[::3, ::3]
            
            #ploting streamlines
            
            strm = axes[i].streamplot(X[skip], Y[skip], u[skip], v[skip], transform=projection, color="black", density=1)
            

            
        if title ==True:
            axes[i].set_title(season_label[i], fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="left")
        elif title ==False:
            axes[i].set_title("", fontdict= {"fontsize": 22, "fontweight":"bold"})
            
    if fig_title is not None:
        fig.suptitle(fig_title, fontsize= 20, weight = "bold")
        
        
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.88, top=0.95, bottom=0.06)
        
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
            
        
            


def plot_monthly_mean(variable, data_mlt, cmap, units, months, axes=None, fig=None, vmax=None, vmin=None, levels=None, domain=None, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None, data_v10=None, 
                     data_u10=None, left_labels= True, bottom_labels=True):
    """
    

    Parameters
    ----------
    variable : TYPE: str
        DESCRIPTION. The variable to be plotted. Note, it will be display as colorbar name
    data_slt : TYPE: datarray
        DESCRIPTION. The processed data to be visualized (must contain the season time coordinate)
    cmap : TYPE: plt.cmap 
        DESCRIPTION. Color map handle from matplotlib
    units : TYPE: str
        DESCRIPTION. The unit of the dataset to be visualized 
    months : TYPE: str
        DESCRIPTION. The range of months to visualise eg. Jan-Jun or Ju-Dec
    axes : TYPE, optional
        DESCRIPTION. The default is None.
    fig : TYPE, optional
        DESCRIPTION. The default is None.
    vmax : TYPE: float, optional
        DESCRIPTION. The default is None. maximum value limit of the variable to be ploted 
    vmin : TYPE: float, optional
        DESCRIPTION. The default is None. minimum value limit of the variable to be ploted 
    levels : TYPE: float, optional
        DESCRIPTION. The default is None. the number of levels for colorbar scale
    domain : TYPE: str, optional
        DESCRIPTION. The default is None. eg. Africa, Asia, Europe
    output_name : TYPE: str, optional
        DESCRIPTION. The default is None. Filename of generated figures
    output_format : TYPE: str, optional
        DESCRIPTION. The default is None. Format to save figure eg. pdf, svg, tiff
    level_ticks : TYPE: float, optional
        DESCRIPTION. The default is None. Interval of ticks for colorbar
    title : TYPE: str, optional
        DESCRIPTION. The default is None. Title of plots
    path_to_store : TYPE: str, optional
        DESCRIPTION. The default is None. Directory to store data

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    
    norm = MidpointNormalize(midpoint = 0)
    projection = ccrs.PlateCarree()
    
    if axes is None:
        fig, ((ax1, ax2), (ax3,ax4), (ax5,ax6)) = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(13,13),
                                                               subplot_kw={"projection":projection})
        axes = [ax1, ax2, ax3, ax4, ax5, ax6] # handles for subplots for looping 
        
    cbar_axis = axes[-1]
    
    if months in ["Jan-Jun", "J-J", "January-June"]:
        #index for extracting monthly data
        months_num = [0, 1, 2, 3, 4, 5]
        mnames = ["January", "February", "March", "April", "May", "June"]
        
    elif months in ["July-December", "J-D", "Ju-Dec"]:
        months_num = [6, 7, 8, 9, 10, 11]
        mnames = ["July", "August", "September", "October","November", "December"]
    else:
        raise ValueError("Define the months as a range between Jan-Jun or Ju-Dec")
        
    
    for i,month in enumerate(mnames):
        
        cbar_ax = fig.add_axes([0.90, 0.30, 0.03, 0.45])   # axis for subplot colorbar # left, bottom, width, height
        cbar_ax.get_xaxis().set_visible(False)
        cbar_ax.yaxis.set_ticks_position('right')
        cbar_ax.set_yticklabels([])
        
        if axes[i]==cbar_axis:
            if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
                ticks = np.linspace(vmin, vmax, level_ticks)
                if vmin < 0:
                    p = data_mlt[months_num[i]].plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, center=0, 
                                 levels=levels, transform = projection, norm=norm, 
                                 cbar_kwargs= {"pad":0.05, "drawedges": True, "orientation": "vertical", 
                                               "shrink": 0.30, "format": "%.0f", "ticks":ticks}, extend= "neither",
                                 add_colorbar=True, cbar_ax = cbar_ax)
                else:
                    p = data_mlt[months_num[i]].plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, 
                                 levels=levels, transform = projection, 
                                 cbar_kwargs= {"pad":0.05, "drawedges": True, "orientation": "vertical", 
                                               "shrink": 0.30, "format": "%.0f", "ticks":ticks}, extend= "neither",
                                 add_colorbar=True, cbar_ax=cbar_ax) 
            else:
                p = data_mlt[months_num[i]].plot.imshow(ax =axes[i], cmap=cmap, transform = projection, 
                                 cbar_kwargs= {"pad":0.05, "drawedges": True, "orientation": "vertical", 
                                               "shrink": 0.30, "format": "%.0f", "ticks":ticks}, extend= "neither")
            p.colorbar.set_label(label=variable + " [" + units + "]", size= 20)
            p.colorbar.ax.tick_params(labelsize=20)
    
            # ploting background extent
            plot_background(p, domain= domain, left_labels=left_labels, bottom_labels=bottom_labels)
            
            if all(data is not None for data in [data_v10, data_u10]):
                # extracting variables for quiver 
                x = data_v10.coords["lon"].data
                y = data_v10.coords["lat"].data
            
                u = data_u10[months_num[i]].data
                v = data_v10[months_num[i]].data
            
                X,Y = np.meshgrid(x,y)
                skip = (slice(None, None, 3), slice(None, None, 3))  #for extracting the data on interval or use data[::3, ::3]
                
                # ploting winds using quiver 
                q = axes[i].quiver(X[skip], Y[skip], u[skip], v[skip], transform=projection,  pivot= "mid", scale= 100,
                              headwidth=3, headlength=5, headaxislength=4.5)
                qk = axes[i].quiverkey(q, 1.02, -0.02, 2, r'$1 \frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties=
                          {"size": 20, "weight":"bold"})
            
        else:
            if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
                ticks = np.linspace(vmin, vmax, level_ticks)
                if vmin < 0:
                    p = data_mlt[months_num[i]].plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, center=0, 
                                 levels=levels, transform = projection, norm=norm, extend= "neither", add_colorbar=False)
                else:
                    p = data_mlt[months_num[i]].plot.imshow(ax =axes[i], cmap=cmap, vmin=vmin, vmax=vmax, 
                                 levels=levels, transform = projection, extend= "neither", add_colorbar=False) 
            else:
                p = data_mlt[months_num[i]].plot.imshow(ax =axes[i], cmap=cmap, transform = projection, extend= "neither")
                
            # ploting background extent
            plot_background(p, domain= domain, left_labels=left_labels, bottom_labels=bottom_labels)
            
            if all(data is not None for data in [data_v10, data_u10]):
                # extracting variables for quiver 
                x = data_v10.coords["lon"].data
                y = data_v10.coords["lat"].data
            
                u = data_u10[months_num[i]].data
                v = data_v10[months_num[i]].data
            
                X,Y = np.meshgrid(x,y)
                skip = (slice(None, None, 3), slice(None, None, 3))  #for extracting the data on interval or use data[::3, ::3]
                
                # ploting winds using quiver 
                q = axes[i].quiver(X[skip], Y[skip], u[skip], v[skip], transform=projection,  pivot= "mid", scale= 100,
                              headwidth=3, headlength=5, headaxislength=4.5)
                #qk = ax.quiverkey(q, 0.1, 0.2, 2, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')
            
        if title is not None:
            axes[i].set_title(month, fontdict= {"fontsize": 20, "fontweight":"bold"})
            
    if title is not None:
        fig.suptitle(title, fontsize= 20, weight = "bold")
        
        
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.88, top=0.95, bottom=0.06)
        
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")




def plot_iso_profiles(df_iso, df_geosp, dim, iso_color, iso_label, ax=None, season=None, month=None,  
                      xmax=None, xmin=None, ymax=None, ymin=None, ax_legend=None, isomax=None, isomin=None,
                       output_name=None, output_format=None, title=None, path_to_store=None, left_labels=True,
                       bottom_labels=True, right_labels=True, shade_color=None, shade_alpha=None, edgecolor= "dimgrey"):
    """
    

    Parameters
    ----------
    df_iso : TYPE: DataFrame
        DESCRIPTION. The output from extract_profile functions for isotope
    df_geosp : TYPE: DataFrame
        DESCRIPTION. The output from extract_profile functions for elevation 
    dim : TYPE: str
        DESCRIPTION. The direction of the profile line (whether lat or lon)
    iso_color : TYPE: Matplotlib color handle 
        DESCRIPTION. Color for a specific isotopic profile
    iso_label : TYPE: str
        DESCRIPTION. The lable for module experiment used for constructing isotopic profile
    ax : TYPE: plt axes handle, optional
        DESCRIPTION. The default is None. This must be defined in the control script if multiple experiments are used
    season : TYPE: str, optional
        DESCRIPTION. The default is None. Must be defined if specific season is required 
    month : TYPE: int, optional
        DESCRIPTION. The default is None. he default is None. Must be defined if specific month is required 
    xmax : TYPE: float, optional
        DESCRIPTION. The default is None. The maximum limit of coordinates 
    xmin : TYPE:float, optional
        DESCRIPTION. The default is None. The minimun limit of cordinates
    ymax : TYPE: float, optional
        DESCRIPTION. The default is None. The maximum limit of elevation axis
    ymin : TYPE: float, optional
        DESCRIPTION. The default is None. The minimum limit of elevation axis 
    ax_legend : TYPE: Boolean, optional
        DESCRIPTION. The default is None. True if you want to show legend. Can also be defined as fig.lenged if mutiple data
        are used in the control script. Check the example script
    isomax : TYPE: float, optional
        DESCRIPTION. The default is None. The maximum limit of the iso values 
    isomin : TYPE:float, optional
        DESCRIPTION. The default is None. The minimum limit of the iso values 
    output_name : TYPE: str, optional
        DESCRIPTION. The default is None. Filename of generated figures
    output_format : TYPE: str, optional
        DESCRIPTION. The default is None. Format to save figure eg. pdf, svg, tiff
    level_ticks : TYPE: float, optional
        DESCRIPTION. The default is None. Interval of ticks for colorbar
    title : TYPE: str, optional
        DESCRIPTION. The default is None. Title of plots
    path_to_store : TYPE: str, optional
        DESCRIPTION. The default is None. Directory to store data
        
    left labels, right_labels, bottom_labels, : TYPE: Bol, optional
        DESCRIPTION. To set the left, right, and bottom axis label to None

    shade_color: TYPE: STR
        DESCRIPTION: shade color for plotting fill_between
    
    shade_alpha: TYPE: float, optional
        DESCRIPTION: shade factor for plotting fill_between

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(10, 10))
    if season is not None:
        df_iso = df_iso[season]
        df_geosp = df_geosp[season]
    elif month is not None:
        df_iso = df_iso[month]
        df_geosp = df_geosp[month]
    
    # filling elevation data
    
    if shade_color is not None: 
        if shade_alpha is None:
            ax.fill_between(df_geosp.index, df_geosp,  0, color=shade_color, alpha=0.05, edgecolor = edgecolor, linestyle="-", linewidth=3)
        else:
            ax.fill_between(df_geosp.index, df_geosp,  0, color=shade_color, alpha=shade_alpha, edgecolor = edgecolor, linestyle="-", linewidth=3)
    else:
        
        ax.fill_between(df_geosp.index, df_geosp,  0, color="dimgrey", alpha=0.1, edgecolor = edgecolor, linestyle="-", linewidth=3)
    
    if all(parameter is not None for parameter in [xmax, xmin, ymax, ymin]):
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        
        
    if bottom_labels == True:
        if dim == "lon":
             ax.set_xlabel("Longitude [E°]", fontsize=22)
        elif dim == "lat":
             ax.set_xlabel("Latitude [N°]", fontsize=22)
        else:
            raise ValueError("Define dim as lat or lon")
            
            
    if left_labels == True:
        ax.set_ylabel("Elevation [m]", fontsize=22)
            
    
    if left_labels ==False:
        ax.grid(True)
        #ax.axes.yaxis.set_visible(False)
        ax.set_yticklabels([])
    
        
    if bottom_labels == False:
        ax.grid(True)
        #ax.axes.yaxis.set_visible(False)
        ax.set_xticklabels([])
        
        
        
    ax.tick_params(which="both")
    #creating secondary axis 
    ax2 = ax.twinx()
    ax2.grid(False)
    
    #plotting iso values on the same plot 
    ax2.plot(df_iso.index, df_iso, linestyle= "--", label = iso_label, color = iso_color)
    
    if all(parameter is not None for parameter in [isomax, isomin]):
        ax2.set_ylim(isomin, isomax)
    if right_labels == True:
    
        ax2.set_ylabel(u'$\delta^{18}$O ‰ vs SMOW', fontsize=22)
        
    else:
        ax2.set_yticklabels([])
        
    ax2.tick_params(axis= "y")
    ax2.tick_params(axis= "x")
    
    if ax_legend == True:
        ax2.legend(frameon=True, fontsize=20, loc="upper left")
    
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 20, "fontweight":"bold"}, loc="left")
        
    plt.tight_layout()
    
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
        


def scatter_plot_laspe_rate(reg_params, df_x_y_yhat, color, marker, label, ylabel=None, xlabel=None, ax=None, ax_legend=None,
                            output_name=None, output_format=None, title=None, path_to_store=None,
                            xmax=None, xmin=None, ymax=None, ymin=None, left_labels=True,
                            bottom_labels=True,):
    """
    

    Parameters
    ----------
    reg_params : TYPE: output from stats.linregress
        DESCRIPTION.
    df_x_y_yhat : TYPE: DataFrame output from linear_regression module in Climat_analysis
        DESCRIPTION.
    color : TYPE: plt.color handle 
        DESCRIPTION.
    marker : TYPE: plt.marker handle for scatter
        DESCRIPTION.
    label : TYPE: str
        DESCRIPTION. Additional lable for data aside the equation of line of fitting
    ylabel : TYPE: str
        DESCRIPTION. Y-axis lable name 
    xlabel : TYPE: str
        DESCRIPTION. X-axis label name
     ax : TYPE: plt axes handle, optional
        DESCRIPTION. The default is None. This must be defined in the control script if multiple experiments are used
    ax_legend : TYPE: Boolean, optional
        DESCRIPTION. The default is None. True if you want to show legend. Can also be defined as fig.lenged if mutiple data
        are used in the control script. Check the example script
    output_name : TYPE: str, optional
        DESCRIPTION. The default is None. Filename of generated figures
    output_format : TYPE: str, optional
        DESCRIPTION. The default is None. Format to save figure eg. pdf, svg, tiff
    level_ticks : TYPE: float, optional
        DESCRIPTION. The default is None. Interval of ticks for colorbar
    title : TYPE: str, optional
        DESCRIPTION. The default is None. Title of plots
    path_to_store : TYPE: str, optional
        DESCRIPTION. The default is None. Directory to store data
    
    left labels, right_labels, bottom_labels, : TYPE: Bol, optional
        DESCRIPTION. To set the left, right, and bottom axis label to None

    Returns
    -------
    None.

    """
    
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(10, 10))
    
    #ploting scatter points 
    ax.scatter(df_x_y_yhat["X"], df_x_y_yhat["Y"], color=color, marker=marker)
    ax.plot(df_x_y_yhat["X"], df_x_y_yhat["yhat"], 
            color=color, label= "ILR = {:.2f} [‰/km], r²={:.2f}".format(reg_params.slope*1000, 
                                                                      reg_params.rvalue*-1) + " [" + label + "]")
    
    #ax.set_aspect("equal", "box")
    
    if bottom_labels ==True:
        ax.set_xlabel("Elevation [m]", fontsize=22, fontweight="bold")
        
    elif bottom_labels ==False:
        ax.grid(True)
        ax.set_xticklabels([])
    else:
        raise ValueError("define xlabel or set bottom labels to False")
    
    
    if left_labels ==True:
        ax.set_ylabel(u'$\delta^{18}$O ‰ vs SMOW', fontsize=22, fontweight="bold")
        
    elif left_labels ==False:
        ax.grid(True)
        ax.set_yticklabels([])
    else:
        raise ValueError("define xlabel or set bottom labels to False")
        
    
    
    # if all(parameter is not None for parameter in [xlabel, ylabel]):
    #     ax.set_xlabel(xlabel, fontsize= 20,)
    #     ax.set_ylabel(ylabel, fontsize= 20,)
        
    if all(parameter is not None for parameter in [xmax, xmin, ymax, ymin]):
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    
    if ax_legend =="True":
        ax.legend(frameon=True, fontsize=22, loc="upper right")
       
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="left")
        
    plt.tight_layout()
    
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
        
        
    
def plot_echam_topo(variable, data, cmap, units, ax=None, vmax=None, vmin=None, levels=None, domain=None, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None, cbar = None, cbar_orientation=None, cbar_position = None,
                     fig=None, left_labels= True, bottom_labels=True, norm=None, projection=None, plot_coastlines=True, plot_borders=False, 
                     sea_land_mask=None):
    """
    

    Parameters
    ----------
    variable : TYPE: str
        DESCRIPTION. The variable to be plotted. Note, it will be display as colorbar name
    data : TYPE: datarray
        DESCRIPTION. The processed data to be visualized (Eg. topo input file or can be retrieved from model output)
    cmap : TYPE: plt.cmap 
        DESCRIPTION. Color map handle from matplotlib
    units : TYPE: str
        DESCRIPTION. The unit of the dataset to be visualized 
    ax : TYPE: matplotlib ax handle, optional
        DESCRIPTION. The default is None.
    fig : TYPE: Matplotlib figure handle, optional
        DESCRIPTION. The default is None.
    vmax : TYPE: float, optional
        DESCRIPTION. The default is None. maximum value limit of the variable to be ploted 
    vmin : TYPE: float, optional
        DESCRIPTION. The default is None. minimum value limit of the variable to be ploted 
    levels : TYPE: float, optional
        DESCRIPTION. The default is None. the number of levels for colorbar scale
    domain : TYPE: str, optional
        DESCRIPTION. The default is None. eg. Africa, Asia, Europe
    output_name : TYPE: str, optional
        DESCRIPTION. The default is None. Filename of generated figures
    output_format : TYPE: str, optional
        DESCRIPTION. The default is None. Format to save figure eg. pdf, svg, tiff
    level_ticks : TYPE: float, optional
        DESCRIPTION. The default is None. Interval of ticks for colorbar
    title : TYPE: str, optional
        DESCRIPTION. The default is None. Title of plots
    path_to_store : TYPE: str, optional
        DESCRIPTION. The default is None. Directory to store data
    cbar : TYPE: Boolean, optional
        DESCRIPTION. The default is None. True is the plot require colobar axis
    cbar_orientation : TYPE: , optional
        DESCRIPTION. 
    cbar_position : TYPE: list, optional
        DESCRIPTION. The default is None. The default is None. the list defing the position of the color bar eg. [0.90, 0.30, 0.02, 0.40]
    
    left labels, right_labels, bottom_labels, : TYPE: Bol, optional
        DESCRIPTION. To set the left, right, and bottom axis label to None

    Returns
    -------
    None.

    """
    if norm is None:
        norm = MidpointNormalize(midpoint = 0)
        
    if projection is None:
        projection = ccrs.PlateCarree()
    projection_p = ccrs.PlateCarree()
    #generating plot using geoaxis predefined or from default
    if ax is None:
        fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":projection})
    if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
        ticks = np.linspace(0, vmax, level_ticks)
        if cbar==True:
            
            cbar_pos = cbar_position
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
            
            p = data.plot.imshow(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, transform = projection_p,
                                 cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": cbar_orientation, 
                                               "shrink": 0.5, "format": "%.0f", "ticks":ticks}, extend= "neither", add_colorbar=True, cbar_ax=cbar_ax,
                                 add_labels=False)
            
            
            p.colorbar.set_label(label=variable + " [" + units + "]", size= 20, fontweight="bold")
            p.colorbar.ax.tick_params(labelsize=20, size=0,)
            
        elif cbar == False:
            p = data.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, transform = projection_p, add_colorbar=False, add_labels=False)
                                 
    else:
        p = data.plot.imshow(ax =ax, cmap=cmap, transform = projection_p, 
                                 cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": "horizontal", 
                                               "shrink": 0.70, "format": "%.0f", "ticks":ticks}, extend= "neither")
        
    
    
    # ploting background extent
    
    if plot_coastlines == False:
        sea_land_mask.plot.contour(colors="k", linestyles="-", ax=ax, transform=projection_p, levels=[0], linewidths=3)
        
        
    plot_background(p, domain= domain, left_labels=left_labels, bottom_labels=bottom_labels, 
                    plot_coastlines=plot_coastlines, plot_borders=plot_borders)
    
    if title is not None:
        ax.set_title(title, fontsize=20, weight="bold", loc="left")
        
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
    
def plot_eofsAsCovariance(variable, data, mode_var=None, cmap = None, levels=None, units=None, ax=None, domain=None, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None, cbar = None, cbar_orientation=None, cbar_position = None,
                     fig=None, use_AlberEqualArea=None, vmax=None, vmin=None, left_labels= True, bottom_labels=True):
    """
    

    Parameters
    ----------
    variable : TYPE: str
        DESCRIPTION. The variable to be plotted. Note, it will be display as colorbar name
    data : TYPE: datarray
        DESCRIPTION. The processed data to be visualized (Eg. topo input file or can be retrieved from model output)
    cmap : TYPE: plt.cmap 
        DESCRIPTION. Color map handle from matplotlib
    units : TYPE: str
        DESCRIPTION. The unit of the dataset to be visualized 
    ax : TYPE: matplotlib ax handle, optional
        DESCRIPTION. The default is None.
    fig : TYPE: Matplotlib figure handle, optional
        DESCRIPTION. The default is None.
    vmax : TYPE: float, optional
        DESCRIPTION. The default is None. maximum value limit of the variable to be ploted 
    vmin : TYPE: float, optional
        DESCRIPTION. The default is None. minimum value limit of the variable to be ploted 
    levels : TYPE: float, optional
        DESCRIPTION. The default is None. the number of levels for colorbar scale
    domain : TYPE: str, optional
        DESCRIPTION. The default is None. eg. Africa, Asia, Europe
    output_name : TYPE: str, optional
        DESCRIPTION. The default is None. Filename of generated figures
    output_format : TYPE: str, optional
        DESCRIPTION. The default is None. Format to save figure eg. pdf, svg, tiff
    level_ticks : TYPE: float, optional
        DESCRIPTION. The default is None. Interval of ticks for colorbar
    title : TYPE: str, optional
        DESCRIPTION. The default is None. Title of plots
    path_to_store : TYPE: str, optional
        DESCRIPTION. The default is None. Directory to store data
    cbar : TYPE: Boolean, optional
        DESCRIPTION. The default is None. True is the plot require colobar axis
    cbar_orientation : TYPE: , optional
        DESCRIPTION. 
    cbar_position : TYPE: list, optional
        DESCRIPTION. The default is None. The default is None. the list defing the position of the color bar eg. [0.90, 0.30, 0.02, 0.40]
    mode_var : TYPE: float, optional
        DESCRIPTION. The default is None. The explained variance estimated from the EOF analysis
    use_AlberEqualArea : TYPE: Boolean, optional
        DESCRIPTION. The default is None. To use ccrs.AlberEqualArea() as geoaxis projection
    
    left labels, right_labels, bottom_labels, : TYPE: Bol, optional
        DESCRIPTION. To set the left, right, and bottom axis label to None

    Returns
    -------
    None.

    """
    
    
    norm = MidpointNormalize(midpoint=0)
    projection = ccrs.PlateCarree()
    
    # defining axis
    if ax is None:
        if use_AlberEqualArea == True:
            
            fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":ccrs.AlbersEqualArea(
                central_latitude=40, central_longitude=-35, standard_parallels=(0, 80))})
        else:
             fig, ax = plt.subplots(1, 1, sharex=False, figsize= (15, 13), subplot_kw= {"projection":projection})
             
    if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
        ticks = np.linspace(vmin, vmax, level_ticks)
        if cbar==True:
            
            cbar_pos = cbar_position
            cbar_ax = fig.add_axes(cbar_pos)   # axis for subplot colorbar # left, bottom, width, height
            if cbar_orientation == "vertical":
                cbar_ax.get_xaxis().set_visible(False)
                cbar_ax.yaxis.set_ticks_position('right')
            else:
                cbar_ax.get_yaxis().set_visible(False)
                cbar_ax.xaxis.set_ticks_position('bottom')
                
            cbar_ax.set_yticklabels([])
            cbar_ax.tick_params(size=0)
            
            p = data.plot.imshow(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, transform = projection,
                                 cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": cbar_orientation, 
                                               "shrink": 0.5, "format": "%.0f", "ticks":ticks}, extend= "both", add_colorbar=True, cbar_ax=cbar_ax,
                                 add_labels=False)
            
            
            p.colorbar.set_label(label=variable + " [" + units + "]", size= 22, fontweight="bold")
            p.colorbar.ax.tick_params(labelsize=22, size=0,)
            
        elif cbar == False:
            p = data.plot.imshow(ax =ax, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, transform = projection, add_colorbar=False, add_labels=False)
                                 
    else:
        p = data.plot.imshow(ax =ax, cmap=cmap, transform = projection, 
                                 cbar_kwargs= {"pad":0.1, "drawedges": True, "orientation": "horizontal", 
                                               "shrink": 0.70, "format": "%.0f", "ticks":ticks}, extend= "both")
        
    # ploting background extent
    
    if use_AlberEqualArea ==True: #!!! Its currently not possible to add labels after set boundary, cartopy >0.18 supports labels aside PlateCarree and Mercator but not after clipping boundaries
    
    
    # Alternative is to generate the Map and use any editting software to add the labels
        plot_background(p, domain=domain, use_AlbersEqualArea=True, ax=ax, left_labels=left_labels,
                        bottom_labels=bottom_labels)
    else:
        plot_background(p, domain= domain, bottom_labels=bottom_labels, left_labels=left_labels)
    
    if title is not None:
        if mode_var is not None:
            
            ax.set_title(title + " ({:.2f})".format(mode_var) , fontsize=22, weight="bold", loc="left")
        else:
            ax.set_title(title, fontsize=22, weight="bold", loc="left")
        
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
        
        
def plot_vertical_section(variable, data, cmap, units, season=None, ax=None, fig=None, vmax=None, vmin=None, levels=None, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None, plot_colorbar=True,
                     cbar_pos=None, fig_title=None, season_label=None, geosp_data=None, dim=None, left_labels=False, bottom_labels=False,
                     right_labels=False, use_norm=False, use_cbar_norm=False,
                     data_u =None, data_v=None, plot_winds=False, norm=None):
    
    # extracting coords from data (select up to 200 hPa) [:,:10] y = data.columns.values[:10], data = data.iloc[:, :10]
    x = data.index.values
    y = data.columns.values[:12]
    
    data = data.iloc[:, :12]
    # applying meshgrid 
    X,Y = np.meshgrid(x,y)
    Z = data.values.T
     
    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, sharex= False, sharey= False, figsize=(8, 7))
    
    #ploting with plt.contourf 
    if norm is None:
        norm = MidpointNormalize(midpoint=0)
        
    
    if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
        ticks = np.linspace(vmin, vmax, level_ticks)
        if vmin < 0 and use_norm==True:
            p = ax.contourf(X,Y,Z, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, norm=norm)
        else:
            p = ax.contourf(X,Y,Z, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels,)
    else:
        p = ax.contourf(X,Y,Z, cmap=cmap)
        
    # invert axis bcos of higher values at the surface
    ax.invert_yaxis()
    
    if plot_colorbar == True:
        if use_cbar_norm == True:
            if norm is None:
                norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            
        if cbar_pos is None:
            cbar_pos = [0.90, 0.30, 0.03, 0.45]
        
        if fig is not None:
            
            cbar_ax = fig.add_axes(cbar_pos)   # axis for subplot colorbar # left, bottom, width, height
            
            cbar_ax.get_xaxis().set_visible(False)
            cbar_ax.yaxis.set_ticks_position('right')
            cbar_ax.set_yticklabels([])
            cbar_ax.tick_params(size=0)
        
        if level_ticks is not None:
            if use_cbar_norm == True:
                cb =fig.colorbar(p, cax=cbar_ax, drawedges=True, orientation="vertical", shrink=0.7, 
                         format="%.2f", ticks = ticks, extend = "both", pad = 0.05, norm=norm)
            else:
                cb =fig.colorbar(p, cax=cbar_ax, drawedges=True, orientation="vertical", shrink=0.7, 
                         format="%.2f", ticks = ticks, extend = "both", pad = 0.05,)
        else:
            cb =fig.colorbar(p, cax=cbar_ax, drawedges=True, orientation="vertical", shrink=0.7, 
                     format="%.2f", extend = "neither", pad=0.05, norm=norm)
            
        cb.set_label(label=variable + " [" + units + "]", size= 20, fontweight="bold")
        cb.ax.tick_params(labelsize=20, size=0,)
        
    if geosp_data is not None:
        ax2 = ax.twinx()
        ax2.grid(False)
        
        if season:
            df = geosp_data[season]
        else:
            df = geosp_data
            
        ax2.fill_between(df.index, df/1000, 0, color=black, edgecolor=black, linestyle="-", linewidth=2.5)
        
        # setting limit to match pressure levels (Try with the cdo converted height levels later)
        ax2.set_ylim(0, 11.5)
    
    
    if plot_winds == True:
        
        if all(data is not None for data in [data_v, data_u]):
            # extracting variables for quiver 
            x = data_v.index.values[::1]
            y = data_v.columns.values[:14]
        
            u = data_u.iloc[::1, :14] # meridoinal values (m/s)
            v = data_v.iloc[::1, :14] * -120 #(omega in Pa/s)
        
            X,Y = np.meshgrid(x,y)
            #skip = (slice(None, None, 3), slice(None, None, 3))  #for extracting the data on interval or use data[::3, ::3]
            
            # ploting winds using quiver 
            q = ax.quiver(X, Y, u.values, v.values,  pivot= "mid", scale= 100,
                          headwidth=3, headlength=5, headaxislength=4.5)
            
            qk = ax.quiverkey(q, 0.90, -0.1, 2, r'$2 \frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties=
                              {"size": 20, "weight":"bold"})
        
    if bottom_labels == True:
        if dim == "lon":
             ax.set_xlabel("Longitude [E°]", fontsize=22, fontweight="bold")
        elif dim == "lat":
             ax.set_xlabel("Latitude [N°]", fontsize=22, fontweight="bold")
        else:
            raise ValueError("Define dim as lat or lon")

        
            
    if left_labels == True:
        ax.set_ylabel("Pressure [hPa]", fontsize=22, fontweight="bold")
            
    
    if left_labels ==False:
        ax.grid(False)
        #ax.axes.yaxis.set_visible(False)
        ax.set_yticklabels([])
    
        
    if bottom_labels == False:
        ax.grid(False)
        #ax.axes.yaxis.set_visible(False)
        ax.set_xticklabels([])
        
    if right_labels == True:
    
        ax2.set_ylabel("Height [km]", fontsize=22, fontweight="bold")
        
    else:
        if geosp_data is not None:
            ax2.set_yticklabels([])
    
    
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="left")
        
    plt.tight_layout()
    
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
        
               
def plot_hovmoller_space_time(variable, data, cmap, units, ax=None, fig=None, vmax=None, vmin=None, levels=None, output_name=None, 
                     output_format=None, level_ticks=None, title=None, path_to_store=None, plot_colorbar=True,
                     cbar_pos=None, fig_title=None, left_labels=False, bottom_labels=False,
                     right_labels=False, use_norm=False, use_cbar_norm=False,
                     plot_contour=True):
    
    # extracting coords from data (select up to 200 hPa) [:,:10] y = data.columns.values[:10], data = data.iloc[:, :10]
    y = data.index.values
    x = data.columns.values
    
    # applying meshgrid 
    X,Y = np.meshgrid(x,y)
    Z = data.values
     
    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, sharex= False, sharey= False, figsize=(8, 7))
    
    #ploting with plt.contourf 
    norm = MidpointNormalize(midpoint=0)
    
    if all(parameter is not None for parameter in [vmin, vmax, levels, level_ticks]):
        ticks = np.linspace(vmin, vmax, level_ticks)
        if vmin < 0 and use_norm==True:
            p = ax.contourf(X,Y,Z, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels, norm=norm)
        else:
            p = ax.contourf(X,Y,Z, cmap=cmap, vmin=vmin, vmax=vmax, levels=levels,)
    else:
        p = ax.contourf(X,Y,Z, cmap=cmap)
    
    
    if plot_colorbar == True:
        if use_cbar_norm == True:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            
        if cbar_pos is None:
            cbar_pos = [0.90, 0.30, 0.03, 0.45]
        
        if fig is not None:
            
            cbar_ax = fig.add_axes(cbar_pos)   # axis for subplot colorbar # left, bottom, width, height
            
            cbar_ax.get_xaxis().set_visible(False)
            cbar_ax.yaxis.set_ticks_position('right')
            cbar_ax.set_yticklabels([])
            cbar_ax.tick_params(size=0)
        
        if level_ticks is not None:
            if use_cbar_norm == True:
                cb =fig.colorbar(p, cax=cbar_ax, drawedges=True, orientation="vertical", shrink=0.7, 
                         format="%.2f", ticks = ticks, extend = "neither", pad = 0.05, norm=norm)
            else:
                cb =fig.colorbar(p, cax=cbar_ax, drawedges=True, orientation="vertical", shrink=0.7, 
                         format="%.2f", ticks = ticks, extend = "neither", pad = 0.05,)
        else:
            cb =fig.colorbar(p, cax=cbar_ax, drawedges=True, orientation="vertical", shrink=0.7, 
                     format="%.2f", extend = "neither", pad=0.05, norm=norm)
            
        cb.set_label(label=variable + " [" + units + "]", size= 20, fontweight="bold")
        cb.ax.tick_params(labelsize=20, size=0,)
    
    if plot_contour == True:
        
        c = ax.contour(X,Y,Z, color="black", linewidth=2, levels=10)
        clb = ax.clabel(c, fmt="%2.0f", use_clabeltext=True, colors="black", fontsize=22)
    if bottom_labels == True:
        
        ax.set_xlabel("Months ", fontsize=22, fontweight="bold")
        ax.xaxis.set_ticks(np.arange(min(x), max(x)+1, 1.0))

        
            
    if left_labels == True:
        ax.set_ylabel("Latitude [°N]", fontsize=22, fontweight="bold")
            
    
    if left_labels ==False:
        ax.grid(False)
        #ax.axes.yaxis.set_visible(False)
        ax.set_yticklabels([])
    
        
    if bottom_labels == False:
        ax.grid(False)
        #ax.axes.yaxis.set_visible(False)
        ax.set_xticklabels([])
    
    
    
    if title is not None:
        ax.set_title(title, fontdict= {"fontsize": 22, "fontweight":"bold"}, loc="left")
        
    plt.tight_layout()
    
    #optional if one plot is required, alternatively save from the control script
    if all(parameter is not None for parameter in [output_format, output_name, path_to_store]):
        plt.savefig(os.path.join(path_to_store, output_name + "." + output_format), format= output_format, bbox_inches="tight")
        
        
        
def plot_wind_streamlines():
    pass

def plot_vertical_winds(*args):
    pass



