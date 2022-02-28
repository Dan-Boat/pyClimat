#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 18:50:11 2021

@author: dboateng
This module contains all the utilities used in the Climat_plots
"""
# Import modules

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib.path as mpath
from cycler import cycler
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
from cartopy.util import add_cyclic_point

# defining image sizes and colors 
# A4 paper size: 210 mm X 297 mm

cm = 0.3937   # 1 cm in inch for plot size
pt = 1/72.27  # pt in inch from latex geometry package
textwidth = 345*pt
big_width = textwidth + 2*3*cm

# colors
orange = 'orangered'
lightblue = 'teal'
brown = 'sienna'
red = '#a41a36'
blue = '#006c9e'
green = '#55a868'
purple = '#8172b2'
lightbrown = '#ccb974'
pink = 'fuchsia'
lightgreen = 'lightgreen'
skyblue = "skyblue"
tomato = "tomato"
gold = "gold"
magenta = "magenta"
black = "black"
grey = "grey"
golden = "darkgoldenrod"

#divergence colors 
RdBu_r = plt.cm.RdBu_r
RdBu = plt.cm.RdBu
Blues = plt.cm.Blues
Spectral_r = plt.cm.Spectral_r
terrain = plt.cm.terrain
viridis = plt.cm.viridis
BwR_r = plt.cm.bwr_r
BwR = plt.cm.bwr
Seismic = plt.cm.seismic
Reds = plt.cm.Reds
Greys = plt.cm.Greys
RdGy = plt.cm.RdGy
PiYG = plt.cm.PiYG_r
PuOr = plt.cm.PuOr
BrBG = plt.cm.BrBG
RdYlBu = plt.cm.RdYlBu
RdYlBu_r = plt.cm.RdYlBu_r
RdYlGn = plt.cm.RdYlGn
Spectral = plt.cm.Spectral
Spectral_r = plt.cm.Spectral_r
YlGnBu = plt.cm.YlGnBu
winter = plt.cm.winter
PuBu = plt.cm.PuBu
PuBu_r = plt.cm.PuBu_r

# defining plot styles (which contains fonts and backgrouds)
# plt.style.use (can be seaborn, dark_background, fivethirtyeight, bmh)
def apply_style(fontsize=21):
    font = {'family':'serif','size':20, 'serif': ['computer modern roman']}
    plt.rc('font',**font)
    plt.rc('legend',**{'fontsize':20})
    small_size=fontsize-1
    
    plt.style.use("seaborn-colorblind")    
    plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
    mpl.rc('text', usetex=True)
    #mpl.rc('font', size=20, family='serif')
    mpl.rc('xtick', labelsize=small_size)
    mpl.rc('ytick', labelsize=small_size)
    #mpl.rc('legend', fontsize=small_size)
    mpl.rc('axes', labelsize=fontsize)
    mpl.rc('lines', linewidth=2.5)
    #mpl.rc("font", weight="bold")
   

def apply_style_2(fontsize=21):
    small_size=fontsize-1
    plt.style.use('fivethirtyeight')  
    plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
    # plt.rcParams.update({"text.usetex": True, "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"]})
    mpl.rc('text', usetex=True)
    mpl.rc('font', size=fontsize, family='serif')
    mpl.rc('xtick', labelsize=small_size)
    mpl.rc('ytick', labelsize=small_size)
    mpl.rc('legend', fontsize=small_size)
    mpl.rc('axes', labelsize=fontsize)
    mpl.rc('lines', linewidth=3)

# defining function for selecting background domain for cartopy

def plot_background(p, domain=None, use_AlbersEqualArea=None,ax=None, left_labels=True,
                    bottom_labels=True):
    """
    This funtion defines the plotting domain and also specifies the background. It requires 
    the plot handle from xarray.plot.imshow and other optional arguments 
    Parameters
    -------------

    p: TYPE: plot handle 
    DESCRIPTION: the plot handle after plotting with xarray.plot.imshow
    
    domian = TYPE:str 
    DESCRIPTION: defines the domain size, eg. "Europe", "Asia", "Africa"
                  "South America", "Alaska", "Tibet Plateau" or "Himalaya", "Eurosia",
                  "New Zealand", default: global
    """
    p.axes.set_global()                    # setting global axis 
    p.axes.coastlines(resolution = "50m")  # add coastlines outlines to the current axis
    p.axes.add_feature(cfeature.BORDERS, edgecolor="black", linewidth = 0.3) #adding country boarder lines
    
    #setting domain size
    if domain is not None: 
        if domain == "Europe":   # Europe
            minLon = -15
            maxLon = 40
            minLat = 35
            maxLat = 65
        elif domain == "South America":   # South America
            minLon = -83
            maxLon = -32
            minLat = -25
            maxLat = 15
        elif domain == "Tibetan Plateau" or domain == "Himalayas":  #Tibet Plateau/Himalayas
            minLon = 40
            maxLon = 120
            minLat = 0
            maxLat = 60
        elif domain == "Eurasia":        # Eurasia
            minLon = -18
            maxLon = 164
            minLat = 20
            maxLat = 77
        elif domain == "Cascades":  # Cascades
            minLon = -129
            maxLon = -120
            minLat = 45
            maxLat = 52
        elif domain == "Alaska":  # Alaska
            minLon = -165
            maxLon = -125
            minLat = 52
            maxLat = 68
        elif domain == "Africa":  # Africa
            minLon = -30
            maxLon = 55
            minLat = -35
            maxLat = 40
        elif domain == "New Zealand":  # New Zealand 
            minLon = 165
            maxLon = 180
            minLat = -47
            maxLat = -34
        elif domain == "Olympic Mnts":  # Olympic Mnt's
            minLon = -126
            maxLon = -118
            minLat = 43
            maxLat = 52
        elif domain == "NH": # Northen Hemisphere
            minLon = -80
            maxLon = 60
            minLat = 20
            maxLat = 80
            
        else:
            print("ERROR: invalid geographical domain passed in options")
        p.axes.set_extent([minLon, maxLon, minLat, maxLat], ccrs.PlateCarree())
    if domain is None: 
        p.axes.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        
    if use_AlbersEqualArea == True:
        vertices = [(lon, 0) for lon in  range(-100, 31, 1)]
        vertices += [(lon, 90) for lon in range(30, -101,-1)]
        boundary= mpath.Path(vertices)
        p.axes.set_boundary(boundary, transform=ccrs.PlateCarree())
        #ax.set_global()                    # setting global axi
        gl= ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1,
                     edgecolor = "gray", linestyle = "--", color="gray", alpha=0.5)
    
    else:
        
        # adding gridlines    
        gl= p.axes.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1,
                     edgecolor = "gray", linestyle = "--", color="gray", alpha=0.5)
    
    gl.top_labels = False                  # labesl at top
    gl.right_labels = False
    
    if left_labels == True:
        gl.left_labels = True
    else:
        gl.left_labels = False
    
    if bottom_labels == True:
        gl.bottom_labels =True
    else:
        gl.bottom_labels = False
        
    gl.xformatter = LongitudeFormatter()     # axis formatter
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}   #axis style 
    gl.ylabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}

# for divergence plots
class MidpointNormalize(colors.Normalize):
    """
    At the moment its a bug to use divergence colormap and set the colorbar range midpoint 
    to zero if both vmax and vmin has different magnitude. This might be possible in 
    future development in matplotlib through colors.offsetNorm(). This class was original developed 
    by Joe Kingto and modified by Daniel Boateng. It sets the divergence color bar to a scale of 0-1 by dividing the midpoint to 0.5
    Use this class at your own risk since its non-standard practice for quantitative data.
    """
    
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))      
    