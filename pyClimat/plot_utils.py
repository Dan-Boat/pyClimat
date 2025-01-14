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
from matplotlib import rc
import matplotlib.colors


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
cyan = "cyan"
olive = "olive"

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
PuOr_r = plt.cm.PuOr_r
BrBG = plt.cm.BrBG
BrBG_r = plt.cm.BrBG_r
RdYlBu = plt.cm.RdYlBu
RdYlBu_r = plt.cm.RdYlBu_r
RdYlGn = plt.cm.RdYlGn
Spectral = plt.cm.Spectral
Spectral_r = plt.cm.Spectral_r
YlGnBu = plt.cm.YlGnBu
winter = plt.cm.winter
PuBu = plt.cm.PuBu
PuBu_r = plt.cm.PuBu_r
plasma = plt.cm.plasma
hot = plt.cm.hot_r

#customized plots
GryBr = matplotlib.colors.LinearSegmentedColormap.from_list('my_gradient', (
    (0.000, (0.000, 0.000, 0.000)),
    (0.250, (0.533, 0.533, 0.533)),
    (0.500, (1.000, 1.000, 1.000)),
    (0.750, (0.886, 0.537, 0.137)),
    (1.000, (0.412, 0.129, 0.000))))


GryBr_r = matplotlib.colors.LinearSegmentedColormap.from_list('my_gradient', (
    (0.000, (0.412, 0.129, 0.000)),
    (0.250, (0.886, 0.537, 0.137)),
    (0.500, (1.000, 1.000, 1.000)),
    (0.750, (0.533, 0.533, 0.533)),
    (1.000, (0.000, 0.000, 0.000))))


PrecAno = matplotlib.colors.LinearSegmentedColormap.from_list('precipitation_anomalies', (
    # Edit this gradient at https://eltos.github.io/gradient/#precipitation_anomalies=6F3200-FF953F-FFFFFF-009B30-00551A
    (0.000, (0.435, 0.196, 0.000)),
    (0.250, (1.000, 0.584, 0.247)),
    (0.500, (1.000, 1.000, 1.000)),
    (0.750, (0.000, 0.608, 0.188)),
    (1.000, (0.000, 0.333, 0.102))))

corr = matplotlib.colors.LinearSegmentedColormap.from_list('correlation', (
    # Edit this gradient at https://eltos.github.io/gradient/#correlation=0000D6-6F6FFF-FCFCFF-D7A959-815508
    (0.000, (0.000, 0.000, 0.839)),
    (0.250, (0.435, 0.435, 1.000)),
    (0.500, (0.988, 0.988, 1.000)),
    (0.750, (0.843, 0.663, 0.349)),
    (1.000, (0.506, 0.333, 0.031))))


vegetation = matplotlib.colors.LinearSegmentedColormap.from_list('correlation', (
    # Edit this gradient at https://eltos.github.io/gradient/#correlation=0.1:FFFEFE-0.1:795519-30:D7A959-53.2:F6D743-77.7:2DD62D-100:125912
    (0.000, (1.000, 0.996, 0.996)),
    (0.001, (1.000, 0.996, 0.996)),
    (0.001, (0.475, 0.333, 0.098)),
    (0.300, (0.843, 0.663, 0.349)),
    (0.532, (0.965, 0.843, 0.263)),
    (0.777, (0.176, 0.839, 0.176)),
    (1.000, (0.071, 0.349, 0.071))))

causal = matplotlib.colors.LinearSegmentedColormap.from_list('my_gradient', (
    # Edit this gradient at https://eltos.github.io/gradient/#EEFBFA-00FFEB
    (0.000, (0.933, 0.984, 0.980)),
    (1.000, (0.000, 1.000, 0.922))))


PrecipitationCal = matplotlib.colors.LinearSegmentedColormap.from_list('my_gradient', (
    # Edit this gradient at https://eltos.github.io/gradient/#FEFEFE-9394FF-0005FF-7EFF99-07FF00-FFFC74-FFF800-EA5E5A-FB0604
    (0.000, (0.996, 0.996, 0.996)),
    (0.125, (0.576, 0.580, 1.000)),
    (0.250, (0.000, 0.020, 1.000)),
    (0.375, (0.494, 1.000, 0.600)),
    (0.500, (0.027, 1.000, 0.000)),
    (0.625, (1.000, 0.988, 0.455)),
    (0.750, (1.000, 0.973, 0.000)),
    (0.875, (0.918, 0.369, 0.353)),
    (1.000, (0.984, 0.024, 0.016))))
# defining plot styles (which contains fonts and backgrouds)
# plt.style.use (can be seaborn, dark_background, fivethirtyeight, bmh



def apply_style(fontsize=28, style=None, linewidth=2):
    """
    Add more attributes: https://matplotlib.org/stable/tutorials/introductory/customizing.html
    Parameters
    ----------
    fontsize : TYPE, optional
        DESCRIPTION. The default is 10.
    style : TYPE, optional
        DESCRIPTION. The default is "bmh". ["seaborn", "fivethirtyeight",]
    Returns
    -------
    None.
    """
    if style is not None:
        plt.style.use(style)  #seaborn-talk 'seaborn-paper' 'seaborn-poster'
        
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    mpl.rc('lines', linewidth=linewidth)
    mpl.rc("font", weight="bold")
    
    mpl.rc('text', usetex=True)
    mpl.rc('xtick', labelsize=fontsize)
    mpl.rc('ytick', labelsize=fontsize)
    mpl.rc('ytick.major', size=6)
    mpl.rc('xtick.major', size=6)
    mpl.rc('legend', fontsize=fontsize)
    mpl.rc('axes', labelsize=fontsize)
    mpl.rc('axes', titlesize=fontsize)
    mpl.rc('axes', titleweight="bold")
    mpl.rc('axes', labelweight="bold")
    
    # mpl.rc('figure', labelsize=fontsize)
    # mpl.rc('figure', titlesize=fontsize)
    # mpl.rc('figure', titleweight="bold")
    # mpl.rc('figure', labelweight="bold")
    
def apply_style2(fontsize=20, style=None, linewidth=2, usetex=True):
    """
    

    Parameters
    ----------
    fontsize : TYPE, optional
        DESCRIPTION. The default is 10.
    style : TYPE, optional
        DESCRIPTION. The default is "bmh". ["seaborn", "fivethirtyeight",]

    Returns
    -------
    None.

    """
    if style is not None:
        plt.style.use(style)  
        
    rc('font',**{'family':'sans-serif','sans-serif':['Arial'], 
                 'size': 22, 'style': 'normal', 'weight': 'medium'})
    
    mpl.rc('text', usetex=usetex)
    mpl.rc('xtick', labelsize=fontsize)
    mpl.rc('ytick', labelsize=fontsize)
    mpl.rc('legend', fontsize=fontsize)
    mpl.rc('axes', labelsize=fontsize)
    mpl.rc('lines', linewidth=linewidth)
    mpl.rc("font", weight="bold")
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rcParams['svg.fonttype'] = 'none'
    

# defining function for selecting background domain for cartopy

def plot_background(p, domain=None, use_AlbersEqualArea=False,ax=None, left_labels=True,
                    bottom_labels=True, plot_coastlines=True, plot_borders=False, coast_resolution=None,
                    coast_color="black"):
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
    
    if plot_coastlines ==True:
        if coast_resolution is None:
            p.axes.coastlines(resolution = "50m", linewidth=1.5, color=coast_color)  # add coastlines outlines to the current axis (110m, 50m , 10m)
            
        else:
            p.axes.coastlines(resolution = coast_resolution, linewidth=1.5, color=coast_color)  # add coastlines outlines to the current axis
    
    if plot_borders == True:
        p.axes.add_feature(cfeature.BORDERS, edgecolor="grey", linewidth = 0.5) #adding country boarder lines
    
    #setting domain size
    if domain is not None: 
        if domain == "Europe":   # Europe
            minLon = -20
            maxLon = 35
            minLat = 35
            maxLat = 65
        elif domain == "South America":   # South America
            minLon = -83
            maxLon = -32
            minLat = -35
            maxLat = 5
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
            
        elif domain == "GNIP_view":        # Eurasia
            minLon = -25
            maxLon = 30
            minLat = 40
            maxLat = 72
            
            
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
            
        elif domain == "East Africa":  # Africa
            minLon = 20
            maxLon = 55
            minLat = -20
            maxLat = 18
            
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
            maxLon = 30
            minLat = 20
            maxLat = 80
            
        elif domain == "West Africa":
            minLon = -25
            maxLon = 40
            minLat = -5
            maxLat = 35
            
        elif domain == "West Africa Wide":
            minLon = -45
            maxLon = 40
            minLat = 0
            maxLat = 45
            
        elif domain == "NH Wide":
            minLon = -120
            maxLon = 120
            minLat = 25
            maxLat = 85
            
        elif domain == "Europe Wide":
            minLon = -35
            maxLon = 35
            minLat = 32
            maxLat = 75
            
        elif domain =="West Africa ANDEL":
            minLon = -20
            maxLon = 10
            minLat = 0
            maxLat = 35
        
        elif domain =="Alps":
            minLon = 3
            maxLon = 15
            minLat = 43.5
            maxLat = 49.5
        
            
            
        
            
        else:
            print("ERROR: invalid geographical domain passed in options")
        p.axes.set_extent([minLon, maxLon, minLat, maxLat], ccrs.PlateCarree())
    if domain is None: 
        #uncomment 
        #p.axes.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
        print(None)
        
    if use_AlbersEqualArea == True:
        vertices = [(lon, 20) for lon in  range(-80, 21, 1)]
        vertices += [(lon, 80) for lon in range(20, -81,-1)]
       
        
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
    
    #uncoment
    
    gl.xlabel_style = {"fontsize": 25, "color": "black", "fontweight": "semibold"}   #axis style 
    gl.ylabel_style = {"fontsize": 25, "color": "black", "fontweight": "semibold"}

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


class FixedPointNormalized(matplotlib.colors.Normalize):
    
    def __init__(self, vmin=None, vmax=None, sealevel=0, color_val=0.21875,
                 clip=False):
        
        self.sealevel = sealevel
        self.color_val = color_val
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)
        
    def __call__(self, value, clip=None):
        
        x,y = [self.vmin, self.sealevel, self.vmax], [0, self.color_val, 1]
        
        return np.ma.masked_array(np.interp(value, x, y))



def creat_norm():
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    import matplotlib.colors as col
    
    levels = [i for i in range(-100, 4000, 100)]
    terrain_new = mpl.cm.get_cmap("terrain", 256)
    terrain_adjust = ListedColormap(terrain_new(np.linspace(0.23, 1, 256)))
    new_colors = terrain_adjust(np.linspace(0,1,256))
    
    blue = np.array([135/256, 206/256, 250/256, 1])
    new_colors[:1, :] = blue
    terrain_shift = ListedColormap(new_colors)

    norm_new = col.BoundaryNorm(levels, ncolors=terrain_shift.N, clip=True)
    
    return norm_new, terrain_shift


