# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:47:02 2023

@author: dboateng

Contain all the plotting functions (time series plots, box plots, covariance plots, etc.)
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath





proj = ccrs.AlbersEqualArea(central_longitude=-25,
                            central_latitude=40,
                            standard_parallels=(40, 80))
ax = plt.axes(projection=proj)    
ax.set_extent([-80, 30, 20, 80], crs=ccrs.PlateCarree())
ax.coastlines()

# Make a boundary path in PlateCarree projection, I choose to start in
# the bottom left and go round anticlockwise, creating a boundary point
# every 1 degree so that the result is smooth:
vertices = [(lon, 40) for lon in range(-80, 41, 1)] + \
            [(lon, 80) for lon in range(40, -81, -1)]
boundary = mpath.Path(vertices)
ax.set_boundary(boundary, transform=ccrs.PlateCarree())

plt.show()