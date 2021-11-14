#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 09:38:03 2021

@author: dboateng
"""
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

from eofs.xarray import Eof
import matplotlib.path as mpath


module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
years= "1003_1017"
period = "1m"
data = read_ECHAM_processed(main_path= module_output_main_path , exp_name=exp_name, read_wiso=False, 
                            period="1m", add_name="slp")

# extracting slp

slp = extract_var(Dataset=data , varname="slp")


# Pa --> mbar

slp = slp / 100
slp.attrs["units"] = "mbar"

slp_season = slp.groupby("time.season")
slp_winter = slp_season["DJF"]


# domian for EOF
maxlat = 80
minlat= 20
maxlon= 40
minlon=-80

# extracting transects 

#convert lon to -180 to 180
slp_winter = slp_winter.assign_coords({"lon": (((slp_winter.lon + 180) % 360) - 180)})

slp_winter = slp_winter.where((slp_winter.lat >= minlat)&(slp_winter.lat <=maxlat), drop=True)
slp_winter= slp_winter.where((slp_winter.lon >= minlon)&(slp_winter.lon <=maxlon), drop=True)

# estimating anomalies 
slp_mean = slp_winter.mean(dim="time")
slp_ano = slp_winter - slp_mean

# standardized anomalies 

# slp_ano_std = slp_ano.std(dim="time")

# slp_ano /= slp_ano_std

wtgs = np.sqrt(np.abs(np.cos(slp_ano.lat*np.pi/180)))

anomalies_wtgs = slp_ano * wtgs
# loading patterns to eof 

solver = Eof(anomalies_wtgs)

#solver = Eof(anomalies, weights=wtgs)

eofs = solver.eofsAsCovariance(neofs=3, pcscaling=1)
eofs = eofs.sortby(eofs.lon)

var_frac = solver.varianceFraction(neigs=3)
print(var_frac)

pcs = solver.pcs(pcscaling=1, npcs=3)

# visualise
projection = ccrs.AlbersEqualArea()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize= (13, 13), subplot_kw={"projection":ccrs.AlbersEqualArea(central_latitude=40, central_longitude=-35)})
#ax = plt.axes(projection=ccrs.PlateCarree())
#ax = plt.axes(projection= ccrs.Orthographic(central_longitude=10.0, central_latitude=30.0,))
p=eofs[0].plot.imshow(ax=ax, transform=ccrs.PlateCarree(), vmax=10, vmin=-10, levels = 22, center=0, cmap=RdBu_r)

ax.set_extent([-100, 30, 0, 90], crs=ccrs.PlateCarree())
#ax.coastlines()
# vertices = [(lon, 0) for lon in  range(-100, 31, 1)]
# vertices += [(lon, 90) for lon in range(30, -101,-1)]
# boundary= mpath.Path(vertices)
# ax.set_boundary(boundary, transform=ccrs.PlateCarree())
# ax.set_global()                    # setting global axis 
ax.coastlines(resolution = "50m")  # add coastlines outlines to the current axis
ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth = 0.3) #adding country boarder lines
gl=ax.gridlines(draw_labels = True, linewidth = 1,
                     edgecolor = "gray", linestyle = "--", color="gray", alpha=0.5, dms=True)

gl.top_labels = False                  # labesl at top
gl.right_labels = False
gl.xformatter = LongitudeFormatter()     # axis formatter
gl.yformatter = LatitudeFormatter()
gl.xlabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}   #axis style 
gl.ylabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}

#plot_background(p, domain= "NH")

plt.show()



