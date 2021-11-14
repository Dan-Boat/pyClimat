#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 12:47:37 2021


@author: dboateng
"""
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *



module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
years= "1003_1017"
period = "1m"
data = read_ECHAM_processed(main_path= module_output_main_path , exp_name=exp_name, read_wiso=False, 
                            period="1m", add_name="geosp")

data_slp = read_ECHAM_processed(main_path= module_output_main_path , exp_name=exp_name, read_wiso=False, 
                            period="1m", add_name="slp")

# extracting geop500

geopoth = extract_var(Dataset=data , varname="geopoth", lev_units="hPa")

slp = extract_var(Dataset=data_slp, varname="slp", units="mbar")

eofs = EOF_analysis(data=geopoth, maxlon=40, minlon=-80, maxlat=80, minlat=20, season="DJF",
                    neofs=3, pcscaling=1, apply_coslat_weights=True, lev=500)

eofs_slp, pcs, var_frac = EOF_analysis(data=slp, maxlon=40, minlon=-80, maxlat=80, minlat=20, season="DJF",
                    neofs=3, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=3,
                    return_pcs=True, npcs=3)

NAO, mode_var1 = (eofs_slp[0]*-1), var_frac[0]

EA, mode_var2 = (eofs_slp[1]*-1), var_frac[1]

SCA, mode_var3 = (eofs_slp[2]), var_frac[2]

projection = ccrs.AlbersEqualArea(
                central_latitude=40, central_longitude=-35, standard_parallels=(0, 80))

fig, (ax1,ax2, ax3) = plt.subplots(nrows = 3, ncols = 1, figsize=(15, 17), subplot_kw={"projection":
                                                                                                                  projection})
plot_eofsAsCovariance(variable= "slp", data=NAO, mode_var=mode_var1, units="mbar", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.90, 0.30, 0.02, 0.40], cbar_orientation="vertical",use_AlberEqualArea=(True), ax=ax1, fig=fig, title="NAO")

plot_eofsAsCovariance(variable= "slp", data=EA, mode_var=mode_var2, units="mbar", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=(True), ax=ax2, fig=fig, title="EA")

plot_eofsAsCovariance(variable= "slp", data=SCA, mode_var=mode_var3, units="mbar", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=(True), ax=ax3, fig=fig, title="SCA")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.06)


# Write this script into the example folder for all the experiments...save in pdf to add coordinates before adding to manuscript.

print(eofs)

# # vertical levels in Pa --> hPa

# geopoth["lev"] = geopoth.lev / 100
# geopoth["lev"].attrs["units"] = "hPa"

# #selecting pressure level of 500 hPa

# geopoth500 = geopoth.sel(lev=500)

# # grouping into season

# geopoth500 = geopoth500.groupby("time.season")

# # extracting winter 

# geopoth500_DJF = geopoth500["DJF"]

# # extracting transect for EOF

# maxlat = 80
# minlat= 20
# maxlon= 40
# minlon=-80

# #convert lon to -180 to 180
# geopoth500_DJF = geopoth500_DJF.assign_coords({"lon": (((geopoth500_DJF.lon + 180) % 360) - 180)})

# geopoth500_NH = geopoth500_DJF.where((geopoth500_DJF.lat >= minlat)&(geopoth500_DJF.lat <=maxlat), drop=True)
# geopoth500_NH = geopoth500_NH.where((geopoth500_NH.lon >= minlon)&(geopoth500_NH.lon <=maxlon), drop=True)


# # anomalies 

# data_mean = geopoth500_NH.mean(dim="time")

# anomalies = geopoth500_NH - data_mean

# #standardize with the std

# data_std = anomalies.std(dim="time")

# anomalies /= data_std

# wtgs = np.sqrt(np.abs(np.cos(anomalies.lat*np.pi/180)))

# anomalies_wtgs = anomalies * wtgs
# # loading patterns to eof 

# solver = Eof(anomalies_wtgs)

# #solver = Eof(anomalies, weights=wtgs)

# eofs = solver.eofsAsCovariance(neofs=3, pcscaling=1)
# eofs = eofs.sortby(eofs.lon)

# var_frac = solver.varianceFraction(neigs=3)
# print(var_frac)

# pcs = solver.pcs(pcscaling=1, npcs=3)

# # visualise
# projection = ccrs.PlateCarree()
# ax = plt.axes(projection=projection)
# #ax = plt.axes(projection= ccrs.Orthographic(central_longitude=10.0, central_latitude=30.0,))
# p=eofs[0].plot.imshow(ax=ax, transform=projection)

# plot_background(p, domain= "Europe")
# # gl= p.axes.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1,
# #                      edgecolor = "gray", linestyle = "--", color="gray", alpha=0.5)
    
# # gl.top_labels = False                  # labesl at top
# # gl.right_labels = False
# # gl.xformatter = LongitudeFormatter()     # axis formatter
# # gl.yformatter = LatitudeFormatter()
# # gl.xlabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}   #axis style 
# # gl.ylabel_style = {"fontsize": 20, "color": "black", "fontweight": "semibold"}
# print(geopoth)

