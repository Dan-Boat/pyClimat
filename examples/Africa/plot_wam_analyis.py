#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:30:05 2022

@author: dboateng
"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean

from wam_analysis import *


def plot_monthly_sections(data_sahara, data_sahel, data_guinea, ax, ymax, ymin, varname, title):
    
    ax.plot(data_sahara["mean"], "--", color=red, label="Sahara", linewidth=3)
    ax.plot(data_sahel["mean"], "--", color=black, label="Sahel", linewidth=3)
    ax.plot(data_guinea["mean"], "--", color=blue, label="Coast of Guinea", linewidth=3)
    
    ax.set_ylabel(varname, fontsize= 20, fontweight="bold")
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.set_title(title, fontdict= {"fontsize": 15, "fontweight":"bold"}, loc= "left")
    
    
    
# fonts and ploting stlye 
apply_style(fontsize=22, style=None, linewidth=2) 

projection = ccrs.PlateCarree()
fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(20, 13), subplot_kw={"projection":
                                                                                                                  projection})

# add patch (Sahel: 10-20 N, 20W - 30E; coast of Guinea: 5-10N, 20W-30E, Sahara region: 20-30N, 20W-30E)

#sahara --> 30, 20 N 20W, 30 E

lat_sh, h_sh = 20 , 10  # lat and height
lon_sh , w_sh = -20, 50  # long and width 

#sahel --> 20, 10 N 20W, 30 E

lat_sl, h_sl = 10 , 10  # lat and height
lon_sl , w_sl = -20, 50  # long and width 

#guinea -->5, 10 N 20W, 30 E

lat_g, h_g = 5 , 5  # lat and height
lon_g , w_g = -20, 50  # long and width 


plot_annual_mean(ax=ax1, variable="Precipitation", data_alt=PI_prec_alt, cmap=YlGnBu, units="mm/month", vmax=350, vmin=50, domain="West Africa", 
                  levels=22, level_ticks=6, title="[A]", left_labels=True, bottom_labels=False, use_colorbar_default=True)

ax1.add_patch(patches.Rectangle(xy =(lon_w, lat_w), width= w_w, height=h_w, ls= "--", color= red, transform = projection, 
                                fc="None", lw=2.5,))


plot_annual_mean(ax=ax2, variable="Temperature", data_alt=PI_t2m_alt, cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                  levels=22, level_ticks=11, title="[B]", left_labels=True, bottom_labels=True, use_colorbar_default=True)

fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(20, 13))

plot_monthly_sections(data_sahara=PI_month_sahara_prec, data_sahel=PI_month_sahel_prec, data_guinea=PI_month_guinea_prec,
                      ax=ax1, ymax=150, ymin=0, varname="Precipitation [mm/month]", title="[A]")

plot_monthly_sections(data_sahara=PI_month_sahara_t2m, data_sahel=PI_month_sahel_t2m, data_guinea=PI_month_guinea_t2m,
                      ax=ax2, ymax=30, ymin=5, varname="Temperature [°C]", title="[B]")

ax2.legend(bbox_to_anchor=(1.04,1),frameon=True, fontsize=25, loc="upper left")
plt.tight_layout()

plt.show()