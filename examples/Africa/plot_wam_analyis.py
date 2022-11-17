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
from pyClimat.plots import plot_annual_mean, plot_seasonal_mean

from wam_analysis import *

path_to_store = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"


def plot_monthly_sections(data_sahara, data_sahel, data_guinea, ax, ymax, ymin, varname, title):
    
    ax.plot(data_sahara["mean"], "--", color=red, label="Sahara", linewidth=3)
    ax.plot(data_sahel["mean"], "--", color=black, label="Sahel", linewidth=3)
    ax.plot(data_guinea["mean"], "--", color=blue, label="Coast of Guinea", linewidth=3)
    
    ax.set_ylabel(varname, fontsize= 20, fontweight="bold")
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.set_title(title, fontdict= {"fontsize": 15, "fontweight":"bold"}, loc= "left")
    
    
def plot_monthly_period_per_section(PI, MH, LGM, PLIO, ax, title): # extend with ERA (gold), SSP2.6(purple), 4,5(brown), 8.5(orange)  
    
    ax.plot(PI["mean"], "--", color=black, label="PI", linewidth=3)
    ax.plot(MH["mean"], "--", color=red, label="MH", linewidth=3)
    ax.plot(LGM["mean"], "--", color=blue, label="LGM", linewidth=3)
    ax.plot(PLIO["mean"], "--", color=green, label="PLIO", linewidth=3)
    
    ax.set_ylabel(varname, fontsize= 20, fontweight="bold")
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.set_title(title, fontdict= {"fontsize": 15, "fontweight":"bold"}, loc= "left")
      
# fonts and ploting stlye 

def plot_PI_JJAS():
    apply_style(fontsize=22, style=None, linewidth=2) 
    
    projection = ccrs.PlateCarree()
    fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 13), subplot_kw={"projection":
                                                                                                                      projection})
    
    # add patch (Sahel: 10-20 N, 20W - 30E; coast of Guinea: 5-10N, 20W-30E, Sahara region: 20-30N, 20W-30E)
    
    #sahara --> 30, 20 N 20W, 30 E
    
    lat_sh, h_sh = 20 , 10  # lat and height
    lon_sh , w_sh = -20, 50  # long and width 
    
    #sahel --> 20, 10 N 20W, 30 E
    
    lat_sl, h_sl = 10 , 9.5  # lat and height
    lon_sl , w_sl = -20, 50  # long and width 
    
    #guinea -->5, 10 N 20W, 30 E
    
    lat_g, h_g = 5 , 4.5  # lat and height
    lon_g , w_g = -20, 50  # long and width 
    
    
    plot_annual_mean(ax=ax1, variable="Precipitation", data_alt=PI_prec_alt, cmap=YlGnBu, units="mm/month", vmax=350, vmin=50, domain="West Africa", 
                      levels=22, level_ticks=6, title="[A]", left_labels=True, bottom_labels=False, use_colorbar_default=True,
                      )
    
    ax1.add_patch(patches.Rectangle(xy =(lon_sh, lat_sh), width= w_sh, height=h_sh, ls= "--", color= red, transform = projection, 
                                    fc="None", lw=2.5,))
    ax1.add_patch(patches.Rectangle(xy =(lon_sl, lat_sl), width= w_sl, height=h_sl, ls= "--", color= black, transform = projection, 
                                    fc="None", lw=2.5,))
    ax1.add_patch(patches.Rectangle(xy =(lon_g, lat_g), width= w_g, height=h_g, ls= "--", color= blue, transform = projection, 
                                    fc="None", lw=2.5,))
    
    
    plot_annual_mean(ax=ax2, variable="Temperature", data_alt=PI_t2m_alt, cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                      levels=22, level_ticks=11, title="[B]", left_labels=True, bottom_labels=True, use_colorbar_default=True)
    
    plot_annual_mean(ax=ax3, variable="Sea Level Pressure", data_alt=PI_slp_alt, cmap=RdYlBu_r, units="hPa", vmax=1020, vmin=1000, domain="West Africa", 
                      levels=22, level_ticks=6, title="[C]", left_labels=True, bottom_labels=True, use_colorbar_default=True,
                      plot_winds=True, data_u=PI_u850_alt, data_v=PI_v850_alt)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
    plt.tight_layout() 
    plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.06)
    
    plt.savefig(os.path.join(path_to_store, "PI_tp_t2m_slp.png"), format= "png", bbox_inches="tight", dpi=300)
    
    
apply_style(fontsize=22, style=None, linewidth=2) 

projection = ccrs.PlateCarree()
fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 13), subplot_kw={"projection":
                                                                                                                  projection})
plot_annual_mean(variable="Precipitation anomalies", data_alt=MH_prec_alt_diff , cmap=BrBG, units="mm/month", ax=ax1, fig=fig, vmax=150, vmin=-150,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                           title= ["[A]  MH - PI"], bottom_labels=True, left_labels=True, compare_data1=PI_prec,
                           compare_data2=MH_prec, max_pvalue=0.01, plot_stats=True, time="JJAS")

plot_annual_mean(variable="Precipitation anomalies", data_alt=LGM_prec_alt_diff , cmap=BrBG, units="mm/month", ax=ax2, fig=fig, vmax=150, vmin=-150,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=True, cbar_pos = [0.35, 0.25, 0.25, 0.02],
                           title= ["[B]  LGM - PI"], bottom_labels=True, left_labels=False, compare_data1=PI_prec,
                           compare_data2=LGM_prec, max_pvalue=0.01, plot_stats=True, orientation= "horizontal", time="JJAS") 

plot_annual_mean(variable="Precipitation anomalies", data_alt=PLIO_prec_alt_diff , cmap=BrBG, units="mm/month", ax=ax3, fig=fig, vmax=150, vmin=-150,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                           title= ["[C]  PLIO - PI"], bottom_labels=True, left_labels=False, compare_data1=PI_prec,
                           compare_data2=PLIO_prec, max_pvalue=0.01, plot_stats=True, time="JJAS") 
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.15)
plt.savefig(os.path.join(path_to_store, "prec_anomalies.png"), format= "png", bbox_inches="tight", dpi=300)

projection = ccrs.PlateCarree()
fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 13), subplot_kw={"projection":
                                                                                                                  projection})
plot_annual_mean(variable="Temperature anomalies", data_alt=MH_t2m_alt_diff , cmap=RdBu_r, units="°C", ax=ax1, fig=fig, vmax=12, vmin=-12,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                           title= ["[A]  MH - PI"], bottom_labels=True, left_labels=True, compare_data1=PI_t2m,
                           compare_data2=MH_t2m, max_pvalue=0.01, plot_stats=True, time="JJAS")

plot_annual_mean(variable="Temperature anomalies", data_alt=LGM_t2m_alt_diff , cmap=RdBu_r, units="°C", ax=ax2, fig=fig, vmax=12, vmin=-12,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=True, cbar_pos = [0.35, 0.25, 0.25, 0.02],
                           title= ["[B]  LGM - PI"], bottom_labels=True, left_labels=False, compare_data1=PI_t2m,
                           compare_data2=LGM_t2m, max_pvalue=0.01, plot_stats=True, orientation= "horizontal", time="JJAS") 

plot_annual_mean(variable="Temperature anomalies", data_alt=PLIO_t2m_alt_diff , cmap=RdBu_r, units="°C", ax=ax3, fig=fig, vmax=12, vmin=-12,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                           title= ["[C]  PLIO - PI"], bottom_labels=True, left_labels=False, compare_data1=PI_t2m,
                           compare_data2=PLIO_t2m, max_pvalue=0.01, plot_stats=True, time="JJAS") 
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.15)
plt.savefig(os.path.join(path_to_store, "t2m_anomalies.png"), format= "png", bbox_inches="tight", dpi=300)

projection = ccrs.PlateCarree()
fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 13), subplot_kw={"projection":
                                                                                                                  projection})
plot_annual_mean(variable="Sea Level Pressure anomalies", data_alt=MH_slp_alt_diff , cmap=RdBu, units="hPa", ax=ax1, fig=fig, vmax=15, vmin=-10,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                           title= ["[A]  MH - PI"], bottom_labels=True, left_labels=True, plot_winds=True, data_u=MH_u850_alt, data_v=MH_v850_alt, time="JJAS")

plot_annual_mean(variable="Sea Level Pressure anomalies", data_alt=LGM_slp_alt_diff , cmap=RdBu, units="hPa", ax=ax2, fig=fig, vmax=15, vmin=-10,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=True, cbar_pos = [0.35, 0.25, 0.25, 0.02],
                           title= ["[B]  LGM - PI"], bottom_labels=True, left_labels=False, plot_winds=True, data_u=LGM_u850_alt, data_v=LGM_v850_alt, 
                           orientation= "horizontal", time="JJAS") 

plot_annual_mean(variable="Sea Level Pressure anomalies", data_alt=PLIO_slp_alt_diff , cmap=RdBu, units="hPa", ax=ax3, fig=fig, vmax=15, vmin=-10,
                 levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                           title= ["[C]  PLIO - PI"], bottom_labels=True, left_labels=False, plot_winds=True, data_u=PLIO_u850_alt, data_v=PLIO_v850_alt, 
                           time="JJAS") 
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.15)
plt.savefig(os.path.join(path_to_store, "slp_anomalies.png"), format= "png", bbox_inches="tight", dpi=300)
plt.show()
    
    

# left, bottom, width, height
# fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(20, 13))

# plot_monthly_sections(data_sahara=PI_month_sahara_prec, data_sahel=PI_month_sahel_prec, data_guinea=PI_month_guinea_prec,
#                       ax=ax1, ymax=250, ymin=0, varname="Precipitation [mm/month]", title="[A]")

# plot_monthly_sections(data_sahara=PI_month_sahara_t2m, data_sahel=PI_month_sahel_t2m, data_guinea=PI_month_guinea_t2m,
#                       ax=ax2, ymax=35, ymin=5, varname="Temperature [°C]", title="[B]")

# ax2.legend(bbox_to_anchor=(1.04,1),frameon=True, fontsize=25, loc="upper left")
# plt.tight_layout()
# plt.savefig(os.path.join(path_to_store, "PI_monthly_sections.png"), format= "png", bbox_inches="tight", dpi=300)
# plt.show()