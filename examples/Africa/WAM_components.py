# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:14:50 2024

@author: dboateng
"""

import os 
import pandas as pd 

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean 
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff

# define paths 
ERA5_path = "D:/Datasets/ERA5/monthly_1950_2021/"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

ERA5_tp_path = os.path.join(ERA5_path, "tp_monthly.nc")
ERA5_t2m_path = os.path.join(ERA5_path, "t2m_monthly.nc")
ERA5_v_path = os.path.join(ERA5_path, "v10_monthly.nc")
ERA5_u_path = os.path.join(ERA5_path, "u10_monthly.nc")

#read in postprocessed and analysed data 
ERA5_t2m = read_ERA_processed(path=ERA5_t2m_path, varname="t2m")   - 273.15 #°C
ERA5_tp = read_ERA_processed(path=ERA5_tp_path, varname="tp") * 1000 * 30  #mm/month
ERA5_v10 = read_ERA_processed(path=ERA5_v_path, varname="v10") #m/s
ERA5_u10 = read_ERA_processed(path=ERA5_u_path, varname="u10") #m/s


from1980to2014 = pd.date_range(start="1980-01-01", end="2014-12-31", freq="MS")

ERA5_t2m_alt = compute_lterm_mean(data=ERA5_t2m, time="month", month="JJAS", time_range=from1980to2014)
ERA5_tp_alt = compute_lterm_mean(data=ERA5_tp, time="month", month="JJAS", time_range=from1980to2014)
ERA5_v10_alt = compute_lterm_mean(data=ERA5_v10, time="month", month="JJAS", time_range=from1980to2014)
ERA5_u10_alt = compute_lterm_mean(data=ERA5_u10, time="month", month="JJAS", time_range=from1980to2014)




projection = ccrs.PlateCarree()
fig,(ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(24, 15), 
                             subplot_kw={"projection": projection})
    
    
plot_annual_mean(ax=ax1, fig=fig, variable="Precipitation", data_alt=ERA5_tp_alt, cmap=YlGnBu, units="mm/month", vmax=400, vmin=50, domain="West Africa", 
                  levels=22, level_ticks=6, title="ERA5", left_labels=True, bottom_labels=False, 
                  add_colorbar=True, orientation= "horizontal", plot_winds=True,data_u=ERA5_u10_alt,
                  data_v=ERA5_v10_alt, cbar_pos = [0.05, 0.05, 0.25, 0.02])


plot_annual_mean(ax=ax2, fig=fig, variable="Temperature", data_alt=ERA5_t2m_alt, cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                  levels=22, level_ticks=11, title="ERA5", left_labels=True, bottom_labels=True, 
                  add_colorbar=True, orientation= "horizontal", cbar_pos = [0.50, 0.05, 0.25, 0.02])


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, "WAM_temp_prec.pdf"), format= "pdf", bbox_inches="tight", dpi=300)