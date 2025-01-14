# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:55:16 2024

@author: dboateng
"""

import os 
import pandas as pd 
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates


from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean 
from pyClimat.data import read_ERA_processed, read_ECHAM_processed, read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff
from pyClimat.variables import extract_var

main_path = "D:/Datasets/Model_output_pst/PD"
gnip_path = "D:/Datasets/GNIP_data/world/scratch/station_world_overview_5years.csv" 

path_to_store = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Alps/Miocene/plots/armelle"


df = pd.read_csv(gnip_path)


#load datasets 
PD_data = read_from_path(main_path, "PD_1980_2014_monthly.nc", decode=True)
PD_wiso = read_from_path(main_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
d18Op = extract_var(Dataset=PD_data, varname="d18op", units="per mil", Dataset_wiso=PD_wiso,
                    )

d18Op_alt = compute_lterm_mean(data=d18Op, time="annual")

apply_style(fontsize=25, style=None, linewidth=2) 

projection = ccrs.Robinson(central_longitude=0, globe=None)

lat_sh, h_sh = 46.66 , 0.8  # lat and height
lon_sh , w_sh = 6.66, 3.14  # long and width 

lat_vh, h_vh = 42.7 , 2  # lat and height
lon_vh , w_vh = 4.9, 2.5  # long and width 
    
    
#projection_p = ccrs.PlateCarree()
fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(13, 13), subplot_kw={"projection":  
                                                                        projection})




plot_annual_mean(ax=ax1, variable='$\delta^{18}$Op VSMOW', data_alt=d18Op_alt, cmap="Spectral_r", 
                  units="‰", vmax=0, vmin=-15, domain="Alps", 
                  levels=22, level_ticks=9, GNIP_data=df, title=None, left_labels=True, bottom_labels=True, 
                  use_colorbar_default=True, center=False,
                  plot_projection=projection, plot_borders=True)

ax1.add_patch(patches.Rectangle(xy =(lon_sh, lat_sh), width= w_sh, height=h_sh, ls= "--", color= red, transform = ccrs.PlateCarree(), 
                                    fc="None", lw=2.5,))

ax1.add_patch(patches.Rectangle(xy =(lon_vh, lat_vh), width= w_vh, height=h_vh, ls= "--", color= blue, transform = ccrs.PlateCarree(), 
                                    fc="None", lw=2.5,))

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "global_echam_gnip.pdf"), format= "pdf", bbox_inches="tight")


