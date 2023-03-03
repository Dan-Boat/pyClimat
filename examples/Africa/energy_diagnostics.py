# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:34:43 2023

@author: dboateng

ploting the sensible and latent heat flux anomalies across the WAM region
Calculate the MSE and plot their meridoinal gradient
"""

import os
import xarray as xr 
import numpy as np
import pandas as pd 


# from pyClimat
from pyClimat.data import read_from_path
from pyClimat.analysis import extract_var, compute_lterm_diff, compute_lterm_mean
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean


main_path = "D:/Datasets/Model_output_pst"
path_to_store = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

lgm_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
plio_path = os.path.join(main_path, "PLIO", "MONTHLY_MEANS")
mh_path = os.path.join(main_path, "MH", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")
pd_path = os.path.join(main_path, )

filename_lterm = "1003_1017_1m_mlterm.nc"


PI_data = read_from_path(pi_path, filename_lterm)
LGM_data = read_from_path(lgm_path, filename_lterm)
PLIO_data = read_from_path(plio_path, filename_lterm)
MH_data = read_from_path(mh_path, filename_lterm)


# extract variables
def extract_energy_vars(data):
    sf = extract_var(Dataset=data, varname="sensible heat")
    lf = extract_var(Dataset=data, varname="latent heat")
    ep = extract_var(Dataset=data, varname="E-P", units="mm/month")
    
    sf_alt = compute_lterm_mean(data=sf, time="month", month="JJAS")
    lf_alt = compute_lterm_mean(data=lf, time="month", month="JJAS")
    ep_alt = compute_lterm_mean(data=ep, time="month", month="JJAS")
    
        
    
    return sf_alt, lf_alt, ep_alt

PI_sf, PI_lf, PI_ep = extract_energy_vars(data=PI_data)
MH_sf, MH_lf, MH_ep = extract_energy_vars(data=MH_data)
LGM_sf, LGM_lf, LGM_ep = extract_energy_vars(data=LGM_data)
PLIO_sf, PLIO_lf, PLIO_ep = extract_energy_vars(data=PLIO_data)


#plot 
apply_style(fontsize=22, style=None, linewidth=2) 
    
projection = ccrs.PlateCarree()
# fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 13), 
#                                    subplot_kw={"projection": projection})

# plot_annual_mean(variable="Latent heat flux anomaly", data_alt=MH_lf - PI_lf , cmap="PRGn", units="W/m²", ax=ax1, fig=fig, vmax=80, vmin=-80,
#                     levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
#                               title= ["(a) MH-PI"], bottom_labels=True, left_labels=True, time="JJAS")

# plot_annual_mean(variable="Latent Heat Flux", data_alt=LGM_lf - PI_lf , cmap="PRGn", units="W/m²", ax=ax2, fig=fig, vmax=80, vmin=-80,
#                     levels=22, domain="West Africa", level_ticks=11, add_colorbar=True, cbar_pos = [0.35, 0.25, 0.25, 0.02],
#                               title= ["(b) LGM-PI"], bottom_labels=True, left_labels=False, orientation= "horizontal", time="JJAS") 

# plot_annual_mean(variable="Latent heat flux anomaly", data_alt=PLIO_lf - PI_lf , cmap="PRGn", units="W/m²", ax=ax3, fig=fig, vmax=80, vmin=-80,
#                     levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
#                               title= ["(c) mPLIO-PI"], bottom_labels=True, left_labels=True, time="JJAS")
    
    
# fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
# plt.tight_layout()
# plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.15)
# plt.savefig(os.path.join(path_to_store, "lf_anomalies.svg"), format= "svg", bbox_inches="tight", dpi=300)

fig, ((ax1,ax2),(ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(24, 18), 
                                    subplot_kw={"projection": projection})

plot_annual_mean(variable="Latent heat flux", data_alt=PI_lf , cmap=Spectral_r, units="W/m²", ax=ax1, fig=fig, vmax=30, vmin=-150,
                    levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                              title= ["(a) PI"], bottom_labels=True, left_labels=True, time="JJAS", center=False)

plot_annual_mean(variable="Latent heat flux ", data_alt=MH_lf, cmap=Spectral_r, units="W/m²", ax=ax2, fig=fig, vmax=30, vmin=-150,
                    levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                              title= ["(b) MH"], bottom_labels=True, left_labels=True, time="JJAS", center=False)

plot_annual_mean(variable="Latent Heat Flux", data_alt=LGM_lf, cmap=Spectral_r, units="W/m²", ax=ax3, fig=fig, vmax=30, vmin=-150,
                    levels=22, domain="West Africa", level_ticks=11, add_colorbar=True, cbar_pos = [0.35, 0.01, 0.25, 0.02],
                              title= ["(c) LGM"], bottom_labels=True, left_labels=False, orientation= "horizontal", time="JJAS", 
                              center=False) 

plot_annual_mean(variable="Latent heat flux", data_alt=PLIO_lf, cmap=Spectral_r, units="W/m²", ax=ax4, fig=fig, vmax=10, vmin=-150,
                    levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                              title= ["(d) mPLIO"], bottom_labels=True, left_labels=True, time="JJAS", 
                              center=False)
    
    
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.15)
plt.savefig(os.path.join(path_to_store, "lf_magnitudes.png"), format= "png", bbox_inches="tight", dpi=300)

fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(28, 13), 
                                   subplot_kw={"projection": projection})

plot_annual_mean(variable="E-P", data_alt=MH_ep , cmap="bwr", units="mm/month", ax=ax1, fig=fig, vmax=150, vmin=-150,
                    levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                              title= ["(a) MH"], bottom_labels=True, left_labels=True, time="JJAS")

plot_annual_mean(variable="E-P", data_alt=LGM_ep, cmap="bwr", units="mm/month", ax=ax2, fig=fig, vmax=150, vmin=-150,
                    levels=22, domain="West Africa", level_ticks=11, add_colorbar=True, cbar_pos = [0.35, 0.25, 0.25, 0.02],
                              title= ["(b) LGM"], bottom_labels=True, left_labels=False, orientation= "horizontal", time="JJAS") 

plot_annual_mean(variable="E-P", data_alt=PLIO_ep, cmap="bwr", units="mm/month", ax=ax3, fig=fig, vmax=150, vmin=-150,
                    levels=22, domain="West Africa", level_ticks=11, add_colorbar=False,
                              title= ["(c) mPLIO"], bottom_labels=True, left_labels=True, time="JJAS")
    
    
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.94, bottom=0.15)
plt.savefig(os.path.join(path_to_store, "ep.svg"), format= "svg", bbox_inches="tight", dpi=300)

