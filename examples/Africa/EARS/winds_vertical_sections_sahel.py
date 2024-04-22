# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 13:49:19 2024

@author: dboateng
"""

import os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl 
import cartopy.crs as ccrs
from cycler import cycler

import geocat.comp as gc


# import pyClimat models 
from pyClimat.plot_utils import *
from pyClimat.plots import plot_annual_mean
from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, compute_lterm_diff, extract_transect
from pyClimat.variables import extract_var

projection_trans = ccrs.PlateCarree()


path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"


# define paths
pi_path = "D:/Datasets/Model_output_pst/PI/MONTHLY_MEANS/"
miocene_path= "D:/Datasets/Model_output_pst/a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"
pliocene_path= "D:/Datasets/Model_output_pst/PLIO/MONTHLY_MEANS/"
holocene_path = "D:/Datasets/Model_output_pst/MH/MONTHLY_MEANS/"
miocene_path = "D:/Datasets/Model_output_pst/a014_hpc-bw_e5w2.3_t159_MIO_W1E1_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"
miocene_ea_high_path = "D:/Datasets/Model_output_pst/a025_dkrz-levante_e5w2.3_t159_MIO_EA_high_450ppm_t159l31.6h/output_processed/MONTHLY_MEANS/"

echam_dynamics_filename = "1003_1017_1m_mlterm_dynamics.nc"

mmco_slm_path = "D:/Datasets/topo/Miotopofiles/CTL_Mio_Herold/topo_Herold_Miocene_2160x1080_SLM.nc"

mio_slm_data = xr.open_dataset(mmco_slm_path)
mio_slm = mio_slm_data.SLM



def read_winds_q(path):

    data = read_from_path(path=path, filename=echam_dynamics_filename)
    u_600 = extract_var(Dataset=data, varname="u", lev=500, lev_units="hPa")
    u_925 = extract_var(Dataset=data, varname="u", lev=925, lev_units="hPa")
    q = extract_var(Dataset=data, varname="q", lev=600, lev_units="hPa", units="g/kg")
    
    
    u_alt=compute_lterm_mean(data=u_600, 
                                 time="month", month="JJAS")
    
    u_925_alt=compute_lterm_mean(data=u_925, 
                                 time="month", month="JJAS")
    
    q_alt=compute_lterm_mean(data=q, 
                                 time="month", month="JJAS")
    
    region_u = extract_transect(data=u_925_alt, maxlon=40, maxlat=30, minlon=-20, minlat=10)
    region_u = region_u.sortby(region_u.lon)
    
    return q_alt, region_u, u_alt


mh_q, mh_u_925, mh_u_600 = read_winds_q(path=holocene_path)
plio_q, plio_u_925, plio_u_600 = read_winds_q(path=pliocene_path)
mio_q, mio_u_925, mio_u_600 = read_winds_q(path=miocene_path)

# region_u = extract_transect(data=u_alt, maxlon=20, maxlat=15, minlon=-40, minlat=-25)
# region_u = region_v.sortby(region_u.lon)

apply_style2(fontsize=24, style=None, linewidth=2.5, usetex=True) 
            
projection = ccrs.Robinson(central_longitude=0, globe=None)

projection_trans = ccrs.PlateCarree()


fig,(ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(26,13), subplot_kw={"projection":projection})

plot_annual_mean(variable="q", data_alt=mio_q, ax=ax1,
                  cmap="Spectral", units="g/kg [JJAS]", vmax=8, vmin=0, 
                levels=22, level_ticks=11, add_colorbar=True, plot_coastlines=False, bottom_labels=True,
                left_labels=True, fig=fig, plot_borders=False, plot_projection=projection, title="(a) Mid-Miocene - PI", 
                orientation="horizontal",  cbar_pos= [0.30, 0.01, 0.45, 0.03], sea_land_mask=mio_slm,
                domain="West Africa", center=False)

c1 = mio_u_600.plot.contour(colors="k", linestyles="-", ax=ax1, transform=projection_trans, linewidths=2.0,
                                       add_labels=False, levels=[-10.0])
#c1.clabel(fmt="%.1f", use_clabeltext=True, colors="black", fontsize=22)

c2 = mio_u_925.plot.contour(linestyles="--", ax=ax1, transform=projection_trans, linewidths=2.0,
                                       add_labels=False, levels=[0.00])
#c2.clabel(fmt="%.1f", use_clabeltext=True, fontsize=22)


plot_annual_mean(variable="q", data_alt=plio_q, ax=ax2,
                  cmap="Spectral", units="g/kg [JJAS]", vmax=8, vmin=0, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=True,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(b) Mid-Pliocene - PI", 
                domain="West Africa", center=False)

c1 = plio_u_600.plot.contour(colors="k", linestyles="-", ax=ax2, transform=projection_trans, linewidths=2.0,
                                       add_labels=False, levels=[-10.0])
#c1.clabel(fmt="%.1f", use_clabeltext=True, colors="black", fontsize=22)

c2 = plio_u_925.plot.contour(linestyles="--", ax=ax2, transform=projection_trans, linewidths=2.0,
                                       add_labels=False, levels=[0.00])
#c2.clabel(fmt="%.1f", use_clabeltext=True, fontsize=22)

plot_annual_mean(variable="q", data_alt=mh_q, ax=ax3,
                  cmap="Spectral", units="g/kg [JJAS]", vmax=8, vmin=0, 
                levels=22, level_ticks=11, add_colorbar=False, plot_coastlines=True, bottom_labels=True,
                left_labels=False, fig=fig, plot_borders=False, plot_projection=projection, title="(c) Mid-Holocene - PI", 
                orientation="horizontal", domain="West Africa", center=False)

c1 = mh_u_600.plot.contour(colors="k", linestyles="-", ax=ax3, transform=projection_trans, linewidths=2.0,
                                       add_labels=False, levels=[-10.0])
#c1.clabel(fmt="%.1f", use_clabeltext=True, colors="black", fontsize=22)

c2 = mh_u_925.plot.contour(linestyles="--", ax=ax3, transform=projection_trans, linewidths=2.0,
                                       add_labels=False, levels=[0.00])
#c2.clabel(fmt="%.1f", use_clabeltext=True, fontsize=22)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.10,)
plt.savefig(os.path.join(path_to_plots, "poster_map_winds.pdf"), format= "pdf", bbox_inches="tight", dpi=600)
