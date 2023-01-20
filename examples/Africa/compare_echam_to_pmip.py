#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 08:41:06 2023

@author: dboateng
"""
import os 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, extract_var

from pyClimat.plots import plot_annual_mean
from pyClimat.plot_utils import *


# define path 
path_to_plots = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

main_path_mh = "/home/dboateng/Model_output_pst/PMIP_postprocessed/MH"

awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
cesm_path = os.path.join(main_path_mh, "CESM2")
ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
hadGEM_path = os.path.join(main_path_mh, "HadGEM3-GC31-LL")
ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")

def read_t2m_tp_from_path(path=None, echam=False, data=None):
    
    if echam == True:
        
        t2m = extract_var(Dataset=data, varname="temp2", units="°C")
        pr = extract_var(Dataset=data, varname="prec", units="mm/month")
        
    else:

        pr = read_from_path(path, "pr_1m_lterm.nc", varname="pr", decode=True) *60*60*24*30  #mm/month
        t2m = read_from_path(path, "tas_1m_lterm.nc", varname="tas", decode=True) - 273.15 #°C
    
    pr_alt = compute_lterm_mean(data=pr, time="month", month="JJAS")
    t2m_alt = compute_lterm_mean(data=t2m, time="month", month="JJAS")
    
    return pr_alt, t2m_alt


awi_pr_alt, awi_t2m_alt = read_t2m_tp_from_path(awi_path)
cesm_pr_alt, cesm_t2m_alt = read_t2m_tp_from_path(cesm_path)
ec_earth_pr_alt, ec_earth_t2m_alt = read_t2m_tp_from_path(ec_earth_path)
giss_pr_alt, giss_t2m_alt = read_t2m_tp_from_path(giss_path)
hadGEM_pr_alt, hadGEM_t2m_alt = read_t2m_tp_from_path(hadGEM_path)
ipsl_pr_alt, ipsl_t2m_alt = read_t2m_tp_from_path(ipsl_path)
miroc_pr_alt, miroc_t2m_alt = read_t2m_tp_from_path(miroc_path)
mpi_esm_pr_alt, mpi_esm_t2m_alt = read_t2m_tp_from_path(mpi_esm_path)

# echam

from read_data import MH_data
echam_pr_alt, echam_t2m_alt = read_t2m_tp_from_path(echam=True, data=MH_data)
    
apply_style(fontsize=22, style=None, linewidth=2) 

projection = ccrs.PlateCarree()
fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols = 3, figsize=(28, 22), 
                                                                       subplot_kw={"projection": projection})

axes = [ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]

data = [echam_pr_alt, awi_pr_alt, cesm_pr_alt, ec_earth_pr_alt, giss_pr_alt, hadGEM_pr_alt,
        ipsl_pr_alt, miroc_pr_alt, mpi_esm_pr_alt]

labels = ["ECHAM5-wiso", "AWI.AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
          "HadGEM3-GC31-LL", "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"]


for i,label in enumerate(labels):
    if label == "ECHAM5-wiso":
        
        plot_annual_mean(ax=axes[i], fig=fig, variable="Precipitation", data_alt=data[i], cmap=YlGnBu, units="mm/month", vmax=400, vmin=50, domain="West Africa", 
                          levels=22, level_ticks=6, title=label, left_labels=True, bottom_labels=True, 
                          add_colorbar=True, cbar_pos = [0.40, 0.05, 0.25, 0.02], orientation= "horizontal")
    else:
        
        plot_annual_mean(ax=axes[i], fig=fig, variable="Precipitation", data_alt=data[i], cmap=YlGnBu, units="mm/month", vmax=400, vmin=50, domain="West Africa", 
                          levels=22, level_ticks=6, title=label, left_labels=True, bottom_labels=True, 
                          add_colorbar=False)
        
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, "compare_pmip_mh_and_echam_prec.svg"), format= "svg", bbox_inches="tight", dpi=300)


apply_style(fontsize=22, style=None, linewidth=2) 
projection = ccrs.PlateCarree()
fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols = 3, figsize=(28, 22), 
                                                                       subplot_kw={"projection": projection})

axes = [ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]

data = [echam_t2m_alt, awi_t2m_alt, cesm_t2m_alt, ec_earth_t2m_alt, giss_t2m_alt, hadGEM_t2m_alt,
        ipsl_t2m_alt, miroc_t2m_alt, mpi_esm_t2m_alt]

labels = ["ECHAM5-wiso", "AWI.AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
          "HadGEM3-GC31-LL", "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"]


for i,label in enumerate(labels):
    if label == "ECHAM5-wiso":
        
        plot_annual_mean(ax=axes[i], fig=fig, variable="Temperature", data_alt=data[i], cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                          levels=22, level_ticks=11, title=label, left_labels=True, bottom_labels=True, 
                          add_colorbar=True, cbar_pos = [0.40, 0.05, 0.25, 0.02], orientation= "horizontal")
    else:
        
        plot_annual_mean(ax=axes[i], fig=fig, variable="Temperature", data_alt=data[i], cmap=Spectral_r, units="°C", vmax=40, vmin=10, domain="West Africa", 
                          levels=22, level_ticks=11, title=label, left_labels=True, bottom_labels=True, 
                          add_colorbar=False)
        
fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.06, hspace=0.01)
plt.savefig(os.path.join(path_to_plots, "compare_pmip_mh_and_echam_t2m.svg"), format= "svg", bbox_inches="tight", dpi=300)


plt.show()  

# extract vars and compute mean for JJA 
