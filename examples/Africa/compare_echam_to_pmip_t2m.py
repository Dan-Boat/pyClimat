#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 18:00:15 2023

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
flag = 4 # 1 MH, 2 LGM, 3 mPLIO, 4 PI

path_to_plots = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

main_path_mh = "/home/dboateng/Model_output_pst/PMIP_postprocessed/MH"

main_path_lgm = "/home/dboateng/Model_output_pst/PMIP_postprocessed/LGM"

main_path_plio = "/home/dboateng/Model_output_pst/PMIP_postprocessed/mPILO"

main_path_pi = "/home/dboateng/Model_output_pst/PMIP_postprocessed/PI"

if flag == 1:
    
    #MH
    awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
    cesm_path = os.path.join(main_path_mh, "CESM2")
    ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
    giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
    hadGEM_path = os.path.join(main_path_mh, "HadGEM3-GC31-LL")
    ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
    miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
    mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")

elif flag == 2:
    #LGM
    awi_path = os.path.join(main_path_lgm, "AWI-ESM-1-1-LR")
    cesm_waccm_path = os.path.join(main_path_lgm, "CESM2-WACCM-FV2")
    inm_cm_path = os.path.join(main_path_lgm, "INM-CM4-8")
    miroc_path = os.path.join(main_path_lgm, "MIROC-ES2L")
    mpi_esm_path = os.path.join(main_path_lgm, "MPI-ESM1-2-LR")

elif flag ==3:
    #PLIO
    cesm_path = os.path.join(main_path_plio, "CESM2")
    ec_earth_path = os.path.join(main_path_plio, "EC-Earth3-LR")
    giss_path = os.path.join(main_path_plio, "GISS-E2-1-G")
    hadGEM_path = os.path.join(main_path_plio, "HadGEM3-GC31-LL")
    ipsl_path = os.path.join(main_path_plio, "IPSL-CM6A-LR")
    norESM_path = os.path.join(main_path_plio, "NorESM1-F")

else:
    #PI
    awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
    cesm_path = os.path.join(main_path_mh, "CESM2")
    cesm_waccm_path = os.path.join(main_path_lgm, "CESM2-WACCM-FV2")
    ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
    giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
    hadGEM_path = os.path.join(main_path_mh, "HadGEM3-GC31-LL")
    inm_cm_path = os.path.join(main_path_lgm, "INM-CM4-8")
    ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
    miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
    mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")
    norESM_path = os.path.join(main_path_plio, "NorESM1-F")


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

if flag == 1:
    awi_pr_alt, awi_t2m_alt = read_t2m_tp_from_path(awi_path)
    cesm_pr_alt, cesm_t2m_alt = read_t2m_tp_from_path(cesm_path)
    ec_earth_pr_alt, ec_earth_t2m_alt = read_t2m_tp_from_path(ec_earth_path)
    giss_pr_alt, giss_t2m_alt = read_t2m_tp_from_path(giss_path)
    hadGEM_pr_alt, hadGEM_t2m_alt = read_t2m_tp_from_path(hadGEM_path)
    ipsl_pr_alt, ipsl_t2m_alt = read_t2m_tp_from_path(ipsl_path)
    miroc_pr_alt, miroc_t2m_alt = read_t2m_tp_from_path(miroc_path)
    mpi_esm_pr_alt, mpi_esm_t2m_alt = read_t2m_tp_from_path(mpi_esm_path)
    
elif flag ==2:
    awi_pr_alt, awi_t2m_alt = read_t2m_tp_from_path(awi_path)
    cesm_waccm_pr_alt, cesm_waccm_t2m_alt = read_t2m_tp_from_path(cesm_waccm_path)
    inm_pr_alt, inm_t2m_alt = read_t2m_tp_from_path(inm_cm_path)
    miroc_pr_alt, miroc_t2m_alt = read_t2m_tp_from_path(miroc_path)
    mpi_esm_pr_alt, mpi_esm_t2m_alt = read_t2m_tp_from_path(mpi_esm_path)
    
elif flag ==3:
    
    cesm_pr_alt, cesm_t2m_alt = read_t2m_tp_from_path(cesm_path)
    ec_earth_pr_alt, ec_earth_t2m_alt = read_t2m_tp_from_path(ec_earth_path)
    giss_pr_alt, giss_t2m_alt = read_t2m_tp_from_path(giss_path)
    hadGEM_pr_alt, hadGEM_t2m_alt = read_t2m_tp_from_path(hadGEM_path)
    ipsl_pr_alt, ipsl_t2m_alt = read_t2m_tp_from_path(ipsl_path)
    norESM_pr_alt, norESM_t2m_alt = read_t2m_tp_from_path(norESM_path)
    
    
elif flag ==4:
    
    awi_pr_alt, awi_t2m_alt = read_t2m_tp_from_path(awi_path)
    cesm_pr_alt, cesm_t2m_alt = read_t2m_tp_from_path(cesm_path)
    cesm_waccm_pr_alt, cesm_t2m_alt = read_t2m_tp_from_path(cesm_waccm_path)
    ec_earth_pr_alt, ec_earth_t2m_alt = read_t2m_tp_from_path(ec_earth_path)
    giss_pr_alt, giss_t2m_alt = read_t2m_tp_from_path(giss_path)
    hadGEM_pr_alt, hadGEM_t2m_alt = read_t2m_tp_from_path(hadGEM_path)
    ipsl_pr_alt, ipsl_t2m_alt = read_t2m_tp_from_path(ipsl_path)
    miroc_pr_alt, miroc_t2m_alt = read_t2m_tp_from_path(miroc_path)
    mpi_esm_pr_alt, mpi_esm_t2m_alt = read_t2m_tp_from_path(mpi_esm_path)
    norESM_pr_alt, norESM_t2m_alt = read_t2m_tp_from_path(norESM_path)
    
else:
    print("Define the flag .....")

    
    
def plot_t2m(labels, data, axes, figname):
    
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
    plt.savefig(os.path.join(path_to_plots, figname), format= "svg", bbox_inches="tight", dpi=300)
    
from read_data import MH_data, LGM_data, PLIO_data, PI_data


apply_style(fontsize=22, style=None, linewidth=2) 

projection = ccrs.PlateCarree()


if flag == 1:
    mh_echam_pr_alt, mh_echam_t2m_alt = read_t2m_tp_from_path(echam=True, data=MH_data)

    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols = 3, figsize=(28, 22), 
                                                                           subplot_kw={"projection": projection})
    axes = [ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    
    data_mh = [mh_echam_t2m_alt, awi_t2m_alt, cesm_t2m_alt, ec_earth_t2m_alt, giss_t2m_alt, hadGEM_t2m_alt,
            ipsl_t2m_alt, miroc_t2m_alt, mpi_esm_t2m_alt]
    
    labels_mh = ["ECHAM5-wiso", "AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
              "HadGEM3-GC31-LL", "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"]
    
    figname = "compare_pmip_mh_and_echam_t2m.svg"
    
    plot_t2m(labels=labels_mh, data=data_mh, axes=axes, figname=figname)
    
    
elif flag == 2:
    
    lgm_echam_pr_alt, lgm_echam_t2m_alt = read_t2m_tp_from_path(echam=True, data=LGM_data)
    
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(28, 22), 
                                                                           subplot_kw={"projection": projection})
    axes = [ax1,ax2, ax3, ax4, ax5, ax6]
    
    
    data_lgm = [lgm_echam_t2m_alt, awi_t2m_alt, cesm_waccm_t2m_alt, inm_t2m_alt, miroc_t2m_alt, mpi_esm_t2m_alt]
    
    labels_lgm = ["ECHAM5-wiso", "AWI-ESM-1-1-LR", "CESM2-WACCM-FV2", "INM-CM4-8", "MIROC-ES2L", "MPI-ESM1-2-LR"]
    
    figname = "compare_pmip_lgm_and_echam_t2m.svg"
    
    plot_t2m(labels=labels_lgm, data=data_lgm, axes=axes, figname=figname)

elif flag == 3:
    
    plio_echam_pr_alt, plio_echam_t2m_alt = read_t2m_tp_from_path(echam=True, data=PLIO_data)
    
    data_plio = [plio_echam_t2m_alt, cesm_t2m_alt, ec_earth_t2m_alt, giss_t2m_alt, hadGEM_t2m_alt,
            ipsl_t2m_alt, norESM_t2m_alt]
    
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(nrows = 3, ncols = 3, figsize=(28, 22), 
                                                                           subplot_kw={"projection": projection})
    axes = [ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    
    labels_plio = ["ECHAM5-wiso", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
              "HadGEM3-GC31-LL", "IPSL-CM6A-LR", "NorESM1-F"]
    
    
    figname = "compare_pmip_plio_and_echam_t2m.svg"
    
    plot_t2m(labels=labels_plio, data=data_plio, axes=axes, figname=figname)
    
    
    
elif flag ==4:
    
    pi_echam_pr_alt, pi_echam_t2m_alt = read_t2m_tp_from_path(echam=True, data=PI_data)
    
    data_pi = [pi_echam_t2m_alt, awi_t2m_alt, cesm_t2m_alt, cesm_t2m_alt, ec_earth_t2m_alt, giss_t2m_alt, hadGEM_t2m_alt,
            ipsl_t2m_alt, miroc_t2m_alt, mpi_esm_t2m_alt, norESM_t2m_alt]

    labels_pi = ["ECHAM5-wiso", "AWI-ESM-1-1-LR", "CESM2", "CESM2-WACCM-FV2", "EC-Earth3-LR", "GISS-E2-1-G", 
              "HadGEM3-GC31-LL", "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR", "NorESM1-F"]
    
    fig, ((ax1,ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12)) = plt.subplots(nrows = 4, ncols = 3, figsize=(28, 25), 
                                                                           subplot_kw={"projection": projection})
    

    axes = [ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]

    figname = "compare_pmip_pi_and_echam_t2m.svg"
    
    plot_t2m(labels=labels_pi, data=data_pi, axes=axes, figname=figname)
    
    
else:
    print("Define the right flag for the experiment")

plt.show()  

# extract vars and compute mean for JJA 
