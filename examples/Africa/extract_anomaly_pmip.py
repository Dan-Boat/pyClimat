#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:52:10 2023

@author: dboateng
"""

import os 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, extract_var, extract_profile, compute_lterm_diff, extract_transect, extract_profile

from pyClimat.plots import plot_annual_mean
from pyClimat.plot_utils import *


# define path 
path_to_plots = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

main_path_mh = "/home/dboateng/Model_output_pst/PMIP_postprocessed/MH"

main_path_pi = "/home/dboateng/Model_output_pst/PMIP_postprocessed/PI"

main_path_lgm = "/home/dboateng/Model_output_pst/PMIP_postprocessed/LGM"

main_path_plio = "/home/dboateng/Model_output_pst/PMIP_postprocessed/mPILO"

#PI

pi_awi_path = os.path.join(main_path_pi, "AWI-ESM-1-1-LR")
pi_cesm_path = os.path.join(main_path_pi, "CESM2")
pi_cesm_waccm_path = os.path.join(main_path_pi, "CESM2-WACCM-FV2")
pi_ec_earth_path = os.path.join(main_path_pi, "EC-Earth3-LR")
pi_giss_path = os.path.join(main_path_pi, "GISS-E2-1-G")
pi_hadGEM_path = os.path.join(main_path_pi, "HadGEM3-GC31-LL")
pi_inm_cm_path = os.path.join(main_path_pi, "INM-CM4-8")
pi_ipsl_path = os.path.join(main_path_pi, "IPSL-CM6A-LR")
pi_miroc_path = os.path.join(main_path_pi, "MIROC-ES2L")
pi_mpi_esm_path = os.path.join(main_path_pi, "MPI-ESM1-2-LR")
pi_norESM_path = os.path.join(main_path_pi, "NorESM1-F")


#MH
mh_awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
mh_cesm_path = os.path.join(main_path_mh, "CESM2")
mh_ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
mh_giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
mh_hadGEM_path = os.path.join(main_path_mh, "HadGEM3-GC31-LL")
mh_ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
mh_miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mh_mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")

#LGM 
lgm_awi_path = os.path.join(main_path_lgm, "AWI-ESM-1-1-LR")
lgm_cesm_waccm_path = os.path.join(main_path_lgm, "CESM2-WACCM-FV2")
lgm_inm_cm_path = os.path.join(main_path_lgm, "INM-CM4-8")
lgm_miroc_path = os.path.join(main_path_lgm, "MIROC-ES2L")
lgm_mpi_esm_path = os.path.join(main_path_lgm, "MPI-ESM1-2-LR")


#PLIO
plio_cesm_path = os.path.join(main_path_plio, "CESM2")
plio_ec_earth_path = os.path.join(main_path_plio, "EC-Earth3-LR")
plio_giss_path = os.path.join(main_path_plio, "GISS-E2-1-G")
plio_hadGEM_path = os.path.join(main_path_plio, "HadGEM3-GC31-LL")
plio_ipsl_path = os.path.join(main_path_plio, "IPSL-CM6A-LR")
plio_norESM_path = os.path.join(main_path_plio, "NorESM1-F")

from read_data import PI_data, MH_data, LGM_data, PLIO_data

def read_tp_from_path(path_pi=None, path_main=None, echam=False, data_pi=None, data_main=None):
    
    if echam == True:
        
        pr_pi = extract_var(Dataset=data_pi, varname="prec", units="mm/month")
        pr_main = extract_var(Dataset=data_main, varname="prec", units="mm/month")
        
    else:
        pr_pi = read_from_path(path_pi, "pr_1m_lterm.nc", varname="pr", decode=True) *60*60*24*30  #mm/month
        
        pr_main = read_from_path(path_main, "pr_1m_lterm.nc", varname="pr", decode=True) *60*60*24*30  #mm/month
        
        
    
    pr_diff_alt = compute_lterm_diff(data_control=pr_pi, data_main=pr_main, time="month", month="JJAS")
    
    # sahel anomaly regional means
    minlat = 10
    maxlat = 20
    minlon = -20
    maxlon = 30
    
    
    sahel_anomaly = extract_transect(data=pr_diff_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat)
    sahel_anomaly_mean = sahel_anomaly.data.mean()
    
    # latitudinal cross-section
    minlat = -5
    maxlat = 30
    minlon = -20
    maxlon = 30
    lat_cross_section = extract_profile(data=pr_diff_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
    
    return sahel_anomaly_mean, lat_cross_section

# MH
mh_awi_sahel_anomaly, mh_awi_lat_cross_section = read_tp_from_path(path_pi=pi_awi_path, path_main=mh_awi_path)
mh_cesm_sahel_anomaly, mh_cesm_lat_cross_section = read_tp_from_path(path_pi=pi_cesm_path, path_main=mh_cesm_path)
mh_ec_earth_sahel_anomaly, mh_ec_earth_lat_cross_section = read_tp_from_path(path_pi=pi_ec_earth_path, path_main=mh_ec_earth_path)
mh_giss_sahel_anomaly, mh_giss_lat_cross_section = read_tp_from_path(path_pi=pi_giss_path, path_main=mh_giss_path)
mh_hadGEM_sahel_anomaly, mh_hadGEM_lat_cross_section = read_tp_from_path(path_pi=pi_hadGEM_path, path_main=mh_hadGEM_path)
mh_ipsl_sahel_anomaly, mh_ipsl_lat_cross_section = read_tp_from_path(path_pi=pi_ipsl_path, path_main=mh_ipsl_path)
mh_miroc_sahel_anomaly, mh_miroc_lat_cross_section = read_tp_from_path(path_pi=pi_miroc_path, path_main=mh_miroc_path)
mh_mpi_esm_sahel_anomaly, mh_mpi_esm_lat_cross_section = read_tp_from_path(path_pi=pi_mpi_esm_path, path_main=mh_mpi_esm_path)
mh_echam_sahel_anomaly, mh_echam_lat_cross_section = read_tp_from_path(echam=True, data_pi=PI_data, data_main=MH_data)


#LGM 

lgm_awi_sahel_anomaly, lgm_awi_lat_cross_section = read_tp_from_path(path_pi=pi_awi_path, path_main=lgm_awi_path)
lgm_cesm_waccm_sahel_anomaly, lgm_cesm_waccm_lat_cross_section = read_tp_from_path(path_pi=pi_cesm_waccm_path, path_main=lgm_cesm_waccm_path)
lgm_inm_cm_sahel_anomaly, lgm_inm_cm_lat_cross_section = read_tp_from_path(path_pi=pi_inm_cm_path, path_main=lgm_inm_cm_path)
lgm_miroc_sahel_anomaly, lgm_miroc_lat_cross_section = read_tp_from_path(path_pi=pi_miroc_path, path_main=lgm_miroc_path)
lgm_mpi_esm_sahel_anomaly, lgm_mpi_esm_lat_cross_section = read_tp_from_path(path_pi=pi_mpi_esm_path, path_main=lgm_mpi_esm_path)
lgm_echam_sahel_anomaly, lgm_echam_lat_cross_section = read_tp_from_path(echam=True, data_pi=PI_data, data_main=LGM_data)

#PLIO
plio_cesm_sahel_anomaly, plio_cesm_lat_cross_section = read_tp_from_path(path_pi=pi_cesm_path, path_main=plio_cesm_path)
plio_ec_earth_sahel_anomaly, plio_ec_earth_lat_cross_section = read_tp_from_path(path_pi=pi_ec_earth_path, path_main=plio_ec_earth_path)
plio_giss_sahel_anomaly, plio_giss_lat_cross_section = read_tp_from_path(path_pi=pi_giss_path, path_main=plio_giss_path)
plio_hadGEM_sahel_anomaly, plio_hadGEM_lat_cross_section = read_tp_from_path(path_pi=pi_hadGEM_path, path_main=plio_hadGEM_path)
plio_ipsl_sahel_anomaly, plio_ipsl_lat_cross_section = read_tp_from_path(path_pi=pi_ipsl_path, path_main=plio_ipsl_path)
plio_norESM_sahel_anomaly, plio_norESM_lat_cross_section = read_tp_from_path(path_pi=pi_norESM_path, path_main=plio_norESM_path)
plio_echam_sahel_anomaly, plio_echam_lat_cross_section = read_tp_from_path(echam=True, data_pi=PI_data, data_main=PLIO_data)


mh_model_names = ["ECHAM5-wiso", "AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
          "HadGEM3-GC31-LL", "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"]

lgm_model_names = ["ECHAM5-wiso", "AWI-ESM-1-1-LR", "CESM2-WACCM-FV2", "INM-CM4-8", "MIROC-ES2L", "MPI-ESM1-2-LR"]

plio_model_names = ["ECHAM5-wiso", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
          "HadGEM3-GC31-LL", "IPSL-CM6A-LR", "NorESM1-F"]

def plot_lat_sections_pmip():

    #PMIP models colors
    echam_color = black
    awi_color = blue 
    cesm_color = green
    cesm_waccm_color = red
    ec_earth_color = gold
    giss_color = purple
    hadGEM_color = "red"
    ipsl_color =grey
    miroc_color = pink
    inm_color = "#5025BE"
    norESM_color = olive
    mpi_esm_color = cyan
    
    
    
    mh_cross_sections = [mh_echam_lat_cross_section, mh_awi_lat_cross_section, mh_cesm_lat_cross_section, mh_ec_earth_lat_cross_section, mh_giss_lat_cross_section,
                         mh_hadGEM_lat_cross_section, mh_ipsl_lat_cross_section, mh_miroc_lat_cross_section, mh_mpi_esm_lat_cross_section,
                         ]
    
    mh_colors = [echam_color, awi_color, cesm_color, ec_earth_color, giss_color, hadGEM_color, ipsl_color, miroc_color, mpi_esm_color]
    
    lgm_cross_sections = [lgm_echam_lat_cross_section, lgm_awi_lat_cross_section, lgm_cesm_waccm_lat_cross_section, lgm_inm_cm_lat_cross_section, lgm_miroc_lat_cross_section,
                         lgm_mpi_esm_lat_cross_section, ]
    
    lgm_colors = [echam_color, awi_color, cesm_waccm_color, inm_color, miroc_color, mpi_esm_color]
    
    plio_cross_sections = [plio_echam_lat_cross_section, plio_cesm_lat_cross_section, plio_ec_earth_lat_cross_section, plio_giss_lat_cross_section, plio_hadGEM_lat_cross_section,
                          plio_ipsl_lat_cross_section, plio_norESM_lat_cross_section, ]
    
    plio_colors = [echam_color, cesm_color, ec_earth_color, giss_color, hadGEM_color, ipsl_color, norESM_color]
    
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(25, 16,), sharey=True)
    
    for i,data in enumerate(mh_cross_sections):
        ax1.plot(data, data.index, linewidth=3, color=mh_colors[i], label= mh_model_names[i])
        
        
        
    for i,data in enumerate(lgm_cross_sections):
        ax2.plot(data, data.index, linewidth=3, color=lgm_colors[i], label= lgm_model_names[i])
        
        
    for i,data in enumerate(plio_cross_sections):
        ax3.plot(data, data.index, linewidth=3, color=plio_colors[i], label= plio_model_names[i])
        
        
        
    axes = [ax1, ax2, ax3]
    
    for ax in axes: 
        ax.legend(bbox_to_anchor=(0.01, -0.35, 1., 0.102), loc=3, ncol=1, borderaxespad=0., frameon = True, 
                      fontsize=20)
        ax.axvline(x=0, linestyle="--", color=grey, linewidth=2)
        ax.set_xlabel('Precipitation anomaly (mm/month)', fontsize=24, fontweight="bold")
        
        if ax == ax1:
            ax.set_ylabel("Latitude (°N)", fontsize=24, fontweight="bold")
            
        ax.xaxis.tick_top()   
        ax.xaxis.set_label_position('top') 
        
    ax1.set_title("(a) MH - PI", fontsize=20, fontweight="bold", loc="left")
    ax2.set_title("(b) LGM - PI", fontsize=20, fontweight="bold", loc="left")
    ax3.set_title("(c) mPLIO - PI", fontsize=20, fontweight="bold", loc="left")
    plt.savefig(os.path.join(path_to_plots, "pmip_lat_sections_anomalies.svg"), bbox_inches="tight", format= "svg")

# try plots (before converting into function)
def plot_the_anomalies_from_pmip():
    #MH
    mh_means = [mh_echam_sahel_anomaly, mh_awi_sahel_anomaly, mh_cesm_sahel_anomaly, mh_ec_earth_sahel_anomaly, mh_giss_sahel_anomaly,
             mh_hadGEM_sahel_anomaly, mh_ipsl_sahel_anomaly, mh_miroc_sahel_anomaly, mh_mpi_esm_sahel_anomaly]
    
    
    #LGM
    
    lgm_means = [lgm_echam_sahel_anomaly, lgm_awi_sahel_anomaly, lgm_cesm_waccm_sahel_anomaly, lgm_inm_cm_sahel_anomaly,
                 lgm_miroc_sahel_anomaly, lgm_mpi_esm_sahel_anomaly]
    
    #PLIO
    
    plio_means = [plio_echam_sahel_anomaly, plio_cesm_sahel_anomaly, plio_ec_earth_sahel_anomaly, plio_giss_sahel_anomaly, plio_hadGEM_sahel_anomaly,
                  plio_ipsl_sahel_anomaly, plio_norESM_sahel_anomaly]
    
    
    
    df_mh_anomaly = pd.DataFrame(index=mh_model_names, columns=["Anomalies"])
    df_lgm_anomaly = pd.DataFrame(index=lgm_model_names, columns=["Anomalies"])
    df_plio_anomaly = pd.DataFrame(index=plio_model_names, columns=["Anomalies"])
    
    
    
    
    for i,model in enumerate(mh_model_names):
        df_mh_anomaly.loc[model] = mh_means[i]
        
    
    for i,model in enumerate(lgm_model_names):
        df_lgm_anomaly.loc[model] = lgm_means[i]
        
    
    for i,model in enumerate(plio_model_names):
        df_plio_anomaly.loc[model] = plio_means[i]
        
    
    apply_style(fontsize=22, style=None, linewidth=2, ) 
    
    
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(26, 13), sharey=False)
    df_mh_anomaly.plot(kind="bar", color="black", fontsize=22, ax=ax1, legend=False, width=0.8)
    df_lgm_anomaly.plot(kind="bar", color="black", fontsize=22, ax=ax2, legend=False, width=0.8)
    df_plio_anomaly.plot(kind="bar", color="black", fontsize=22, ax=ax3, legend=False, width=0.8)
    
    ax1.set_ylabel("Precipitation anomaly (mm/month)", fontsize=24, fontweight="bold")
    ax1.tick_params(which="minor")
    ax1.set_title("(a) MH - PI", fontsize=20, fontweight="bold", loc="left")
    ax2.set_title("(b) LGM - PI", fontsize=20, fontweight="bold", loc="left")
    ax3.set_title("(c) mPLIO - PI", fontsize=20, fontweight="bold", loc="left")
    plt.savefig(os.path.join(path_to_plots, "pmip_sahel_anomalies.svg"), bbox_inches="tight", format= "svg")
    
    
plot_lat_sections_pmip()
plot_the_anomalies_from_pmip()