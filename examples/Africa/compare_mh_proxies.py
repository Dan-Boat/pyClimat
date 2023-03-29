# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 14:08:06 2023

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
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

main_path_mh = "D:/Datasets/CMIP6/PMIP/postprocessed/MH"

main_path_pi = "D:/Datasets/CMIP6/PMIP/postprocessed/PI"

path_to_proxies = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/map_delta_06ka_ALL_grid_2x2.csv"


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

MH_proxies_data = pd.read_csv(path_to_proxies)


from read_data import PI_data, MH_data

def read_tp_from_path(path_pi=None, path_main=None, echam=False, data_pi=None, data_main=None,
                      stats="mean"):
    
    if echam == True:
        
        pr_pi = extract_var(Dataset=data_pi, varname="prec", units="mm/a")
        pr_main = extract_var(Dataset=data_main, varname="prec", units="mm/a")
        
    else:
        
        #convert to annual means (mm/a)
        pr_pi = read_from_path(path_pi, "pr_1m_lterm.nc", varname="pr", decode=True) *60*60*24*365  #mm/month
        
        pr_main = read_from_path(path_main, "pr_1m_lterm.nc", varname="pr", decode=True) *60*60*24*365  #mm/month
        
        
    
    pr_diff_alt = compute_lterm_diff(data_control=pr_pi, data_main=pr_main, time="annual",)
    
    
    # latitudinal cross-section
    minlat = -10
    maxlat = 35
    minlon = -20
    maxlon = 30
    lat_cross_section = extract_profile(data=pr_diff_alt, maxlon=maxlon, minlon=minlon,
                                        maxlat=maxlat, minlat=minlat, dim="lat",
                                        method=stats)
    
    return lat_cross_section

# MH
mh_awi_lat_cross_section = read_tp_from_path(path_pi=pi_awi_path, path_main=mh_awi_path)
mh_cesm_lat_cross_section = read_tp_from_path(path_pi=pi_cesm_path, path_main=mh_cesm_path)
mh_ec_earth_lat_cross_section = read_tp_from_path(path_pi=pi_ec_earth_path, path_main=mh_ec_earth_path)
mh_giss_lat_cross_section = read_tp_from_path(path_pi=pi_giss_path, path_main=mh_giss_path)
mh_hadGEM_lat_cross_section = read_tp_from_path(path_pi=pi_hadGEM_path, path_main=mh_hadGEM_path)
mh_ipsl_lat_cross_section = read_tp_from_path(path_pi=pi_ipsl_path, path_main=mh_ipsl_path)
mh_miroc_lat_cross_section = read_tp_from_path(path_pi=pi_miroc_path, path_main=mh_miroc_path)
mh_mpi_esm_lat_cross_section = read_tp_from_path(path_pi=pi_mpi_esm_path, path_main=mh_mpi_esm_path)
mh_echam_lat_cross_section = read_tp_from_path(echam=True, data_pi=PI_data, data_main=MH_data)
mh_echam_lat_cross_section_std = read_tp_from_path(echam=True, data_pi=PI_data, data_main=MH_data, stats="std")



#plot 
def plot_lat_sections_pmip(show_range=False, data_range=mh_echam_lat_cross_section_std):

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
    
    mh_model_names = ["ECHAM5-wiso", "AWI-ESM-1-1-LR", "CESM2", "EC-Earth3-LR", "GISS-E2-1-G", 
              "HadGEM3-GC31-LL", "IPSL-CM6A-LR", "MIROC-ES2L", "MPI-ESM1-2-LR"]
    
    
    apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 16,))
    
    for i,data in enumerate(mh_cross_sections):
        
        if show_range==True:
            if i == 0:
                ax.fill_betweenx(data.index, x1=data - data_range, x2=data + data_range, color=black, alpha=0.1)
                
        ax.plot(data, data.index, linewidth=3, color=mh_colors[i], label= mh_model_names[i])
        
        
    
    ax.axvline(x=0, linestyle="--", color=grey, linewidth=2)
    ax.set_xlabel('Mean Annual Precipitation (MAP) anomaly (mm/year)', fontsize=24, fontweight="bold")
    
   
    ax.set_ylabel("Latitude (°N)", fontsize=24, fontweight="bold")
        
    ax.xaxis.tick_top()   
    ax.xaxis.set_label_position('top') 
        
    #ax.set_title("(a) MH - PI", fontsize=20, fontweight="bold", loc="left")
    
    ax.scatter(x=MH_proxies_data["map_anm_mean"], y=MH_proxies_data["lat"], s=120, alpha=0.5)
    ax.errorbar(x=MH_proxies_data["map_anm_mean"], y=MH_proxies_data["lat"], xerr=MH_proxies_data["map_se_mean"],
                 yerr=None, ls="none", color="black", capthick=2, capsize=3,
                 label="(Bartlein et al., 2011)")
    ax.set_ylim(-10, 30)
    ax.legend(bbox_to_anchor=(0.01, -0.15, 1., 0.102), loc=3, ncol=3, borderaxespad=0., frameon = True, 
                  fontsize=20)
   
    plt.savefig(os.path.join(path_to_plots, "mh_compare_proxies.svg"), bbox_inches="tight", format= "svg")


plot_lat_sections_pmip(show_range=True)