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

#PI

pi_awi_path = os.path.join(main_path_pi, "AWI-ESM-1-1-LR")


#MH
mh_awi_path = os.path.join(main_path_mh, "AWI-ESM-1-1-LR")
cesm_path = os.path.join(main_path_mh, "CESM2")
ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
hadGEM_path = os.path.join(main_path_mh, "HadGEM3-GC31-LL")
ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")

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


mh_awi_sahel_anomaly, mh_lat_cross_section = read_tp_from_path(path_pi=pi_awi_path, path_main=mh_awi_path)