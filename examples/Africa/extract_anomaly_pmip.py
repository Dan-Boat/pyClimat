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
from pyClimat.analysis import compute_lterm_mean, extract_var, extract_profile

from pyClimat.plots import plot_annual_mean
from pyClimat.plot_utils import *


# define path 
path_to_plots = "/home/dboateng/Python_scripts/ClimatPackage_repogit/examples/Africa/plots"

main_path_mh = "/media/dboateng/Boateng/esd02/data/PMIP4/MidHolocene/postprocessed"

awi_path = os.path.join(main_path_mh, "AWI.AWI-ESM-1-1-LR")
cesm_path = os.path.join(main_path_mh, "CESM2")
ec_earth_path = os.path.join(main_path_mh, "EC-Earth3-LR")
giss_path = os.path.join(main_path_mh, "GISS-E2-1-G")
hadGEM_path = os.path.join(main_path_mh, "HadGEM3-GC31-LL")
ipsl_path = os.path.join(main_path_mh, "IPSL-CM6A-LR")
miroc_path = os.path.join(main_path_mh, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_mh, "MPI-ESM1-2-LR")

def read_tp_from_path(path=None, echam=False, data=None):
    
    if echam == True:
        
        pr = extract_var(Dataset=data, varname="prec", units="mm/month")
        
    else:
        pr = read_from_path(path, "pr_1m_lterm.nc", varname="pr", decode=True) *60*60*24*30  #mm/month
    
    pr_alt = compute_lterm_mean(data=pr, time="month", month="JJAS")
    
    minlat = 0
    maxlat = 35
    minlon = -20
    maxlon = 30
    
    
    pr_alt_lat_cross_section = extract_profile(data=pr_alt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lat")
    
    return pr_alt_lat_cross_section

awi_pr_alt_cross_section = read_tp_from_path(awi_path)