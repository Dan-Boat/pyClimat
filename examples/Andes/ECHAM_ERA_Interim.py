#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 16:28:01 2021

@author: dboateng
"""
# importimng package 

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

t2m_path = "/home/dboateng/Model_output_pst/era/t2m_monthly.nc"
tp_path = "/home/dboateng/Model_output_pst/era/tp_monthly.nc"

# ERA-interim 
t2m = read_ERA_processed(path=t2m_path, varname="t2m") - 273.15 #°C
tp = read_ERA_processed(path=tp_path, varname="tp") * 1000 * 12 #mm/a

t2m = t2m.rename({"longitude": "lon", "latitude":"lat"})

#setting paths 
module_output_main_path = "/home/dboateng/Model_output_pst"
years = "1979_1994"

e005 = "e005_hpc-bw_e5w2.3_PD_t63l31.1m"

#reading 
e005_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name=e005, years=years, read_wiso=False, period="1a")

t2m = t2m.interp(lat=e005_data.lat).interp(lon=e005_data.lon)
tp = tp.interp(lat=e005_data.lat).interp(lon=e005_data.lon)

#extracting temp and prec
e005_temp2 = extract_var(Dataset=e005_data , varname="temp2", units="°C")
e005_prec = extract_var(Dataset= e005_data , varname="prec", units="mm/a")

# computing long-term difference
e005_temp2_diff = compute_lterm_diff(data_control=t2m , data_main=e005_temp2, time="annual")
e005_prec_diff = compute_lterm_diff(data_control=tp , data_main=e005_prec, time="annual")

#visualising variables and saving
projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

plot_annual_mean(variable="Temperature", data_alt=e005_temp2_diff , cmap=RdBu_r, units="°C", vmax=15, vmin=-15, domain="South America", 
                 levels=22, level_ticks=11, title="Andes 100% - ERA-Interim [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-era_annual_diff_temp.png")

plot_annual_mean(variable="Precipitation", data_alt=e005_prec_diff , cmap=RdBu, units="mm/a", vmax=2000, vmin=-2000, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 100% - ERA-Interim [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-era_annual_diff_prec.png")


