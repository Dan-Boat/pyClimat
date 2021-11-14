#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 14:06:24 2021

@author: dboateng
"""
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *


module_output_main_path = "/home/dboateng/Model_output_pst"

exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_150 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"

years= "1003_1017"
period = "1m"


# reading dataset
control_data, control_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period)
east_150_data, east_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_150, years=years,
                                                  period=period)
east_0_data, east_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_0, years=years,
                                                  period=period)

lapse_rate = -0.00612

temp = control_data["temp2"] - 273.15
geosp = control_data["geosp"] / 9.81

non_adiatic = temp - (geosp*lapse_rate)

temp_03 = east_150_data["temp2"] - 273.15
geosp_03 = east_150_data["geosp"] / 9.81

non_adiatic_03 = temp_03 - (geosp_03*lapse_rate)

alt_03 = compute_lterm_mean(data=non_adiatic_03, time= "annual")
alt = compute_lterm_mean(data=non_adiatic, time="annual")

alt_diff= compute_lterm_diff(data_control=non_adiatic , data_main=non_adiatic_03 , time="annual", )#season_calendar="standard")

plot_annual_mean(variable="Temperature", data_alt=alt_diff, cmap=RdBu_r, units="°C", vmax=1.5, vmin=-1.5, domain="Europe", 
                 levels=22, level_ticks=11,)