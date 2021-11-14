#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 11:07:12 2021

@author: dboateng
his script generates plots for the long-term annual difference of surface uplift experiments conducted on the Andes
(Andes 100% - {75%, 50%, 25% and 0%})
"""
# importimng package 

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

#setting paths 
module_output_main_path = "/home/dboateng/Model_output_pst"

e005 = "e005_hpc-bw_e5w2.3_PD_t63l31.1m"
e016 = "e016_hpc-bw_e5w2.3_PD_sam075_t63l31.1m"
e013 = "e013_hpc-bw_e5w2.3_PD_sam050_t63l31.1m"
e018 = "e018_hpc-bw_e5w2.3_PD_sam025_t63l31.1m"
e020 = "e020_hpc-bw_e5w2.3_PD_sam000_t63l31.1m"
years = "1979_1994"


e005_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name=e005, years=years, read_wiso=False, period="1a")
e016_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name=e016, years=years, read_wiso=False, period="1a")
e013_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name=e013, years=years, read_wiso=False, period="1a")
e018_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name=e018, years=years, read_wiso=False, period="1a")
e020_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name=e020, years=years, read_wiso=False, period="1a")

#extracting temp and prec
e005_temp2 = extract_var(Dataset=e005_data , varname="temp2", units="°C")
e005_prec = extract_var(Dataset= e005_data , varname="prec", units="mm/a")

e016_temp2 = extract_var(Dataset=e016_data , varname="temp2", units="°C")
e016_prec = extract_var(Dataset= e016_data , varname="prec", units="mm/a")

e013_temp2 = extract_var(Dataset=e013_data , varname="temp2", units="°C")
e013_prec = extract_var(Dataset= e013_data , varname="prec", units="mm/a")

e018_temp2 = extract_var(Dataset=e018_data , varname="temp2", units="°C")
e018_prec = extract_var(Dataset= e018_data , varname="prec", units="mm/a")

e020_temp2 = extract_var(Dataset=e020_data , varname="temp2", units="°C")
e020_prec = extract_var(Dataset= e020_data , varname="prec", units="mm/a")

# computing long-term difference
e016_temp2_diff = compute_lterm_diff(data_control=e016_temp2 , data_main=e005_temp2, time="annual")
e016_prec_diff = compute_lterm_diff(data_control=e016_prec , data_main=e005_prec, time="annual")

e013_temp2_diff = compute_lterm_diff(data_control=e013_temp2 , data_main=e005_temp2, time="annual")
e013_prec_diff = compute_lterm_diff(data_control=e013_prec , data_main=e005_prec, time="annual")

e018_temp2_diff = compute_lterm_diff(data_control=e018_temp2 , data_main=e005_temp2, time="annual")
e018_prec_diff = compute_lterm_diff(data_control=e018_prec , data_main=e005_prec, time="annual")

e020_temp2_diff = compute_lterm_diff(data_control=e020_temp2 , data_main=e005_temp2, time="annual")
e020_prec_diff = compute_lterm_diff(data_control=e020_prec , data_main=e005_prec, time="annual")

#visualising variables and saving
projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

# 100 -75 %
plot_annual_mean(variable="Temperature", data_alt=e016_temp2_diff , cmap=RdBu_r, units="°C", vmax=15, vmin=-15, domain="South America", 
                 levels=22, level_ticks=11, title="Andes 100% - Andes 75% [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-e016_annual_diff_temp.png")

plot_annual_mean(variable="Precipitation", data_alt=e016_prec_diff , cmap=RdBu, units="mm/a", vmax=2000, vmin=-2000, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 100% - Andes 75% [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-e016_annual_diff_prec.png")

#100 -50 %
plot_annual_mean(variable="Temperature", data_alt=e013_temp2_diff , cmap=RdBu_r, units="°C", vmax=15, vmin=-15, domain="South America", 
                 levels=22, level_ticks=11, title="Andes 100% - Andes 50% [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-e013_annual_diff_temp.png")

plot_annual_mean(variable="Precipitation", data_alt=e013_prec_diff , cmap=RdBu, units="mm/a", vmax=2000, vmin=-2000, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 100% - Andes 50% [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-e013_annual_diff_prec.png")

#100-25 %
plot_annual_mean(variable="Temperature", data_alt=e018_temp2_diff , cmap=RdBu_r, units="°C", vmax=15, vmin=-15, domain="South America", 
                 levels=22, level_ticks=11, title="Andes 100% - Andes 25% [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-e018_annual_diff_temp.png")

plot_annual_mean(variable="Precipitation", data_alt=e018_prec_diff , cmap=RdBu, units="mm/a", vmax=2000, vmin=-2000, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 100% - Andes 25% [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-e018_annual_diff_prec.png")

#100 -0 %
plot_annual_mean(variable="Temperature", data_alt=e020_temp2_diff , cmap=RdBu_r, units="°C", vmax=15, vmin=-15, domain="South America", 
                 levels=22, level_ticks=11, title="Andes 100% - Andes 0% [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-e020_annual_diff_temp.png")

plot_annual_mean(variable="Precipitation", data_alt=e020_prec_diff , cmap=RdBu, units="mm/a", vmax=2000, vmin=-2000, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 100% - Andes 0% [long-term annual difference]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005-e020_annual_diff_prec.png")

