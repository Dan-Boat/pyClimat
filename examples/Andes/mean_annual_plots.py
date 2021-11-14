#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 11:06:11 2021

@author: dboateng
This script generates plots for the long-term annual means of surface uplift experiments conducted on the Andes
(Andes 100%, 75%, 50%, 25% and 0%)
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

#reading data

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

# compute long-term mean or reshape if annual values are loaded

e005_temp2_alt = compute_lterm_mean(data=e005_temp2, time="annual")
e005_prec_alt = compute_lterm_mean(data=e005_prec, time="annual")

e016_temp2_alt = compute_lterm_mean(data=e016_temp2, time="annual")
e016_prec_alt = compute_lterm_mean(data=e016_prec, time="annual")

e013_temp2_alt = compute_lterm_mean(data=e013_temp2, time="annual")
e013_prec_alt = compute_lterm_mean(data=e013_prec, time="annual")

e018_temp2_alt = compute_lterm_mean(data=e018_temp2, time="annual")
e018_prec_alt = compute_lterm_mean(data=e018_prec, time="annual")

e020_temp2_alt = compute_lterm_mean(data=e020_temp2, time="annual")
e020_prec_alt = compute_lterm_mean(data=e020_prec, time="annual")

#visualising variables and saving
projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

#e005-100%
plot_annual_mean(variable="Temperature", data_alt=e005_temp2_alt , cmap=RdBu_r, units="°C", vmax=30, vmin=-15, domain="South America", 
                 levels=22, level_ticks=10, title="Andes 100% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005_temp2_annual_mean.png")

plot_annual_mean(variable="Precipitation", data_alt=e005_prec_alt , cmap=Blues, units="mm/a", vmax=4000, vmin=0, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 100% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005_prec_annual_mean.png")

#e016-75%
plot_annual_mean(variable="Temperature", data_alt=e016_temp2_alt , cmap=RdBu_r, units="°C", vmax=30, vmin=-15, domain="South America", 
                 levels=22, level_ticks=10, title="Andes 75% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e016_temp2_annual_mean.png")

plot_annual_mean(variable="Precipitation", data_alt=e016_prec_alt , cmap=Blues, units="mm/a", vmax=4000, vmin=0, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 75% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e016_prec_annual_mean.png")

#e013-50%
plot_annual_mean(variable="Temperature", data_alt=e013_temp2_alt , cmap=RdBu_r, units="°C", vmax=30, vmin=-15, domain="South America", 
                 levels=22, level_ticks=10, title="Andes 50% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e013_temp2_annual_mean.png")

plot_annual_mean(variable="Precipitation", data_alt=e013_prec_alt , cmap=Blues, units="mm/a", vmax=4000, vmin=0, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 50% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e013_prec_annual_mean.png")

#e018-25%
plot_annual_mean(variable="Temperature", data_alt=e013_temp2_alt , cmap=RdBu_r, units="°C", vmax=30, vmin=-15, domain="South America", 
                 levels=22, level_ticks=10, title="Andes 25% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e018_temp2_annual_mean.png")

plot_annual_mean(variable="Precipitation", data_alt=e018_prec_alt , cmap=Blues, units="mm/a", vmax=4000, vmin=0, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 25% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e018_prec_annual_mean.png")

#e020-25%
plot_annual_mean(variable="Temperature", data_alt=e020_temp2_alt , cmap=RdBu_r, units="°C", vmax=30, vmin=-15, domain="South America", 
                 levels=22, level_ticks=10, title="Andes 0% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e020_temp2_annual_mean.png")

plot_annual_mean(variable="Precipitation", data_alt=e020_prec_alt , cmap=Blues, units="mm/a", vmax=4000, vmin=0, domain= "South America", 
                 levels=22, level_ticks=6,title="Andes 0% [Long-term annual mean]", path_to_store=path_to_store, output_format="png",
                 output_name= "e020_prec_annual_mean.png")