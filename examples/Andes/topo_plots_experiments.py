#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 11:07:38 2021

@author: dboateng
This script generates plots for the topographic configurations used for sensitivity experiments of Andes
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
filename = "T63_jan_surf.nc"

#reading data
e005_data = read_ECHAM_input(main_path=module_output_main_path , exp_name=e005, filename=filename,
                             read_var=True, varname="OROMEA")
e016_data = read_ECHAM_input(main_path=module_output_main_path , exp_name=e016, filename=filename,
                             read_var=True, varname="OROMEA")
e013_data = read_ECHAM_input(main_path=module_output_main_path , exp_name=e013, filename=filename,
                             read_var=True, varname="OROMEA")
e018_data = read_ECHAM_input(main_path=module_output_main_path , exp_name=e018, filename=filename,
                             read_var=True, varname="OROMEA")
e020_data = read_ECHAM_input(main_path=module_output_main_path , exp_name=e020, filename=filename,
                             read_var=True, varname="OROMEA")

#visualisation
projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

plot_annual_mean(variable="Elevation", data_alt=e005_data , cmap=Greys, units="m", vmax=3500, vmin=0, domain="South America", 
                 levels=22, level_ticks=6, title="Andes 100% [Elevation]", path_to_store=path_to_store, output_format="png",
                 output_name= "e005_elevation.png")
plot_annual_mean(variable="Elevation", data_alt=e016_data , cmap=Greys, units="m", vmax=3500, vmin=0, domain="South America", 
                 levels=22, level_ticks=6, title="Andes 75% [Elevation]", path_to_store=path_to_store, output_format="png",
                 output_name= "e016_elevation.png")
plot_annual_mean(variable="Elevation", data_alt=e013_data , cmap=Greys, units="m", vmax=3500, vmin=0, domain="South America", 
                 levels=22, level_ticks=6, title="Andes 50% [Elevation]", path_to_store=path_to_store, output_format="png",
                 output_name= "e013_elevation.png")
plot_annual_mean(variable="Elevation", data_alt=e018_data , cmap=Greys, units="m", vmax=3500, vmin=0, domain="South America", 
                 levels=22, level_ticks=6, title="Andes 25% [Elevation]", path_to_store=path_to_store, output_format="png",
                 output_name= "e018_elevation.png")
plot_annual_mean(variable="Elevation", data_alt=e020_data , cmap=Greys, units="m", vmax=3500, vmin=0, domain="South America", 
                 levels=22, level_ticks=6, title="Andes 0% [Elevation]", path_to_store=path_to_store, output_format="png",
                 output_name= "e020_elevation.png")