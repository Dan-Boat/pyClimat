#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 18:29:56 2021

@author: dboateng
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)


import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name_aw100e100 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_aw100e0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_aw100e200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"

# west set up 
exp_name_aw200e100 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
exp_name_aw200e0 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
exp_name_aw200e200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"

# for supplementary (same but for annual)
years= "1003_1017"
period = "1m"


# reading dataset
aw100e100_data, aw100e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e100, years=years,
                                                  period=period)
aw100e0_data, aw100e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e0, years=years,
                                                  period=period)
aw100e200_data, aw100e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e200, years=years,
                                                  period=period)
aw200e100_data, aw200e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e100, years=years,
                                                  period=period)
aw200e0_data, aw200e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e0, years=years,
                                                  period=period)
aw200e200_data, aw200e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e200, years=years,
                                                  period=period)




#  extracting variables
 
#aw100e100
temp2_aw100e100 = extract_var(Dataset=aw100e100_data , varname="temp2", units="°C")
prec_aw100e100 = extract_var(Dataset= aw100e100_data , varname="prec", units="mm/month")
d18op_aw100e100 = extract_var(Dataset=aw100e100_data , varname="d18op", units="per mil", Dataset_wiso= aw100e100_wiso)
u10_aw100e100 = extract_var(Dataset=aw100e100_data , varname="u10")
v10_aw100e100 = extract_var(Dataset=aw100e100_data , varname="v10")


#aw100e0
temp2_aw100e0 = extract_var(Dataset=aw100e0_data , varname="temp2", units="°C")
prec_aw100e0 = extract_var(Dataset= aw100e0_data , varname="prec", units="mm/month")
d18op_aw100e0 = extract_var(Dataset=aw100e0_data , varname="d18op", units="per mil", Dataset_wiso= aw100e0_wiso)
u10_aw100e0 = extract_var(Dataset=aw100e0_data , varname="u10")
v10_aw100e0 = extract_var(Dataset=aw100e0_data , varname="v10")


#aw100e200
temp2_aw100e200= extract_var(Dataset=aw100e200_data , varname="temp2", units="°C")
prec_aw100e200 = extract_var(Dataset= aw100e200_data , varname="prec", units="mm/month")
d18op_aw100e200 = extract_var(Dataset=aw100e200_data , varname="d18op", units="per mil", Dataset_wiso= aw100e200_wiso)
u10_aw100e200 = extract_var(Dataset=aw100e200_data , varname="u10")
v10_aw100e200 = extract_var(Dataset=aw100e200_data , varname="v10")


#aw200e100
temp2_aw200e100 = extract_var(Dataset=aw200e100_data , varname="temp2", units="°C")
prec_aw200e100 = extract_var(Dataset= aw200e100_data , varname="prec", units="mm/month")
d18op_aw200e100 = extract_var(Dataset=aw200e100_data , varname="d18op", units="per mil", Dataset_wiso= aw200e100_wiso)
u10_aw200e100 = extract_var(Dataset=aw200e100_data , varname="u10")
v10_aw200e100 = extract_var(Dataset=aw200e100_data , varname="v10")


#aw200e0
temp2_aw200e0 = extract_var(Dataset=aw200e0_data , varname="temp2", units="°C")
prec_aw200e0 = extract_var(Dataset= aw200e0_data , varname="prec", units="mm/month")
d18op_aw200e0 = extract_var(Dataset=aw200e0_data , varname="d18op", units="per mil", Dataset_wiso= aw200e0_wiso)
u10_aw200e0 = extract_var(Dataset=aw200e0_data , varname="u10")
v10_aw200e0 = extract_var(Dataset=aw200e0_data , varname="v10")


#aw200e200
temp2_aw200e200 = extract_var(Dataset=aw200e200_data , varname="temp2", units="°C")
prec_aw200e200 = extract_var(Dataset= aw200e200_data , varname="prec", units="mm/month")
d18op_aw200e200 = extract_var(Dataset=aw200e200_data , varname="d18op", units="per mil", Dataset_wiso= aw200e200_wiso)
u10_aw200e200 = extract_var(Dataset=aw200e200_data , varname="u10")
v10_aw200e200 = extract_var(Dataset=aw200e200_data , varname="v10")



# compute long-term means for control experiment

# seasonal and annual lterm means

#aw100e100
temp2_aw100e100_slt = compute_lterm_mean(data=temp2_aw100e100 , time="season", season_calendar="standard")
prec_aw100e100_slt = compute_lterm_mean(data=prec_aw100e100, time="season", season_calendar="standard")
d18op_aw100e100_slt = compute_lterm_mean(data=d18op_aw100e100, time="season", season_calendar="standard")
u10_aw100e100_slt = compute_lterm_mean(data=u10_aw100e100, time="season", season_calendar="standard")
v10_aw100e100_slt = compute_lterm_mean(data=v10_aw100e100, time="season", season_calendar="standard")


#aw100e0
temp2_aw100e0_slt = compute_lterm_mean(data=temp2_aw100e0 , time="season", season_calendar="standard")
prec_aw100e0_slt = compute_lterm_mean(data=prec_aw100e0, time="season", season_calendar="standard")
d18op_aw100e0_slt = compute_lterm_mean(data=d18op_aw100e0, time="season", season_calendar="standard")
u10_aw100e0_slt = compute_lterm_mean(data=u10_aw100e0, time="season", season_calendar="standard")
v10_aw100e0_slt = compute_lterm_mean(data=v10_aw100e0, time="season", season_calendar="standard")



#aw100e200
temp2_aw100e200_slt = compute_lterm_mean(data=temp2_aw100e200 , time="season", season_calendar="standard")
prec_aw100e200_slt = compute_lterm_mean(data=prec_aw100e200, time="season", season_calendar="standard")
d18op_aw100e200_slt = compute_lterm_mean(data=d18op_aw100e200, time="season", season_calendar="standard")
u10_aw100e200_slt = compute_lterm_mean(data=u10_aw100e200, time="season", season_calendar="standard")
v10_aw100e200_slt = compute_lterm_mean(data=v10_aw100e200, time="season", season_calendar="standard")



#aw200e100
temp2_aw200e100_slt = compute_lterm_mean(data=temp2_aw200e100 , time="season", season_calendar="standard")
prec_aw200e100_slt = compute_lterm_mean(data=prec_aw200e100, time="season", season_calendar="standard")
d18op_aw200e100_slt = compute_lterm_mean(data=d18op_aw200e100, time="season", season_calendar="standard")
u10_aw200e100_slt = compute_lterm_mean(data=u10_aw200e100, time="season", season_calendar="standard")
v10_aw200e100_slt = compute_lterm_mean(data=v10_aw200e100, time="season", season_calendar="standard")


#aw200e0
temp2_aw200e0_slt = compute_lterm_mean(data=temp2_aw200e0 , time="season", season_calendar="standard")
prec_aw200e0_slt = compute_lterm_mean(data=prec_aw200e0, time="season", season_calendar="standard")
d18op_aw200e0_slt = compute_lterm_mean(data=d18op_aw200e0, time="season", season_calendar="standard")
u10_aw200e0_slt = compute_lterm_mean(data=u10_aw200e0, time="season", season_calendar="standard")
v10_aw200e0_slt = compute_lterm_mean(data=v10_aw200e0, time="season", season_calendar="standard")



#aw200e200
temp2_aw200e200_slt = compute_lterm_mean(data=temp2_aw200e200 , time="season", season_calendar="standard")
prec_aw200e200_slt = compute_lterm_mean(data=prec_aw200e200, time="season", season_calendar="standard")
d18op_aw200e200_slt = compute_lterm_mean(data=d18op_aw200e200, time="season", season_calendar="standard")
u10_aw200e200_slt = compute_lterm_mean(data=u10_aw200e200, time="season", season_calendar="standard")
v10_aw200e200_slt = compute_lterm_mean(data=v10_aw200e200, time="season", season_calendar="standard")

