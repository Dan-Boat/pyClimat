#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 15:23:06 2021

@author: dboateng
This is example script for using Climat in analysing and visualising ECHAM module output (isotopic lapse rate)
The script contains directories of module outputs and path to save plot
Note: it is structured solely for the personal needs of the author, therefore, it must be adapted advisably.
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the function

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *


# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"

exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east_200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_Alps_150 = "a005_hpc-bw_e5w2.3_t159_PI_Alps_150_t159l31.6h"
exp_name_east_150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"
exp_name_Alps_0 = "a006_hpc-bw_e5w2.3_t159_PI_Alps_0_t159l31.6h"
#exp_name_Alps_0 = "t015_dkrz-mistral_e5w2.3_PI_Alps0_t159l31.6h"
exp_name_east_50 = "a004_hpc-bw_e5w2.3_t159_PI_Alps_east_50_t159l31.6h"
exp_name_Alps_200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"


years= "1003_1017"
period = "1m"


# reading dataset
control_data, control_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period)
east_0_data, east_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_0, years=years,
                                                  period=period)
east_150_data, east_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_150, years=years,
                                                  period=period)

east_50_data, east_50_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_50, years=years,
                                                  period=period)
east_200_data, east_200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_200, years=years,
                                                  period=period)


#extracting variables and computing long-term means

#control
d18op = extract_var(Dataset=control_data , varname="d18op", units="per mil", Dataset_wiso= control_wiso)
elev = extract_var(Dataset=control_data , varname="elev", units="m")

d18op_slt = compute_lterm_mean(data=d18op, time="season", season_calendar="standard")
elev_slt = compute_lterm_mean(data=elev, time="season", season_calendar="standard")

d18op_alt = compute_lterm_mean(data=d18op, time="annual",)
elev_alt = compute_lterm_mean(data=elev, time="annual",)


#east_0
d18op_east_0 = extract_var(Dataset=east_0_data , varname="d18op", units="per mil", Dataset_wiso= east_0_wiso)
elev_east_0 = extract_var(Dataset=east_0_data , varname="elev", units="m")

east_0_d18op_slt = compute_lterm_mean(data=d18op_east_0, time="season", season_calendar="standard")
east_0_elev_slt = compute_lterm_mean(data=elev_east_0, time="season", season_calendar="standard")

east_0_d18op_alt = compute_lterm_mean(data=d18op_east_0, time="annual",)
east_0_elev_alt = compute_lterm_mean(data=elev_east_0, time="annual",)



#east_200
d18op_east_200 = extract_var(Dataset=east_200_data , varname="d18op", units="per mil", Dataset_wiso= east_200_wiso)
elev_east_200 = extract_var(Dataset=east_200_data , varname="elev", units="m")

east_200_d18op_slt = compute_lterm_mean(data=d18op_east_200, time="season", season_calendar="standard")
east_200_elev_slt = compute_lterm_mean(data=elev_east_200, time="season", season_calendar="standard")

east_200_d18op_alt = compute_lterm_mean(data=d18op_east_200, time="annual",)
east_200_elev_alt = compute_lterm_mean(data=elev_east_200, time="annual",)



# defining coordinates 

maxlat_west, minlat_west, maxlon_west, minlon_west = 47, 44, 8, 1
maxlat_south, minlat_south, maxlon_south, minlon_south = 47, 43, 15, 7.5
maxlat_north, minlat_north, maxlon_north, minlon_north = 50, 46.5, 16, 5

# extracting transects 


# control 

elev_control_north_slt = extract_transect(data=elev_slt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, Dataset=control_data)
elev_control_west_slt = extract_transect(data=elev_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=control_data)
elev_control_south_slt = extract_transect(data=elev_slt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=control_data)

elev_control_north_alt = extract_transect(data=elev_alt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, Dataset=control_data)
elev_control_west_alt = extract_transect(data=elev_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=control_data)
elev_control_south_alt = extract_transect(data=elev_alt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=control_data)

d18op_control_north_slt = extract_transect(data=d18op_slt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=control_data)
d18op_control_west_slt = extract_transect(data=d18op_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=control_data)
d18op_control_south_slt = extract_transect(data=d18op_slt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=control_data)

d18op_control_north_alt = extract_transect(data=d18op_alt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=control_data)
d18op_control_west_alt = extract_transect(data=d18op_alt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=control_data)
d18op_control_south_alt = extract_transect(data=d18op_alt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=control_data)




#east 200

elev_east_200_north_slt = extract_transect(data=east_200_elev_slt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=east_200_data)
elev_east_200_west_slt = extract_transect(data=east_200_elev_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=east_200_data)
elev_east_200_south_slt = extract_transect(data=east_200_elev_slt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=east_200_data)

elev_east_200_north_alt = extract_transect(data=east_200_elev_alt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=east_200_data)
elev_east_200_west_alt = extract_transect(data=east_200_elev_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=east_200_data)
elev_east_200_south_alt = extract_transect(data=east_200_elev_alt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=east_200_data)

d18op_east_200_north_slt = extract_transect(data=east_200_d18op_slt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=east_200_data)
d18op_east_200_west_slt = extract_transect(data=east_200_d18op_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=east_200_data)
d18op_east_200_south_slt = extract_transect(data=east_200_d18op_slt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=east_200_data)

d18op_east_200_north_alt = extract_transect(data=east_200_d18op_alt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=east_200_data)
d18op_east_200_west_alt = extract_transect(data=east_200_d18op_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=east_200_data)
d18op_east_200_south_alt = extract_transect(data=east_200_d18op_alt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=east_200_data)




#east 0

elev_east_0_north_slt = extract_transect(data=east_0_elev_slt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=east_0_data)
elev_east_0_west_slt = extract_transect(data=east_0_elev_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=east_0_data)
elev_east_0_south_slt = extract_transect(data=east_0_elev_slt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=east_0_data)

elev_east_0_north_alt = extract_transect(data=east_0_elev_alt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=east_0_data)
elev_east_0_west_alt = extract_transect(data=east_0_elev_alt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=east_0_data)
elev_east_0_south_alt = extract_transect(data=east_0_elev_alt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=east_0_data)

d18op_east_0_north_slt = extract_transect(data=east_0_d18op_slt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=east_0_data)
d18op_east_0_west_slt = extract_transect(data=east_0_d18op_slt ,maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=east_0_data)
d18op_east_0_south_slt = extract_transect(data=east_0_d18op_slt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=east_0_data)

d18op_east_0_north_alt = extract_transect(data=east_0_d18op_alt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=east_0_data)
d18op_east_0_west_alt = extract_transect(data=east_0_d18op_alt ,maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=east_0_data)
d18op_east_0_south_alt = extract_transect(data=east_0_d18op_alt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=east_0_data)




season = "JJA"

# performing regression (change DJF to JJA for summer plots)

#Alps 100
control_west_reg_slt, control_west_df_slt = linregression(data_x=elev_control_west_slt , data_y=d18op_control_west_slt, season=season, return_yhat=True)
control_north_reg_slt, control_north_df_slt = linregression(data_x=elev_control_north_slt , data_y=d18op_control_north_slt, season=season, return_yhat=True)
control_south_reg_slt, control_south_df_slt = linregression(data_x=elev_control_south_slt , data_y=d18op_control_south_slt, season=season, return_yhat=True)

control_west_reg_alt, control_west_df_alt = linregression(data_x=elev_control_west_alt , data_y=d18op_control_west_alt, return_yhat=True)
control_north_reg_alt, control_north_df_alt = linregression(data_x=elev_control_north_alt , data_y=d18op_control_north_alt, return_yhat=True)
control_south_reg_alt, control_south_df_alt = linregression(data_x=elev_control_south_alt , data_y=d18op_control_south_alt, return_yhat=True)




#Alps_east200
east_200_west_reg_slt, east_200_west_df_slt = linregression(data_x=elev_east_200_west_slt , data_y=d18op_east_200_west_slt, season=season, return_yhat=True)
east_200_north_reg_slt, east_200_north_df_slt = linregression(data_x=elev_east_200_north_slt , data_y=d18op_east_200_north_slt, season=season, return_yhat=True)
east_200_south_reg_slt, east_200_south_df_slt = linregression(data_x=elev_east_200_south_slt , data_y=d18op_east_200_south_slt, season=season, return_yhat=True)

east_200_west_reg_alt, east_200_west_df_alt = linregression(data_x=elev_east_200_west_alt , data_y=d18op_east_200_west_alt,  return_yhat=True)
east_200_north_reg_alt, east_200_north_df_alt = linregression(data_x=elev_east_200_north_alt , data_y=d18op_east_200_north_alt, return_yhat=True)
east_200_south_reg_alt, east_200_south_df_alt = linregression(data_x=elev_east_200_south_alt , data_y=d18op_east_200_south_alt, return_yhat=True)


#Alps east0
east_0_west_reg_slt, east_0_west_df_slt = linregression(data_x=elev_east_0_west_slt , data_y=d18op_east_0_west_slt, season=season, return_yhat=True)
east_0_north_reg_slt, east_0_north_df_slt = linregression(data_x=elev_east_0_north_slt , data_y=d18op_east_0_north_slt, season=season, return_yhat=True)
east_0_south_reg_slt, east_0_south_df_slt = linregression(data_x=elev_east_0_south_slt , data_y=d18op_east_0_south_slt, season=season, return_yhat=True)

east_0_west_reg_alt, east_0_west_df_alt = linregression(data_x=elev_east_0_west_alt , data_y=d18op_east_0_west_alt, return_yhat=True)
east_0_north_reg_alt, east_0_north_df_alt = linregression(data_x=elev_east_0_north_alt , data_y=d18op_east_0_north_alt,  return_yhat=True)
east_0_south_reg_alt, east_0_south_df_alt = linregression(data_x=elev_east_0_south_alt , data_y=d18op_east_0_south_alt, return_yhat=True)




# calculation for bulk experiments
def calculate_bulk_exp():
    Alps_150_data, Alps_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_150, years=years,
                                                      period=period)
    Alps_0_data, Alps_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_0, years=years,
                                                      period=period)
    Alps_200_data, Alps_200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_200, years=years,
                                                      period=period)
    #Alps_0
    d18op_Alps_0 = extract_var(Dataset=Alps_0_data , varname="d18op", units="per mil", Dataset_wiso= Alps_0_wiso)
    elev_Alps_0 = extract_var(Dataset=Alps_0_data , varname="elev", units="m")
    
    Alps_0_d18op_slt = compute_lterm_mean(data=d18op_Alps_0, time="season", season_calendar="standard")
    Alps_0_elev_slt = compute_lterm_mean(data=elev_Alps_0, time="season", season_calendar="standard")
    
    Alps_0_d18op_alt = compute_lterm_mean(data=d18op_Alps_0, time="annual",)
    Alps_0_elev_alt = compute_lterm_mean(data=elev_Alps_0, time="annual",)
    
    
    #Alps_200
    d18op_Alps_200 = extract_var(Dataset=Alps_200_data , varname="d18op", units="per mil", Dataset_wiso= Alps_200_wiso)
    elev_Alps_200 = extract_var(Dataset=Alps_200_data , varname="elev", units="m")
    
    Alps_200_d18op_slt = compute_lterm_mean(data=d18op_Alps_200, time="season", season_calendar="standard")
    Alps_200_elev_slt = compute_lterm_mean(data=elev_Alps_200, time="season", season_calendar="standard")
    
    Alps_200_d18op_alt = compute_lterm_mean(data=d18op_Alps_200, time="annual",)
    Alps_200_elev_alt = compute_lterm_mean(data=elev_Alps_200, time="annual",)
    
    maxlat_west, minlat_west, maxlon_west, minlon_west = 47, 44, 8, 1
    maxlat_south, minlat_south, maxlon_south, minlon_south = 47, 43, 15, 7.5
    maxlat_north, minlat_north, maxlon_north, minlon_north = 50, 46.5, 16, 5
    
    #Alps 200
    
    elev_200_north_slt = extract_transect(data=Alps_200_elev_slt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=Alps_200_data)
    elev_200_west_slt = extract_transect(data=Alps_200_elev_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=Alps_200_data)
    elev_200_south_slt = extract_transect(data=Alps_200_elev_slt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=Alps_200_data)
    
    elev_200_north_alt = extract_transect(data=Alps_200_elev_alt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=Alps_200_data)
    elev_200_west_alt = extract_transect(data=Alps_200_elev_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=Alps_200_data)
    elev_200_south_alt = extract_transect(data=Alps_200_elev_alt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=Alps_200_data)
    
    d18op_200_north_slt = extract_transect(data=Alps_200_d18op_slt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=Alps_200_data)
    d18op_200_west_slt = extract_transect(data=Alps_200_d18op_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=Alps_200_data)
    d18op_200_south_slt = extract_transect(data=Alps_200_d18op_slt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=Alps_200_data)
    
    d18op_200_north_alt = extract_transect(data=Alps_200_d18op_alt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=Alps_200_data)
    d18op_200_west_alt = extract_transect(data=Alps_200_d18op_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=Alps_200_data)
    d18op_200_south_alt = extract_transect(data=Alps_200_d18op_alt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=Alps_200_data)
    
    # Alps 0
    elev_0_north_slt = extract_transect(data=Alps_0_elev_slt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=Alps_0_data)
    elev_0_west_slt = extract_transect(data=Alps_0_elev_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=Alps_0_data)
    elev_0_south_slt = extract_transect(data=Alps_0_elev_slt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=Alps_0_data)
    
    elev_0_north_alt = extract_transect(data=Alps_0_elev_alt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=Alps_0_data)
    elev_0_west_alt = extract_transect(data=Alps_0_elev_alt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=Alps_0_data)
    elev_0_south_alt = extract_transect(data=Alps_0_elev_alt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=Alps_0_data)
    
    d18op_0_north_slt = extract_transect(data=Alps_0_d18op_slt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=Alps_0_data)
    d18op_0_west_slt = extract_transect(data=Alps_0_d18op_slt ,maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=Alps_0_data)
    d18op_0_south_slt = extract_transect(data=Alps_0_d18op_slt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=Alps_0_data)
    
    d18op_0_north_alt = extract_transect(data=Alps_0_d18op_alt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=Alps_0_data)
    d18op_0_west_alt = extract_transect(data=Alps_0_d18op_alt ,maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=Alps_0_data)
    d18op_0_south_alt = extract_transect(data=Alps_0_d18op_alt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=Alps_0_data)
    
    season = "JJA"
    
    #Alps 200
    Alps_200_west_reg_slt, Alps_200_west_df_slt = linregression(data_x=elev_200_west_slt , data_y=d18op_200_west_slt, season=season, return_yhat=True)
    Alps_200_north_reg_slt, Alps_200_north_df_slt = linregression(data_x=elev_200_north_slt , data_y=d18op_200_north_slt, season=season, return_yhat=True)
    Alps_200_south_reg_slt, Alps_200_south_df_slt = linregression(data_x=elev_200_south_slt , data_y=d18op_200_south_slt, season=season, return_yhat=True)
    
    Alps_200_west_reg_alt, Alps_200_west_df_alt = linregression(data_x=elev_200_west_alt , data_y=d18op_200_west_alt, return_yhat=True)
    Alps_200_north_reg_alt, Alps_200_north_df_alt = linregression(data_x=elev_200_north_alt , data_y=d18op_200_north_alt, return_yhat=True)
    Alps_200_south_reg_alt, Alps_200_south_df_alt = linregression(data_x=elev_200_south_alt , data_y=d18op_200_south_alt, return_yhat=True)
    
    #Alps 0
    Alps_0_west_reg_slt, Alps_0_west_df_slt = linregression(data_x=elev_0_west_slt , data_y=d18op_0_west_slt, season=season, return_yhat=True)
    Alps_0_north_reg_slt, Alps_0_north_df_slt = linregression(data_x=elev_0_north_slt , data_y=d18op_0_north_slt, season=season, return_yhat=True)
    Alps_0_south_reg_slt, Alps_0_south_df_slt = linregression(data_x=elev_0_south_slt , data_y=d18op_0_south_slt, season=season, return_yhat=True)
    
    Alps_0_west_reg_alt, Alps_0_west_df_alt = linregression(data_x=elev_0_west_alt , data_y=d18op_0_west_alt, return_yhat=True)
    Alps_0_north_reg_alt, Alps_0_north_df_alt = linregression(data_x=elev_0_north_alt , data_y=d18op_0_north_alt, return_yhat=True)
    Alps_0_south_reg_alt, Alps_0_south_df_alt = linregression(data_x=elev_0_south_alt , data_y=d18op_0_south_alt, return_yhat=True)



# plotting 
import matplotlib as mpl
path_to_store = os.path.join(module_output_main_path, "plots")
plt.style.use("bmh")
plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
mpl.rc('text', usetex=True)
mpl.rc('font', size=20, family='serif')
mpl.rc('xtick', labelsize=20)
mpl.rc('ytick', labelsize=20)
mpl.rc('legend', fontsize=20)
mpl.rc('axes', labelsize=20)
mpl.rc('lines', linewidth=3)

def plot_all_lapse_rate():


    fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), )
    
    #ax1
    scatter_plot_laspe_rate(ax=ax1, reg_params= control_west_reg_slt , df_x_y_yhat=control_west_df_slt , color=red, marker= "*", label= "west",
                           title="[A] Alps 100 (" + season + ")", xmax=1500, xmin=0,
                            ymax=0, ymin= -12, )
    scatter_plot_laspe_rate(ax=ax1, reg_params= control_north_reg_slt , df_x_y_yhat=control_north_df_slt , color=black, marker= "D", label= "north")
    scatter_plot_laspe_rate(ax=ax1, reg_params= control_south_reg_slt , df_x_y_yhat=control_south_df_slt , color=green, marker= "^", label= "south")
                          
    ax1.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    #ax2
    scatter_plot_laspe_rate(ax=ax2, reg_params= control_west_reg_alt , df_x_y_yhat=control_west_df_alt , color=red, marker= "*", label= "west",
                           title="[A] Alps 100" , xmax=1500, xmin=0,
                            ymax=0, ymin= -12, left_labels=False)
    scatter_plot_laspe_rate(ax=ax2, reg_params= control_north_reg_alt , df_x_y_yhat=control_north_df_alt , color=black, marker= "D", label= "north",
                            left_labels=False)
    scatter_plot_laspe_rate(ax=ax2, reg_params= control_south_reg_alt , df_x_y_yhat=control_south_df_alt , color=green, marker= "^", label= "south",
                            left_labels=False)
    ax2.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    #ax3
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_200_west_reg_slt , df_x_y_yhat=east_200_west_df_slt , color=red, marker= "*", label= "west",
                           title="[A] Alps east 200 (" + season + ")", xmax=3500, xmin=0,
                            ymax=0, ymin= -18, bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_200_north_reg_slt , df_x_y_yhat=east_200_north_df_slt , color=black, marker= "D", label= "north",
                            bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_200_south_reg_slt , df_x_y_yhat=east_200_south_df_slt , color=green, marker= "^", label= "south",
                           bottom_labels=False)
    ax3.legend(frameon=True, fontsize=14, loc="upper right", framealpha=0.5, ncol=1,)
    
    #ax4
    scatter_plot_laspe_rate(ax=ax4, reg_params= east_200_west_reg_alt , df_x_y_yhat=east_200_west_df_alt , color=red, marker= "*", label= "west",
                            title="[A] Alps east 200" , xmax=3500, xmin=0,
                            ymax=0, ymin= -18, bottom_labels=False, left_labels=False)
    scatter_plot_laspe_rate(ax=ax4, reg_params= east_200_north_reg_alt , df_x_y_yhat=east_200_north_df_alt , color=black, marker= "D", label= "north",
                            bottom_labels=False, left_labels=False)
    scatter_plot_laspe_rate(ax=ax4, reg_params= east_200_south_reg_alt , df_x_y_yhat=east_200_south_df_alt , color=green, marker= "^", label= "south",
                            bottom_labels=False, left_labels=False)
    ax4.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    #ax5
    scatter_plot_laspe_rate(ax=ax5, reg_params= Alps_200_west_reg_slt , df_x_y_yhat=Alps_200_west_df_slt , color=red, marker= "*", label= "west",
                            title="[A] Alps 200 (" + season + ")", xmax=3500, xmin=0,
                            ymax=0, ymin= -18, )
    scatter_plot_laspe_rate(ax=ax5, reg_params= Alps_200_north_reg_slt , df_x_y_yhat=Alps_200_north_df_slt , color=black, marker= "D", label= "north",
                            )
    scatter_plot_laspe_rate(ax=ax5, reg_params= Alps_200_south_reg_slt , df_x_y_yhat=Alps_200_south_df_slt , color=green, marker= "^", label= "south",
                           )
    ax5.legend(frameon=True, fontsize=14, loc="upper right", framealpha=0.5, ncol=1,)
    
    #ax6
    scatter_plot_laspe_rate(ax=ax6, reg_params= Alps_200_west_reg_alt , df_x_y_yhat=Alps_200_west_df_alt , color=red, marker= "*", label= "west",
                            title="[A] Alps 200" , xmax=3500, xmin=0,
                            ymax=0, ymin= -18, left_labels=False)
    scatter_plot_laspe_rate(ax=ax6, reg_params= Alps_200_north_reg_alt , df_x_y_yhat=Alps_200_north_df_alt , color=black, marker= "D", label= "north",
                           left_labels=False)
    scatter_plot_laspe_rate(ax=ax6, reg_params= Alps_200_south_reg_alt , df_x_y_yhat=Alps_200_south_df_alt , color=green, marker= "^", label= "south",
                            left_labels=False)
    ax6.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.88, top=0.97, bottom=0.05)
    
    plt.savefig(os.path.join(path_to_store, "fig7.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
    fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), )
    
    #ax1
    scatter_plot_laspe_rate(ax=ax1, reg_params= control_west_reg_slt , df_x_y_yhat=control_west_df_slt , color=red, marker= "*", label= "west",
                           title="[A] Alps 100 (" + season + ")", xmax=1500, xmin=0,
                            ymax=0, ymin= -12, bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax1, reg_params= control_north_reg_slt , df_x_y_yhat=control_north_df_slt , color=black, marker= "D", label= "north",
                            bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax1, reg_params= control_south_reg_slt , df_x_y_yhat=control_south_df_slt , color=green, marker= "^", label= "south",
                            bottom_labels=False)
                          
    ax1.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    #ax2
    scatter_plot_laspe_rate(ax=ax2, reg_params= control_west_reg_alt , df_x_y_yhat=control_west_df_alt , color=red, marker= "*", label= "west",
                           title="[A] Alps 100" , xmax=1500, xmin=0,
                            ymax=0, ymin= -12, left_labels=False, bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax2, reg_params= control_north_reg_alt , df_x_y_yhat=control_north_df_alt , color=black, marker= "D", label= "north",
                            left_labels=False, bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax2, reg_params= control_south_reg_alt , df_x_y_yhat=control_south_df_alt , color=green, marker= "^", label= "south",
                            left_labels=False, bottom_labels=False)
    ax2.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    #ax3
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_west_reg_slt , df_x_y_yhat=east_0_west_df_slt , color=red, marker= "*", label= "west",
                           title="[A] Alps east 0 (" + season + ")", xmax=1500, xmin=0,
                            ymax=0, ymin= -12, bottom_labels=True)
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_north_reg_slt , df_x_y_yhat=east_0_north_df_slt , color=black, marker= "D", label= "north",
                            bottom_labels=True)
    
    
    ax3.legend(frameon=True, fontsize=14, loc="upper right", framealpha=0.5, ncol=1,)
    
    #ax4
    scatter_plot_laspe_rate(ax=ax4, reg_params= east_0_west_reg_alt , df_x_y_yhat=east_0_west_df_alt , color=red, marker= "*", label= "west",
                            title="[A] Alps east 0" , xmax=1500, xmin=0,
                            ymax=0, ymin= -12, bottom_labels=True, left_labels=False)
    scatter_plot_laspe_rate(ax=ax4, reg_params= east_0_north_reg_alt , df_x_y_yhat=east_0_north_df_alt , color=black, marker= "D", label= "north",
                            bottom_labels=True, left_labels=False)
    scatter_plot_laspe_rate(ax=ax4, reg_params= east_0_south_reg_alt , df_x_y_yhat=east_0_south_df_alt , color=green, marker= "^", label= "south",
                            bottom_labels=True, left_labels=False)
    ax4.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    #ax5
    scatter_plot_laspe_rate(ax=ax5, reg_params= Alps_0_west_reg_slt , df_x_y_yhat=Alps_0_west_df_slt , color=red, marker= "*", label= "west",
                            title="[A] Alps 0 (" + season + ")", xmax=600, xmin=0,
                            ymax=0, ymin= -10, )
    scatter_plot_laspe_rate(ax=ax5, reg_params= Alps_0_north_reg_slt , df_x_y_yhat=Alps_0_north_df_slt , color=black, marker= "D", label= "north",
                            )
    scatter_plot_laspe_rate(ax=ax5, reg_params= Alps_0_south_reg_slt , df_x_y_yhat=Alps_0_south_df_slt , color=green, marker= "^", label= "south",
                           )
    ax5.legend(frameon=True, fontsize=14, loc="upper right", framealpha=0.5, ncol=1,)
    
    #ax6
    scatter_plot_laspe_rate(ax=ax6, reg_params= Alps_0_west_reg_alt , df_x_y_yhat=Alps_0_west_df_alt , color=red, marker= "*", label= "west",
                            title="[A] Alps 0" , xmax=600, xmin=0,
                            ymax=0, ymin= -10, left_labels=False)
    scatter_plot_laspe_rate(ax=ax6, reg_params= Alps_0_north_reg_alt , df_x_y_yhat=Alps_0_north_df_alt , color=black, marker= "D", label= "north",
                           left_labels=False)
    scatter_plot_laspe_rate(ax=ax6, reg_params= Alps_0_south_reg_alt , df_x_y_yhat=Alps_0_south_df_alt , color=green, marker= "^", label= "south",
                            left_labels=False)
    ax6.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.88, top=0.97, bottom=0.05)
    
    plt.savefig(os.path.join(path_to_store, "fig8.svg"), format= "svg", bbox_inches="tight", dpi=600)


# ploting per section for different exps. 

def plot_lape_rate_per_section():

    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols = 1, figsize=(8, 15), )
    
    #ax1 (west)
    
    scatter_plot_laspe_rate(ax=ax1, reg_params= control_west_reg_slt , df_x_y_yhat=control_west_df_slt , color=black, marker= "*", label= "Alps 100",
                           title="[A] West", xmax=3500, xmin=0,
                            ymax=0, ymin= -18, bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax1, reg_params= east_0_west_reg_slt , df_x_y_yhat=east_0_west_df_slt , color=red, marker= "D", label= "Alps east 0",
                           bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax1, reg_params= east_200_west_reg_slt , df_x_y_yhat=east_200_west_df_slt , color=green, marker= "^", label= "Alps east 200",
                           bottom_labels=False)
    
    ax1.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    
    #ax2 (north)
    scatter_plot_laspe_rate(ax=ax2, reg_params= control_north_reg_slt , df_x_y_yhat=control_north_df_slt , color=black, marker= "*", label= "Alps 100",
                            bottom_labels=False, xmax=3500, xmin=0, title= "[B] North",
                             ymax=0, ymin= -18,)
    scatter_plot_laspe_rate(ax=ax2, reg_params= east_0_north_reg_slt , df_x_y_yhat=east_0_north_df_slt , color=red, marker= "D", label= "Alps east 0",
                            bottom_labels=False)
    
    scatter_plot_laspe_rate(ax=ax2, reg_params= east_200_north_reg_slt , df_x_y_yhat=east_200_north_df_slt , color=green, marker= "^", label= "Alps east 200",
                            bottom_labels=False)
    
    ax2.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    #ax3 (south)
    scatter_plot_laspe_rate(ax=ax3, reg_params= control_south_reg_slt , df_x_y_yhat=control_south_df_slt , color=black, marker= "*", label= "Alps 100",
                            bottom_labels=True, xmax=3500, xmin=0, title= "[C] South",
                             ymax=0, ymin= -18,)
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_south_reg_slt , df_x_y_yhat=east_0_south_df_slt , color=red, marker= "^", label= "Alps east 0",
                           bottom_labels=True)
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_200_south_reg_slt , df_x_y_yhat=east_200_south_df_slt , color=green, marker= "D", label= "Alps east 200",
                           bottom_labels=True)
    
    ax3.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.88, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_store, "fig7.svg"), format= "svg", bbox_inches="tight", dpi=600)
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols = 1, figsize=(8, 15), )
    
    #ax1 (west)
    
    scatter_plot_laspe_rate(ax=ax1, reg_params= control_west_reg_alt , df_x_y_yhat=control_west_df_alt , color=black, marker= "*", label= "Alps 100",
                           title="[A] West", xmax=3500, xmin=0,
                            ymax=0, ymin= -18, bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax1, reg_params= east_0_west_reg_alt , df_x_y_yhat=east_0_west_df_alt , color=red, marker= "D", label= "Alps east 0",
                           bottom_labels=False)
    scatter_plot_laspe_rate(ax=ax1, reg_params= east_200_west_reg_alt , df_x_y_yhat=east_200_west_df_alt , color=green, marker= "^", label= "Alps east 200",
                           bottom_labels=False)
    
    ax1.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    
    #ax2 (north)
    
    scatter_plot_laspe_rate(ax=ax2, reg_params= control_north_reg_alt , df_x_y_yhat=control_north_df_alt , color=black, marker= "*", label= "Alps 100",
                            bottom_labels=False, xmax=3500, xmin=0, title= "[B] North",
                             ymax=0, ymin= -18,)
    scatter_plot_laspe_rate(ax=ax2, reg_params= east_0_north_reg_alt , df_x_y_yhat=east_0_north_df_alt , color=red, marker= "D", label= "Alps east 0",
                            bottom_labels=False)
    
    scatter_plot_laspe_rate(ax=ax2, reg_params= east_200_north_reg_alt , df_x_y_yhat=east_200_north_df_alt , color=green, marker= "^", label= "Alps east 200",
                            bottom_labels=False)
    
    ax2.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    
    #ax3 (south)
    
    scatter_plot_laspe_rate(ax=ax3, reg_params= control_south_reg_alt , df_x_y_yhat=control_south_df_alt , color=black, marker= "*", label= "Alps 100",
                            bottom_labels=True, xmax=3500, xmin=0, title= "[C] South",
                             ymax=0, ymin= -18,)
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_south_reg_alt , df_x_y_yhat=east_0_south_df_alt , color=red, marker= "^", label= "Alps east 0",
                           bottom_labels=True)
    scatter_plot_laspe_rate(ax=ax3, reg_params= east_200_south_reg_alt , df_x_y_yhat=east_200_south_df_alt , color=green, marker= "D", label= "Alps east 200",
                           bottom_labels=True)
    
    ax3.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5, ncol=1)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.88, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_store, "fig8.svg"), format= "svg", bbox_inches="tight", dpi=600)


plot_lape_rate_per_section()




