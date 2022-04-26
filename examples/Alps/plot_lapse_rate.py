#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 15:23:06 2021

@author: dboateng
This is example script for using Climat in analysing and visualising ECHAM module output (isotopic lapse rate)
The script contains directories of module outputs and path to save plot
Note: it is structured solely for the personal needs of the author, therefore, it must be adapted advisably.
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
exp_name_aw100e150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"

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
aw100e150_data, aw100e150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e150, years=years,
                                                  period=period)

aw200e100_data, aw200e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e100, years=years,
                                                  period=period)
aw200e0_data, aw200e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e0, years=years,
                                                  period=period)
aw200e200_data, aw200e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e200, years=years,
                                                  period=period)

#extracting variables and computing long-term means

#aw100e100
d18op_aw100e100 = extract_var(Dataset=aw100e100_data , varname="d18op", units="per mil", Dataset_wiso= aw100e100_wiso)
elev_aw100e100 = extract_var(Dataset=aw100e100_data , varname="elev", units="m")

d18op_aw100e100_slt = compute_lterm_mean(data=d18op_aw100e100, time="season", season_calendar="standard")
elev_aw100e100_slt = compute_lterm_mean(data=elev_aw100e100 , time="season", season_calendar="standard")

#aw100e0
d18op_aw100e0 = extract_var(Dataset=aw100e0_data , varname="d18op", units="per mil", Dataset_wiso= aw100e0_wiso)
elev_aw100e0 = extract_var(Dataset=aw100e0_data , varname="elev", units="m")

d18op_aw100e0_slt = compute_lterm_mean(data=d18op_aw100e0, time="season", season_calendar="standard")
elev_aw100e0_slt = compute_lterm_mean(data=elev_aw100e0 , time="season", season_calendar="standard")

#aw100e200
d18op_aw100e200 = extract_var(Dataset=aw100e200_data , varname="d18op", units="per mil", Dataset_wiso= aw100e200_wiso)
elev_aw100e200 = extract_var(Dataset=aw100e200_data , varname="elev", units="m")

d18op_aw100e200_slt = compute_lterm_mean(data=d18op_aw100e200, time="season", season_calendar="standard")
elev_aw100e200_slt = compute_lterm_mean(data=elev_aw100e200 , time="season", season_calendar="standard")

#aw100e150
d18op_aw100e150 = extract_var(Dataset=aw100e150_data , varname="d18op", units="per mil", Dataset_wiso= aw100e150_wiso)
elev_aw100e150 = extract_var(Dataset=aw100e150_data , varname="elev", units="m")

d18op_aw100e150_slt = compute_lterm_mean(data=d18op_aw100e150, time="season", season_calendar="standard")
elev_aw100e150_slt = compute_lterm_mean(data=elev_aw100e150 , time="season", season_calendar="standard")


#aw200e100
d18op_aw200e100 = extract_var(Dataset=aw200e100_data , varname="d18op", units="per mil", Dataset_wiso= aw200e100_wiso)
elev_aw200e100 = extract_var(Dataset=aw200e100_data , varname="elev", units="m")

d18op_aw200e100_slt = compute_lterm_mean(data=d18op_aw200e100, time="season", season_calendar="standard")
elev_aw200e100_slt = compute_lterm_mean(data=elev_aw200e100 , time="season", season_calendar="standard")

#aw200e0
d18op_aw200e0 = extract_var(Dataset=aw200e0_data , varname="d18op", units="per mil", Dataset_wiso= aw200e0_wiso)
elev_aw200e0 = extract_var(Dataset=aw200e0_data , varname="elev", units="m")

d18op_aw200e0_slt = compute_lterm_mean(data=d18op_aw200e0, time="season", season_calendar="standard")
elev_aw200e0_slt = compute_lterm_mean(data=elev_aw200e0 , time="season", season_calendar="standard")

#aw200e200
d18op_aw200e200 = extract_var(Dataset=aw200e200_data , varname="d18op", units="per mil", Dataset_wiso= aw200e200_wiso)
elev_aw200e200 = extract_var(Dataset=aw200e200_data , varname="elev", units="m")

d18op_aw200e200_slt = compute_lterm_mean(data=d18op_aw200e200, time="season", season_calendar="standard")
elev_aw200e200_slt = compute_lterm_mean(data=elev_aw200e200 , time="season", season_calendar="standard")

# defining coordinates 

maxlat_west, minlat_west, maxlon_west, minlon_west = 47, 44, 8, 1
maxlat_south, minlat_south, maxlon_south, minlon_south = 47, 43, 15, 7.5
maxlat_north, minlat_north, maxlon_north, minlon_north = 50, 46.5, 16, 5

# extracting transects 


# aw100e100 

elev_aw100e100_north_slt = extract_transect(data=elev_aw100e100_slt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, Dataset=aw100e100_data)
elev_aw100e100_west_slt = extract_transect(data=elev_aw100e100_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw100e100_data)
elev_aw100e100_south_slt = extract_transect(data=elev_aw100e100_slt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw100e100_data)


d18op_aw100e100_north_slt = extract_transect(data=d18op_aw100e100_slt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=aw100e100_data)
d18op_aw100e100_west_slt = extract_transect(data=d18op_aw100e100_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw100e100_data)
d18op_aw100e100_south_slt = extract_transect(data=d18op_aw100e100_slt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw100e100_data)


#aw100e0
elev_aw100e0_north_slt = extract_transect(data=elev_aw100e0_slt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, Dataset=aw100e0_data)
elev_aw100e0_west_slt = extract_transect(data=elev_aw100e0_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw100e0_data)
elev_aw100e0_south_slt = extract_transect(data=elev_aw100e0_slt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw100e0_data)


d18op_aw100e0_north_slt = extract_transect(data=d18op_aw100e0_slt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=aw100e0_data)
d18op_aw100e0_west_slt = extract_transect(data=d18op_aw100e0_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw100e0_data)
d18op_aw100e0_south_slt = extract_transect(data=d18op_aw100e0_slt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw100e0_data)


#aw100e200
elev_aw100e200_north_slt = extract_transect(data=elev_aw100e200_slt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, Dataset=aw100e200_data)
elev_aw100e200_west_slt = extract_transect(data=elev_aw100e200_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw100e200_data)
elev_aw100e200_south_slt = extract_transect(data=elev_aw100e200_slt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw100e200_data)


d18op_aw100e200_north_slt = extract_transect(data=d18op_aw100e200_slt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=aw100e200_data)
d18op_aw100e200_west_slt = extract_transect(data=d18op_aw100e200_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw100e200_data)
d18op_aw100e200_south_slt = extract_transect(data=d18op_aw100e200_slt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw100e200_data)


#aw200e100
elev_aw200e100_north_slt = extract_transect(data=elev_aw200e100_slt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, Dataset=aw200e100_data)
elev_aw200e100_west_slt = extract_transect(data=elev_aw200e100_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw200e100_data)
elev_aw200e100_south_slt = extract_transect(data=elev_aw200e100_slt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw200e100_data)


d18op_aw200e100_north_slt = extract_transect(data=d18op_aw200e100_slt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=aw200e100_data)
d18op_aw200e100_west_slt = extract_transect(data=d18op_aw200e100_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw200e100_data)
d18op_aw200e100_south_slt = extract_transect(data=d18op_aw200e100_slt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw200e100_data)

#aw200e0
elev_aw200e0_north_slt = extract_transect(data=elev_aw200e0_slt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, Dataset=aw200e0_data)
elev_aw200e0_west_slt = extract_transect(data=elev_aw200e0_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw200e0_data)
elev_aw200e0_south_slt = extract_transect(data=elev_aw200e0_slt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw200e0_data)


d18op_aw200e0_north_slt = extract_transect(data=d18op_aw200e0_slt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=aw200e0_data)
d18op_aw200e0_west_slt = extract_transect(data=d18op_aw200e0_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw200e0_data)
d18op_aw200e0_south_slt = extract_transect(data=d18op_aw200e0_slt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw200e0_data)

#aw200e200
elev_aw200e200_north_slt = extract_transect(data=elev_aw200e200_slt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask=True, Dataset=aw200e200_data)
elev_aw200e200_west_slt = extract_transect(data=elev_aw200e200_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw200e200_data)
elev_aw200e200_south_slt = extract_transect(data=elev_aw200e200_slt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw200e200_data)


d18op_aw200e200_north_slt = extract_transect(data=d18op_aw200e200_slt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask=True, Dataset=aw200e200_data)
d18op_aw200e200_west_slt = extract_transect(data=d18op_aw200e200_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask=True, Dataset=aw200e200_data)
d18op_aw200e200_south_slt = extract_transect(data=d18op_aw200e200_slt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask=True, Dataset=aw200e200_data)


season = "JJA"

# performing regression (change DJF to JJA for summer plots)

#AW100E100
aw100e100_west_reg_slt, aw100e100_west_df_slt = linregression(data_x=elev_aw100e100_west_slt , data_y=d18op_aw100e100_west_slt, season=season, return_yhat=True)
aw100e100_north_reg_slt, aw100e100_north_df_slt = linregression(data_x=elev_aw100e100_north_slt , data_y=d18op_aw100e100_north_slt, season=season, return_yhat=True)
aw100e100_south_reg_slt, aw100e100_south_df_slt = linregression(data_x=elev_aw100e100_south_slt , data_y=d18op_aw100e100_south_slt, season=season, return_yhat=True)


#AW100E200
aw100e200_west_reg_slt, aw100e200_west_df_slt = linregression(data_x=elev_aw100e200_west_slt , data_y=d18op_aw100e200_west_slt, season=season, return_yhat=True)
aw100e200_north_reg_slt, aw100e200_north_df_slt = linregression(data_x=elev_aw100e200_north_slt , data_y=d18op_aw100e200_north_slt, season=season, return_yhat=True)
aw100e200_south_reg_slt, aw100e200_south_df_slt = linregression(data_x=elev_aw100e200_south_slt , data_y=d18op_aw100e200_south_slt, season=season, return_yhat=True)


#AW100E0
aw100e0_west_reg_slt, aw100e0_west_df_slt = linregression(data_x=elev_aw100e0_west_slt , data_y=d18op_aw100e0_west_slt, season=season, return_yhat=True)
aw100e0_north_reg_slt, aw100e0_north_df_slt = linregression(data_x=elev_aw100e0_north_slt , data_y=d18op_aw100e0_north_slt, season=season, return_yhat=True)
aw100e0_south_reg_slt, aw100e0_south_df_slt = linregression(data_x=elev_aw100e0_south_slt , data_y=d18op_aw100e0_south_slt, season=season, return_yhat=True)

#AW200E100
aw200e100_west_reg_slt, aw200e100_west_df_slt = linregression(data_x=elev_aw200e100_west_slt , data_y=d18op_aw200e100_west_slt, season=season, return_yhat=True)
aw200e100_north_reg_slt, aw200e100_north_df_slt = linregression(data_x=elev_aw200e100_north_slt , data_y=d18op_aw200e100_north_slt, season=season, return_yhat=True)
aw200e100_south_reg_slt, aw200e100_south_df_slt = linregression(data_x=elev_aw200e100_south_slt , data_y=d18op_aw200e100_south_slt, season=season, return_yhat=True)


#AW200E200
aw200e200_west_reg_slt, aw200e200_west_df_slt = linregression(data_x=elev_aw200e200_west_slt , data_y=d18op_aw200e200_west_slt, season=season, return_yhat=True)
aw200e200_north_reg_slt, aw200e200_north_df_slt = linregression(data_x=elev_aw200e200_north_slt , data_y=d18op_aw200e200_north_slt, season=season, return_yhat=True)
aw200e200_south_reg_slt, aw200e200_south_df_slt = linregression(data_x=elev_aw200e200_south_slt , data_y=d18op_aw200e200_south_slt, season=season, return_yhat=True)


#AW200E0
aw200e0_west_reg_slt, aw200e0_west_df_slt = linregression(data_x=elev_aw200e0_west_slt , data_y=d18op_aw200e0_west_slt, season=season, return_yhat=True)
aw200e0_north_reg_slt, aw200e0_north_df_slt = linregression(data_x=elev_aw200e0_north_slt , data_y=d18op_aw200e0_north_slt, season=season, return_yhat=True)
aw200e0_south_reg_slt, aw200e0_south_df_slt = linregression(data_x=elev_aw200e0_south_slt , data_y=d18op_aw200e0_south_slt, season=season, return_yhat=True)



# plotting 

path_to_store = os.path.join(module_output_main_path, "plots")


apply_style(fontsize=22, style=None, linewidth=2)

def plot_lape_rate_per_section():

    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(25, 10), )
    
    #ax1 (west)
    
    scatter_plot_laspe_rate(ax=ax1, reg_params= aw100e100_west_reg_slt , df_x_y_yhat=aw100e100_west_df_slt , color=black, marker= "*", label= "CTL",
                           title="[A] West", xmax=3500, xmin=0,
                            ymax=0, ymin= -16, bottom_labels=True)
    scatter_plot_laspe_rate(ax=ax1, reg_params= aw100e0_west_reg_slt , df_x_y_yhat=aw100e0_west_df_slt , color=red, marker= "D", label= "W1E0",
                           bottom_labels=True)
    scatter_plot_laspe_rate(ax=ax1, reg_params= aw100e200_west_reg_slt , df_x_y_yhat=aw100e200_west_df_slt , color=green, marker= "^", label= "W1E2",
                           bottom_labels=True)
    
    
    
    ax1.legend(frameon=True, fontsize=20, loc="upper right", framealpha=0.5, ncol=1)
    ax1.grid(visible=False)
    
    #ax2 (north)
    scatter_plot_laspe_rate(ax=ax2, reg_params= aw100e100_north_reg_slt , df_x_y_yhat=aw100e100_north_df_slt , color=black, marker= "*", label= "CTL",
                            left_labels=False, xmax=3500, xmin=0, title= "[B] North",
                             ymax=0, ymin= -16,)
    scatter_plot_laspe_rate(ax=ax2, reg_params= aw100e0_north_reg_slt , df_x_y_yhat=aw100e0_north_df_slt , color=red, marker= "D", label= "W1E0",
                            left_labels=False)
    
    scatter_plot_laspe_rate(ax=ax2, reg_params= aw100e200_north_reg_slt , df_x_y_yhat=aw100e200_north_df_slt , color=green, marker= "^", label= "W1E2",
                            left_labels=False)
    
    
    ax2.legend(frameon=True, fontsize=20, loc="upper right", framealpha=0.5, ncol=1)
    ax2.grid(visible=False)
    
    #ax3 (south)
    scatter_plot_laspe_rate(ax=ax3, reg_params= aw100e100_south_reg_slt , df_x_y_yhat=aw100e100_south_df_slt , color=black, marker= "*", label= "CTL",
                            left_labels=False, xmax=3500, xmin=0, title= "[C] South",
                             ymax=0, ymin= -16,)
    scatter_plot_laspe_rate(ax=ax3, reg_params= aw100e0_south_reg_slt , df_x_y_yhat=aw100e0_south_df_slt , color=red, marker= "^", label= "W1E0",
                           left_labels=False)
    scatter_plot_laspe_rate(ax=ax3, reg_params= aw100e200_south_reg_slt , df_x_y_yhat=aw100e200_south_df_slt , color=green, marker= "D", label= "W1E2",
                           left_labels=False)
    
    
    ax3.legend(frameon=True, fontsize=20, loc="upper right", framealpha=0.5, ncol=1)
    ax3.grid(visible=False)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.88, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_store, "fig6.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "fig6.png"), format= "png", bbox_inches="tight", dpi=300)
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, figsize=(25, 10), )
    
    #ax1 (west)
    scatter_plot_laspe_rate(ax=ax1, reg_params= aw200e200_west_reg_slt , df_x_y_yhat=aw200e200_west_df_slt , color=olive, marker= "h", label= "W2E2",
                           title="[A] West", xmax=3500, xmin=0,
                            ymax=0, ymin= -16, )
    
    scatter_plot_laspe_rate(ax=ax1, reg_params= aw200e0_west_reg_slt , df_x_y_yhat=aw200e0_west_df_slt , color=purple, marker= "p", label= "W2E0",
                           )
    scatter_plot_laspe_rate(ax=ax1, reg_params= aw200e100_west_reg_slt , df_x_y_yhat=aw200e100_west_df_slt , color=golden, marker= "s", label= "W2E1",
                           )
    scatter_plot_laspe_rate(ax=ax1, reg_params= aw100e100_west_reg_slt , df_x_y_yhat=aw100e100_west_df_slt , color=black, marker= "*", label= "CTL",
                           )
    
    ax1.legend(frameon=True, fontsize=20, loc="upper right", framealpha=0.5, ncol=1)
    ax1.grid(visible=False)
    
    
    #ax2 (north)
    scatter_plot_laspe_rate(ax=ax2, reg_params= aw200e200_north_reg_slt , df_x_y_yhat=aw200e200_north_df_slt , color=olive, marker= "h", label= "W2E2",
                            left_labels=False, xmax=3500, xmin=0, title= "[B] North",
                             ymax=0, ymin= -16,)
    scatter_plot_laspe_rate(ax=ax2, reg_params= aw200e0_north_reg_slt , df_x_y_yhat=aw200e0_north_df_slt , color=purple, marker= "p", label= "W2E0",
                            left_labels=False)
    scatter_plot_laspe_rate(ax=ax2, reg_params= aw200e100_north_reg_slt , df_x_y_yhat=aw200e100_north_df_slt , color=golden, marker= "s", label= "W2E1",
                            left_labels=False)
    scatter_plot_laspe_rate(ax=ax2, reg_params= aw100e100_north_reg_slt , df_x_y_yhat=aw100e100_north_df_slt , color=black, marker= "*", label= "CTL",
                            left_labels=False,)
    
    ax2.legend(frameon=True, fontsize=20, loc="upper right", framealpha=0.5, ncol=1)
    ax2.grid(visible=False)
    
    
    #ax3 (south)
    scatter_plot_laspe_rate(ax=ax3, reg_params= aw200e200_south_reg_slt , df_x_y_yhat=aw200e200_south_df_slt , color=olive, marker= "h", label= "W2E2",
                            left_labels=False, xmax=3500, xmin=0, title= "[C] South",
                             ymax=0, ymin= -16,)
    
    scatter_plot_laspe_rate(ax=ax3, reg_params= aw200e0_south_reg_slt , df_x_y_yhat=aw200e0_south_df_slt , color=purple, marker= "p", label= "W2E0",
                           left_labels=False)
    scatter_plot_laspe_rate(ax=ax3, reg_params= aw200e100_south_reg_slt , df_x_y_yhat=aw200e100_south_df_slt , color=golden, marker= "s", label= "W2E1",
                           left_labels=False)
    scatter_plot_laspe_rate(ax=ax3, reg_params= aw100e100_south_reg_slt , df_x_y_yhat=aw100e100_south_df_slt , color=black, marker= "*", label= "CTL",
                            left_labels=False,)
    
    ax3.legend(frameon=True, fontsize=20, loc="upper right", framealpha=0.5, ncol=1)
    ax3.grid(visible=False)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.88, top=0.97, bottom=0.05)
    plt.savefig(os.path.join(path_to_store, "fig7.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "fig7.png"), format= "png", bbox_inches="tight", dpi=300)



if __name__ == '__main__':
    plot_lape_rate_per_section()




