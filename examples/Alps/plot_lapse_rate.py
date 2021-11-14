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
exp_name_east_150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"
exp_name_Alps_150 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"
exp_name_Alps_0 = "t015_dkrz-mistral_e5w2.3_PI_Alps0_t159l31.6h"
exp_name_east_200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"

years= "1003_1017"
period = "1m"


# reading dataset
control_data, control_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period)
east_0_data, east_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_0, years=years,
                                                  period=period)
east_150_data, east_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_150, years=years,
                                                  period=period)
east_200_data, east_200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_200, years=years,
                                                  period=period)
Alps_150_data, Alps_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_150, years=years,
                                                  period=period)
Alps_0_data, Alps_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_0, years=years,
                                                  period=period)


#extracting variables and computing long-term means

#control
d18op = extract_var(Dataset=control_data , varname="d18op", units="per mil", Dataset_wiso= control_wiso)
elev = extract_var(Dataset=control_data , varname="elev", units="m")
d18op_slt = compute_lterm_mean(data=d18op, time="season", season_calendar="standard")
elev_slt = compute_lterm_mean(data=elev, time="season", season_calendar="standard")

# east_150
d18op_east_150 = extract_var(Dataset=east_150_data , varname="d18op", units="per mil", Dataset_wiso= east_150_wiso)
elev_east_150 = extract_var(Dataset=east_150_data , varname="elev", units="m")
east_150_d18op_slt = compute_lterm_mean(data=d18op_east_150, time="season", season_calendar="standard")
east_150_elev_slt = compute_lterm_mean(data=elev_east_150, time="season", season_calendar="standard")

# east_200
d18op_east_200 = extract_var(Dataset=east_200_data , varname="d18op", units="per mil", Dataset_wiso= east_200_wiso)
elev_east_200 = extract_var(Dataset=east_200_data , varname="elev", units="m")
east_200_d18op_slt = compute_lterm_mean(data=d18op_east_200, time="season", season_calendar="standard")
east_200_elev_slt = compute_lterm_mean(data=elev_east_200, time="season", season_calendar="standard")

#east_0
d18op_east_0 = extract_var(Dataset=east_0_data , varname="d18op", units="per mil", Dataset_wiso= east_0_wiso)
elev_east_0 = extract_var(Dataset=east_0_data , varname="elev", units="m")
east_0_d18op_slt = compute_lterm_mean(data=d18op_east_0, time="season", season_calendar="standard")
east_0_elev_slt = compute_lterm_mean(data=elev_east_0, time="season", season_calendar="standard")

#Alps_150
d18op_Alps_150 = extract_var(Dataset=Alps_150_data , varname="d18op", units="per mil", Dataset_wiso= Alps_150_wiso)
elev_Alps_150 = extract_var(Dataset=Alps_150_data , varname="elev", units="m")
Alps_150_d18op_slt = compute_lterm_mean(data=d18op_Alps_150, time="season", season_calendar="standard")
Alps_150_elev_slt = compute_lterm_mean(data=elev_Alps_150, time="season", season_calendar="standard")

#Alps_0
d18op_Alps_0 = extract_var(Dataset=Alps_0_data , varname="d18op", units="per mil", Dataset_wiso= Alps_0_wiso)
elev_Alps_0 = extract_var(Dataset=Alps_0_data , varname="elev", units="m")
Alps_0_d18op_slt = compute_lterm_mean(data=d18op_Alps_0, time="season", season_calendar="standard")
Alps_0_elev_slt = compute_lterm_mean(data=elev_Alps_0, time="season", season_calendar="standard")

# defining coordinates 
maxlat_east, minlat_east, maxlon_east, minlon_east = 48.5, 46, 18, 12
maxlat_west, minlat_west, maxlon_west, minlon_west = 46.5, 44, 8, 1
maxlat_south, minlat_south, maxlon_south, minlon_south = 47, 43, 13, 8
maxlat_north, minlat_north, maxlon_north, minlon_north = 50, 47, 14, 5

# extracting transects 
# control 
elev_control_east = extract_transect(data=elev_slt, maxlon=maxlon_east, minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east , sea_land_mask="yes", Dataset=control_data)
elev_control_north = extract_transect(data=elev_slt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask="yes", Dataset=control_data)
elev_control_west = extract_transect(data=elev_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=control_data)
elev_control_south = extract_transect(data=elev_slt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=control_data)

d18op_control_east = extract_transect(data=d18op_slt ,maxlon=maxlon_east, minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east, sea_land_mask="yes", Dataset=control_data)
d18op_control_north = extract_transect(data=d18op_slt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=control_data)
d18op_control_west = extract_transect(data=d18op_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=control_data)
d18op_control_south = extract_transect(data=d18op_slt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=control_data)


#east 150
elev_east_150_east = extract_transect(data=east_150_elev_slt, maxlon=maxlon_east, minlon=minlon_east, maxlat= maxlat_east , minlat=minlat_east, sea_land_mask="yes", Dataset=east_150_data)
elev_east_150_north = extract_transect(data=east_150_elev_slt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_150_data)
elev_east_150_west = extract_transect(data=east_150_elev_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_150_data)
elev_east_150_south = extract_transect(data=east_150_elev_slt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_150_data)


d18op_east_150_east = extract_transect(data=east_150_d18op_slt, maxlon=maxlon_east, minlon=minlon_east, maxlat= maxlat_east, minlat=minlat_east, sea_land_mask="yes", Dataset=east_150_data)
d18op_east_150_north = extract_transect(data=east_150_d18op_slt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_150_data)
d18op_east_150_west = extract_transect(data=east_150_d18op_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_150_data)
d18op_east_150_south = extract_transect(data=east_150_d18op_slt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_150_data)


#east 200
elev_east_200_east = extract_transect(data=east_200_elev_slt, maxlon=maxlon_east, minlon=minlon_east, maxlat= maxlat_east , minlat=minlat_east, sea_land_mask="yes", Dataset=east_200_data)
elev_east_200_north = extract_transect(data=east_200_elev_slt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_200_data)
elev_east_200_west = extract_transect(data=east_200_elev_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_200_data)
elev_east_200_south = extract_transect(data=east_200_elev_slt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_200_data)

d18op_east_200_east = extract_transect(data=east_200_d18op_slt, maxlon=maxlon_east, minlon=minlon_east, maxlat= maxlat_east, minlat=minlat_east, sea_land_mask="yes", Dataset=east_200_data)
d18op_east_200_north = extract_transect(data=east_200_d18op_slt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_200_data)
d18op_east_200_west = extract_transect(data=east_200_d18op_slt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_200_data)
d18op_east_200_south = extract_transect(data=east_200_d18op_slt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_200_data)


#east 0
elev_east_0_east = extract_transect(data=east_0_elev_slt , maxlon=maxlon_east , minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east , sea_land_mask="yes", Dataset=east_0_data)
elev_east_0_north = extract_transect(data=east_0_elev_slt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_0_data)
elev_east_0_west = extract_transect(data=east_0_elev_slt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_0_data)
elev_east_0_south = extract_transect(data=east_0_elev_slt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_0_data)

d18op_east_0_east = extract_transect(data=east_0_d18op_slt , maxlon=maxlon_east , minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east , sea_land_mask="yes", Dataset=east_0_data)
d18op_east_0_north = extract_transect(data=east_0_d18op_slt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_0_data)
d18op_east_0_west = extract_transect(data=east_0_d18op_slt ,maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_0_data)
d18op_east_0_south = extract_transect(data=east_0_d18op_slt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_0_data)

season = "JJA"

# performing regresion (change DJF to JJA for summer plots)
control_east_reg, control_east_df = linregression(data_x=elev_control_east , data_y=d18op_control_east, season=season, return_yhat="yes")
control_west_reg, control_west_df = linregression(data_x=elev_control_west , data_y=d18op_control_west, season=season, return_yhat="yes")
control_north_reg, control_north_df = linregression(data_x=elev_control_north , data_y=d18op_control_north, season=season, return_yhat="yes")
control_south_reg, control_south_df = linregression(data_x=elev_control_south , data_y=d18op_control_south, season=season, return_yhat="yes")

east_150_east_reg, east_150_east_df = linregression(data_x=elev_east_150_east , data_y=d18op_east_150_east, season=season, return_yhat="yes")
east_150_west_reg, east_150_west_df = linregression(data_x=elev_east_150_west , data_y=d18op_east_150_west, season=season, return_yhat="yes")
east_150_north_reg, east_150_north_df = linregression(data_x=elev_east_150_north , data_y=d18op_east_150_north, season=season, return_yhat="yes")
east_150_south_reg, east_150_south_df = linregression(data_x=elev_east_150_south , data_y=d18op_east_150_south, season=season, return_yhat="yes")

east_200_east_reg, east_200_east_df = linregression(data_x=elev_east_200_east , data_y=d18op_east_200_east, season=season, return_yhat="yes")
east_200_west_reg, east_200_west_df = linregression(data_x=elev_east_200_west , data_y=d18op_east_200_west, season=season, return_yhat="yes")
east_200_north_reg, east_200_north_df = linregression(data_x=elev_east_200_north , data_y=d18op_east_200_north, season=season, return_yhat="yes")
east_200_south_reg, east_200_south_df = linregression(data_x=elev_east_200_south , data_y=d18op_east_200_south, season=season, return_yhat="yes")

east_0_east_reg, east_0_east_df = linregression(data_x=elev_east_0_east , data_y=d18op_east_0_east, season=season, return_yhat="yes")
east_0_west_reg, east_0_west_df = linregression(data_x=elev_east_0_west , data_y=d18op_east_0_west, season=season, return_yhat="yes")
east_0_north_reg, east_0_north_df = linregression(data_x=elev_east_0_north , data_y=d18op_east_0_north, season=season, return_yhat="yes")
east_0_south_reg, east_0_south_df = linregression(data_x=elev_east_0_south , data_y=d18op_east_0_south, season=season, return_yhat="yes")

# plotting 

path_to_store = os.path.join(module_output_main_path, "plots")

fig, (ax1,ax2,ax3) = plt.subplots(nrows = 3, ncols = 1, figsize=(8, 17),)

apply_style_2()

scatter_plot_laspe_rate(ax=ax1, reg_params= control_east_reg , df_x_y_yhat=control_east_df , color=blue, marker= "*", label= "east",
                        ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [km]", title="[A] Alps 100 (" + season + ")", xmax=1500, xmin=0,
                        ymax=0, ymin= -9)
scatter_plot_laspe_rate(ax=ax1, reg_params= control_west_reg , df_x_y_yhat=control_west_df , color=red, marker= "o", label= "west",)
scatter_plot_laspe_rate(ax=ax1, reg_params= control_north_reg , df_x_y_yhat=control_north_df , color=black, marker= "D", label= "north",)
scatter_plot_laspe_rate(ax=ax1, reg_params= control_south_reg , df_x_y_yhat=control_south_df , color=green, marker= "^", label= "south",)

ax1.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5)
#ax1.set_aspect("equal", adjustable = "box")

# scatter_plot_laspe_rate(ax=ax2, reg_params= east_150_east_reg , df_x_y_yhat=east_150_east_df , color=blue, marker= "*", label= "east",
#                         ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [m]",title="[B] Alps east 150 (" + season + ")",
#                          xmax=3200, xmin=0, ymin=-14, ymax=0)
# scatter_plot_laspe_rate(ax=ax2, reg_params= east_150_west_reg , df_x_y_yhat=east_150_west_df , color=red, marker= "o", label= "west",)
# scatter_plot_laspe_rate(ax=ax2, reg_params= east_150_north_reg , df_x_y_yhat=east_150_north_df , color=black, marker= "D", label= "north",)
# scatter_plot_laspe_rate(ax=ax2, reg_params= east_150_south_reg , df_x_y_yhat=east_150_south_df , color=green, marker= "^", label= "south",)

scatter_plot_laspe_rate(ax=ax2, reg_params= east_200_east_reg , df_x_y_yhat=east_200_east_df , color=blue, marker= "*", label= "east",
                        ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [m]",title="[B] Alps east 200 (" + season + ")",
                         xmax=3200, xmin=0, ymin=-14, ymax=0)
scatter_plot_laspe_rate(ax=ax2, reg_params= east_200_west_reg , df_x_y_yhat=east_200_west_df , color=red, marker= "o", label= "west",)
scatter_plot_laspe_rate(ax=ax2, reg_params= east_200_north_reg , df_x_y_yhat=east_200_north_df , color=black, marker= "D", label= "north",)
scatter_plot_laspe_rate(ax=ax2, reg_params= east_200_south_reg , df_x_y_yhat=east_200_south_df , color=green, marker= "^", label= "south",)

ax2.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5)

scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_east_reg , df_x_y_yhat=east_0_east_df , color=blue, marker= "*", label= "east",
                        )
scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_west_reg , df_x_y_yhat=east_0_west_df , color=red, marker= "o", label= "west",
                        ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [m]", title="[C] Alps east 0 (" + season + ")",
                        xmax=1500, xmin=0, ymax=-2, ymin=-8)
scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_north_reg , df_x_y_yhat=east_0_north_df , color=black, marker= "D", label= "north",)
scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_south_reg , df_x_y_yhat=east_0_south_df , color=green, marker= "^", label= "south",)

ax3.legend(frameon=True, fontsize=15, loc="upper right",framealpha=0.5)

plt.tight_layout()
plt.subplots_adjust(left=0.15, right=0.88, top=0.97, bottom=0.05)
plt.savefig(os.path.join(path_to_store, "fig8_new.svg"), format= "svg", bbox_inches="tight", dpi=600)
