#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:53:03 2021

@author: dboateng
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"

exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east_150 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_Alps_150 = "a005_hpc-bw_e5w2.3_t159_PI_Alps_150_t159l31.6h"
exp_name_Alps_0 = "a006_hpc-bw_e5w2.3_t159_PI_Alps_0_t159l31.6h"
exp_name_east_50 = "a004_hpc-bw_e5w2.3_t159_PI_Alps_east_50_t159l31.6h"

years= "1003_1017"
period = "1m"


# reading dataset
control_data, control_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period)
east_0_data, east_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_0, years=years,
                                                  period=period)
east_150_data, east_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_150, years=years,
                                                  period=period)
Alps_150_data, Alps_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_150, years=years,
                                                  period=period)
Alps_0_data, Alps_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_0, years=years,
                                                  period=period)
east_50_data, east_50_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_50, years=years,
                                                  period=period)

#extracting variables and computing long-term means

# control
d18op = extract_var(Dataset=control_data , varname="d18op", units="per mil", Dataset_wiso= control_wiso)
elev = extract_var(Dataset=control_data , varname="elev", units="m")
d18op_alt = compute_lterm_mean(data=d18op, time="annual")
elev_alt = compute_lterm_mean(data=elev, time="annual")

d18op_east_150 = extract_var(Dataset=east_150_data , varname="d18op", units="per mil", Dataset_wiso= east_150_wiso)
elev_east_150 = extract_var(Dataset=east_150_data , varname="elev", units="m")
east_150_d18op_alt = compute_lterm_mean(data=d18op_east_150, time="annual")
east_150_elev_alt = compute_lterm_mean(data=elev_east_150, time="annual")

d18op_east_0 = extract_var(Dataset=east_0_data , varname="d18op", units="per mil", Dataset_wiso= east_0_wiso)
elev_east_0 = extract_var(Dataset=east_0_data , varname="elev", units="m")
east_0_d18op_alt = compute_lterm_mean(data=d18op_east_0, time="annual")
east_0_elev_alt = compute_lterm_mean(data=elev_east_0, time="annual")

d18op_east_50 = extract_var(Dataset=east_50_data , varname="d18op", units="per mil", Dataset_wiso= east_50_wiso)
elev_east_50 = extract_var(Dataset=east_50_data , varname="elev", units="m")
east_50_d18op_alt = compute_lterm_mean(data=d18op_east_50, time="annual")
east_50_elev_alt = compute_lterm_mean(data=elev_east_50, time="annual")

d18op_Alps_150 = extract_var(Dataset=Alps_150_data , varname="d18op", units="per mil", Dataset_wiso= Alps_150_wiso)
elev_Alps_150 = extract_var(Dataset=Alps_150_data , varname="elev", units="m")
Alps_150_d18op_alt = compute_lterm_mean(data=d18op_Alps_150, time="annual")
Alps_150_elev_alt = compute_lterm_mean(data=elev_Alps_150, time="annual")

d18op_Alps_0 = extract_var(Dataset=Alps_0_data , varname="d18op", units="per mil", Dataset_wiso= Alps_0_wiso)
elev_Alps_0 = extract_var(Dataset=Alps_0_data , varname="elev", units="m")
Alps_0_d18op_alt = compute_lterm_mean(data=d18op_Alps_0, time="annual")
Alps_0_elev_alt = compute_lterm_mean(data=elev_Alps_0, time="annual")


# defining coordinates 
maxlat_east, minlat_east, maxlon_east, minlon_east = 48.5, 46, 18, 12
maxlat_west, minlat_west, maxlon_west, minlon_west = 46.5, 44, 8, 1
maxlat_south, minlat_south, maxlon_south, minlon_south = 47, 43, 13, 8
maxlat_north, minlat_north, maxlon_north, minlon_north = 50, 47, 14, 5


# extracting transects 
# control 
elev_control_east = extract_transect(data=elev_alt, maxlon=maxlon_east, minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east , sea_land_mask="yes", Dataset=control_data)
elev_control_north = extract_transect(data=elev_alt, maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north , sea_land_mask="yes", Dataset=control_data)
elev_control_west = extract_transect(data=elev_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=control_data)
elev_control_south = extract_transect(data=elev_alt, maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=control_data)

d18op_control_east = extract_transect(data=d18op_alt ,maxlon=maxlon_east, minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east, sea_land_mask="yes", Dataset=control_data)
d18op_control_north = extract_transect(data=d18op_alt , maxlon=maxlon_north, minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=control_data)
d18op_control_west = extract_transect(data=d18op_alt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=control_data)
d18op_control_south = extract_transect(data=d18op_alt ,maxlon=maxlon_south, minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=control_data)


#east 150
elev_east_150_east = extract_transect(data=east_150_elev_alt, maxlon=maxlon_east, minlon=minlon_east, maxlat= maxlat_east , minlat=minlat_east, sea_land_mask="yes", Dataset=east_150_data)
elev_east_150_north = extract_transect(data=east_150_elev_alt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_150_data)
elev_east_150_west = extract_transect(data=east_150_elev_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_150_data)
elev_east_150_south = extract_transect(data=east_150_elev_alt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_150_data)

d18op_east_150_east = extract_transect(data=east_150_d18op_alt, maxlon=maxlon_east, minlon=minlon_east, maxlat= maxlat_east, minlat=minlat_east, sea_land_mask="yes", Dataset=east_150_data)
d18op_east_150_north = extract_transect(data=east_150_d18op_alt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_150_data)
d18op_east_150_west = extract_transect(data=east_150_d18op_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_150_data)
d18op_east_150_south = extract_transect(data=east_150_d18op_alt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_150_data)

#east 0
elev_east_0_east = extract_transect(data=east_0_elev_alt , maxlon=maxlon_east , minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east , sea_land_mask="yes", Dataset=east_0_data)
elev_east_0_north = extract_transect(data=east_0_elev_alt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_0_data)
elev_east_0_west = extract_transect(data=east_0_elev_alt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_0_data)
elev_east_0_south = extract_transect(data=east_0_elev_alt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_0_data)

d18op_east_0_east = extract_transect(data=east_0_d18op_alt , maxlon=maxlon_east , minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east , sea_land_mask="yes", Dataset=east_0_data)
d18op_east_0_north = extract_transect(data=east_0_d18op_alt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=east_0_data)
d18op_east_0_west = extract_transect(data=east_0_d18op_alt ,maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=east_0_data)
d18op_east_0_south = extract_transect(data=east_0_d18op_alt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=east_0_data)


#Alps 150
elev_Alps_150_east = extract_transect(data=Alps_150_elev_alt, maxlon=maxlon_east, minlon=minlon_east, maxlat= maxlat_east , minlat=minlat_east, sea_land_mask="yes", Dataset=Alps_150_data)
elev_Alps_150_north = extract_transect(data=Alps_150_elev_alt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=Alps_150_data)
elev_Alps_150_west = extract_transect(data=Alps_150_elev_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=Alps_150_data)
elev_Alps_150_south = extract_transect(data=Alps_150_elev_alt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=Alps_150_data)

d18op_Alps_150_east = extract_transect(data=Alps_150_d18op_alt, maxlon=maxlon_east, minlon=minlon_east, maxlat= maxlat_east, minlat=minlat_east, sea_land_mask="yes", Dataset=Alps_150_data)
d18op_Alps_150_north = extract_transect(data=Alps_150_d18op_alt,  maxlon=maxlon_north, minlon=minlon_north, maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=Alps_150_data)
d18op_Alps_150_west = extract_transect(data=Alps_150_d18op_alt, maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=Alps_150_data)
d18op_Alps_150_south = extract_transect(data=Alps_150_d18op_alt, maxlon=maxlon_south, minlon=minlon_south, maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=Alps_150_data)

#Alps 0
elev_Alps_0_east = extract_transect(data=Alps_0_elev_alt , maxlon=maxlon_east , minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east , sea_land_mask="yes", Dataset=Alps_0_data)
elev_Alps_0_north = extract_transect(data=Alps_0_elev_alt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=Alps_0_data)
elev_Alps_0_west = extract_transect(data=Alps_0_elev_alt , maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=Alps_0_data)
elev_Alps_0_south = extract_transect(data=Alps_0_elev_alt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=Alps_0_data)

d18op_Alps_0_east = extract_transect(data=Alps_0_d18op_alt , maxlon=maxlon_east , minlon=minlon_east , maxlat= maxlat_east , minlat=minlat_east , sea_land_mask="yes", Dataset=Alps_0_data)
d18op_Alps_0_north = extract_transect(data=Alps_0_d18op_alt ,  maxlon=maxlon_north , minlon=minlon_north , maxlat=maxlat_north , minlat=minlat_north, sea_land_mask="yes", Dataset=Alps_0_data)
d18op_Alps_0_west = extract_transect(data=Alps_0_d18op_alt ,maxlon=maxlon_west,  minlon=minlon_west, maxlat=maxlat_west, minlat=minlat_west, sea_land_mask="yes", Dataset=Alps_0_data)
d18op_Alps_0_south = extract_transect(data=Alps_0_d18op_alt , maxlon=maxlon_south , minlon=minlon_south , maxlat=maxlat_south, minlat=minlat_south, sea_land_mask="yes", Dataset=Alps_0_data)


control_east_reg, control_east_df = linregression(data_x=elev_control_east , data_y=d18op_control_east, return_yhat="yes")
control_west_reg, control_west_df = linregression(data_x=elev_control_west , data_y=d18op_control_west, return_yhat="yes")
control_north_reg, control_north_df = linregression(data_x=elev_control_north , data_y=d18op_control_north, return_yhat="yes")
control_south_reg, control_south_df = linregression(data_x=elev_control_south , data_y=d18op_control_south, return_yhat="yes")

east_150_east_reg, east_150_east_df = linregression(data_x=elev_east_150_east , data_y=d18op_east_150_east,  return_yhat="yes")
east_150_west_reg, east_150_west_df = linregression(data_x=elev_east_150_west , data_y=d18op_east_150_west,  return_yhat="yes")
east_150_north_reg, east_150_north_df = linregression(data_x=elev_east_150_north , data_y=d18op_east_150_north,  return_yhat="yes")
east_150_south_reg, east_150_south_df = linregression(data_x=elev_east_150_south , data_y=d18op_east_150_south, return_yhat="yes")

east_0_east_reg, east_0_east_df = linregression(data_x=elev_east_0_east , data_y=d18op_east_0_east, return_yhat="yes")
east_0_west_reg, east_0_west_df = linregression(data_x=elev_east_0_west , data_y=d18op_east_0_west, return_yhat="yes")
east_0_north_reg, east_0_north_df = linregression(data_x=elev_east_0_north , data_y=d18op_east_0_north, return_yhat="yes")
east_0_south_reg, east_0_south_df = linregression(data_x=elev_east_0_south , data_y=d18op_east_0_south, return_yhat="yes")

# Alps 150

Alps_150_east_reg, Alps_150_east_df = linregression(data_x=elev_Alps_150_east , data_y=d18op_Alps_150_east, return_yhat="yes")
Alps_150_west_reg, Alps_150_west_df = linregression(data_x=elev_Alps_150_west , data_y=d18op_Alps_150_west, return_yhat="yes")
Alps_150_north_reg, Alps_150_north_df = linregression(data_x=elev_Alps_150_north , data_y=d18op_Alps_150_north, return_yhat="yes")
Alps_150_south_reg, Alps_150_south_df = linregression(data_x=elev_Alps_150_south , data_y=d18op_Alps_150_south, return_yhat="yes")

#Alps 0
Alps_0_east_reg, Alps_0_east_df = linregression(data_x=elev_Alps_0_east , data_y=d18op_Alps_0_east, return_yhat="yes")
Alps_0_west_reg, Alps_0_west_df = linregression(data_x=elev_Alps_0_west , data_y=d18op_Alps_0_west, return_yhat="yes")
Alps_0_north_reg, Alps_0_north_df = linregression(data_x=elev_Alps_0_north , data_y=d18op_Alps_0_north, return_yhat="yes")
Alps_0_south_reg, Alps_0_south_df = linregression(data_x=elev_Alps_0_south , data_y=d18op_Alps_0_south, return_yhat="yes")


# plotting 

path_to_store = os.path.join(module_output_main_path, "plots")

fig, (ax1,ax2,ax3) = plt.subplots(nrows = 3, ncols = 1, figsize=(8, 15),)

apply_style_2()

scatter_plot_laspe_rate(ax=ax1, reg_params= control_east_reg , df_x_y_yhat=control_east_df , color=blue, marker= "*", label= "east",
                        ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [km]", title="[A] Alps 100", xmax=1500, xmin=0,
                        ymax=-2, ymin= -12)
scatter_plot_laspe_rate(ax=ax1, reg_params= control_west_reg , df_x_y_yhat=control_west_df , color=red, marker= "o", label= "west",)
scatter_plot_laspe_rate(ax=ax1, reg_params= control_north_reg , df_x_y_yhat=control_north_df , color=black, marker= "D", label= "north",)
scatter_plot_laspe_rate(ax=ax1, reg_params= control_south_reg , df_x_y_yhat=control_south_df , color=green, marker= "^", label= "south",)

ax1.legend(frameon=True, fontsize=15, loc="upper right",framealpha=0.5)
#ax1.set_aspect("equal", adjustable = "box")

scatter_plot_laspe_rate(ax=ax2, reg_params= east_150_east_reg , df_x_y_yhat=east_150_east_df , color=blue, marker= "*", label= "east",
                        ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [m]",title="[B] Alps east 200",
                        xmax=3500, xmin=0, ymin=-18, ymax=-2)
scatter_plot_laspe_rate(ax=ax2, reg_params= east_150_west_reg , df_x_y_yhat=east_150_west_df , color=red, marker= "o", label= "west",)
scatter_plot_laspe_rate(ax=ax2, reg_params= east_150_north_reg , df_x_y_yhat=east_150_north_df , color=black, marker= "D", label= "north",)
scatter_plot_laspe_rate(ax=ax2, reg_params= east_150_south_reg , df_x_y_yhat=east_150_south_df , color=green, marker= "^", label= "south",)

ax2.legend(frameon=True, fontsize=15, loc="upper right", framealpha=0.5)

scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_east_reg , df_x_y_yhat=east_0_east_df , color=blue, marker= "*", label= "east",
                        ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [m]", title="[C] Alps east 0",
                        xmax=1500, xmin=0, ymax=-2, ymin=-12)
scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_west_reg , df_x_y_yhat=east_0_west_df , color=red, marker= "o", label= "west",)
scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_north_reg , df_x_y_yhat=east_0_north_df , color=black, marker= "D", label= "north",)
scatter_plot_laspe_rate(ax=ax3, reg_params= east_0_south_reg , df_x_y_yhat=east_0_south_df , color=green, marker= "^", label= "south",)

ax3.legend(frameon=True, fontsize=15, loc="upper right",framealpha=0.5)

plt.tight_layout()
plt.subplots_adjust(left=0.15, right=0.95, top=0.96, bottom=0.05)
plt.savefig(os.path.join(path_to_store, "figS12_new.svg"), format= "svg", bbox_inches="tight", dpi=600)


# # Alps 150, Alps 0
# fig, (ax1,ax2,ax3) = plt.subplots(nrows = 3, ncols = 1, figsize=(8, 15),)

# apply_style_2()

# scatter_plot_laspe_rate(ax=ax1, reg_params= control_east_reg , df_x_y_yhat=control_east_df , color=blue, marker= "*", label= "east",
#                         ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [km]", title="[A] Alps 100")
# scatter_plot_laspe_rate(ax=ax1, reg_params= control_west_reg , df_x_y_yhat=control_west_df , color=red, marker= "o", label= "west",)
# scatter_plot_laspe_rate(ax=ax1, reg_params= control_north_reg , df_x_y_yhat=control_north_df , color=black, marker= "D", label= "north",)
# scatter_plot_laspe_rate(ax=ax1, reg_params= control_south_reg , df_x_y_yhat=control_south_df , color=green, marker= "^", label= "south",)

# ax1.legend(frameon=True, fontsize=15, loc="upper right",)
# #ax1.set_aspect("equal", adjustable = "box")

# scatter_plot_laspe_rate(ax=ax2, reg_params= Alps_150_east_reg , df_x_y_yhat=Alps_150_east_df , color=blue, marker= "*", label= "east",
#                         ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [m]",title="[B] Alps 150")
# scatter_plot_laspe_rate(ax=ax2, reg_params= Alps_150_west_reg , df_x_y_yhat=Alps_150_west_df , color=red, marker= "o", label= "west",)
# scatter_plot_laspe_rate(ax=ax2, reg_params= Alps_150_north_reg , df_x_y_yhat=Alps_150_north_df , color=black, marker= "D", label= "north",)
# scatter_plot_laspe_rate(ax=ax2, reg_params= Alps_150_south_reg , df_x_y_yhat=Alps_150_south_df , color=green, marker= "^", label= "south",)

# ax2.legend(frameon=True, fontsize=15, loc="upper right",)

# scatter_plot_laspe_rate(ax=ax3, reg_params= Alps_0_east_reg , df_x_y_yhat=Alps_0_east_df , color=blue, marker= "*", label= "east",
#                         ylabel="$\delta^{18}$Op ‰ vs SMOW", xlabel= "Elevation [m]", title="[C] Alps 0")
# scatter_plot_laspe_rate(ax=ax3, reg_params= Alps_0_west_reg , df_x_y_yhat=Alps_0_west_df , color=red, marker= "o", label= "west",)
# scatter_plot_laspe_rate(ax=ax3, reg_params= Alps_0_north_reg , df_x_y_yhat=Alps_0_north_df , color=black, marker= "D", label= "north",)
# scatter_plot_laspe_rate(ax=ax3, reg_params= Alps_0_south_reg , df_x_y_yhat=Alps_0_south_df , color=green, marker= "^", label= "south",)

# ax3.legend(frameon=True, fontsize=15, loc="upper right",)

# plt.tight_layout()
# plt.subplots_adjust(left=0.06, right=0.88, top=0.98, bottom=0.05)
# plt.savefig(os.path.join(path_to_store, "figS13.png"), format= "png", bbox_inches="tight", dpi=600)
