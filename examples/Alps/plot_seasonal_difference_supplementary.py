#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 18:44:44 2021

@author: dboateng
This script visualise the long-term seasonal means for some set of experiments on the Alps
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)


import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_0 = "a006_hpc-bw_e5w2.3_t159_PI_Alps_0_t159l31.6h"
exp_name_150 = "a005_hpc-bw_e5w2.3_t159_PI_Alps_150_t159l31.6h"
years= "1003_1017"
period = "1m"


# reading dataset
Alps_100_data, Alps_100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period)
Alps_0_data, Alps_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_0, years=years,
                                                  period=period)
Alps_150_data, Alps_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_150, years=years,
                                                  period=period)

#  extracting variables 
#control
temp2_Alps_100 = extract_var(Dataset=Alps_100_data , varname="temp2", units="°C")
prec_Alps_100 = extract_var(Dataset= Alps_100_data , varname="prec", units="mm/month")
d18op_Alps_100 = extract_var(Dataset=Alps_100_data , varname="d18op", units="per mil", Dataset_wiso= Alps_100_wiso)
u10_Alps_100 = extract_var(Dataset=Alps_100_data , varname="u10")
v10_Alps_100 = extract_var(Dataset=Alps_100_data , varname="v10")

#east_150
temp2_Alps_150 = extract_var(Dataset=Alps_150_data , varname="temp2", units="°C")
prec_Alps_150 = extract_var(Dataset= Alps_150_data , varname="prec", units="mm/month")
d18op_Alps_150 = extract_var(Dataset=Alps_150_data , varname="d18op", units="per mil", Dataset_wiso= Alps_150_wiso)
u10_Alps_150 = extract_var(Dataset=Alps_150_data , varname="u10")
v10_Alps_150 = extract_var(Dataset=Alps_150_data , varname="v10")

#east_0
temp2_Alps_0 = extract_var(Dataset=Alps_0_data , varname="temp2", units="°C")
prec_Alps_0 = extract_var(Dataset= Alps_0_data , varname="prec", units="mm/month")
d18op_Alps_0 = extract_var(Dataset=Alps_0_data , varname="d18op", units="per mil", Dataset_wiso= Alps_0_wiso)
u10_Alps_0 = extract_var(Dataset=Alps_0_data , varname="u10")
v10_Alps_0 = extract_var(Dataset=Alps_0_data , varname="v10")




# contruct seasonal difference means
#Alps 150%
Alps_150_temp2_slt = compute_lterm_mean(data=temp2_Alps_150, time= "season", season_calendar="standard")
Alps_150_prec_slt = compute_lterm_mean(data=prec_Alps_150, time= "season", season_calendar="standard")
Alps_150_d18op_slt = compute_lterm_mean(data=d18op_Alps_150, time= "season", season_calendar="standard")
Alps_150_u10_slt = compute_lterm_mean(data=u10_Alps_150, time= "season", season_calendar="standard")
Alps_150_v10_slt = compute_lterm_mean(data=v10_Alps_150, time="season", season_calendar="standard")

#Alps 0%
Alps_0_temp2_slt = compute_lterm_mean(data=temp2_Alps_0, time= "season", season_calendar="standard")
Alps_0_prec_slt = compute_lterm_mean(data=prec_Alps_0, time= "season", season_calendar="standard")
Alps_0_d18op_slt = compute_lterm_mean(data=d18op_Alps_0, time= "season", season_calendar="standard")
Alps_0_u10_slt = compute_lterm_mean(data=u10_Alps_0, time= "season", season_calendar="standard")
Alps_0_v10_slt = compute_lterm_mean(data=v10_Alps_0, time="season", season_calendar="standard")

#Alps 100%
Alps_100_temp2_slt = compute_lterm_mean(data=temp2_Alps_100, time= "season", season_calendar="standard")
Alps_100_prec_slt = compute_lterm_mean(data=prec_Alps_100, time= "season", season_calendar="standard")
Alps_100_d18op_slt = compute_lterm_mean(data=d18op_Alps_100, time= "season", season_calendar="standard")
Alps_100_u10_slt = compute_lterm_mean(data=u10_Alps_100, time= "season", season_calendar="standard")
Alps_100_v10_slt = compute_lterm_mean(data=v10_Alps_100, time="season", season_calendar="standard")

# visualising long-term difference

#Alps 150% 
projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection":
                                                                                                                  projection})
plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=Alps_150_d18op_slt , cmap=RdBu, units="‰", seasons=["JJA", "DJF"], 
                   axes=[ax1,ax2], fig=fig, vmax=5, vmin=-20, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.69, 0.02, 0.25], title=True, 
                   season_label= ["[A]       Alps 150% [JJA]", "[B]        Alps 150% [DJF]"])

plot_seasonal_mean(variable="Temperature", data_slt=Alps_150_temp2_slt , cmap=RdBu_r, units="°C", seasons=["JJA", "DJF"], 
                   axes=[ax3,ax4], fig=fig, vmax=30, vmin=-15, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.38, 0.02, 0.25],title=True, season_label = ["[C]", "[D]"])

plot_seasonal_mean(variable="Precipitation", data_slt=Alps_150_prec_slt , cmap=Blues, units="mm/month", seasons=["JJA", "DJF"], 
                   axes=[ax5,ax6], fig=fig, vmax=350, vmin=0, levels=22, domain="Europe", level_ticks=6, cbar_pos = [0.90, 0.06, 0.02, 0.25],
                   title=True, data_u10=Alps_150_u10_slt, data_v10=Alps_150_v10_slt, season_label=["[E]", "[F]"])

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "figS5.png"), format= "png", bbox_inches="tight", dpi=600)

#Alps 0%
fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection":
                                                                                                                  projection})
plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=Alps_0_d18op_slt , cmap=RdBu, units="‰", seasons=["JJA", "DJF"], 
                    axes=[ax1,ax2], fig=fig, vmax=5, vmin=-20, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.69, 0.02, 0.25], title=True,
                    season_label= ["[A]       Alps 0% [JJA]", "[B]        Alps 0% [DJF]"])

plot_seasonal_mean(variable="Temperature", data_slt=Alps_0_temp2_slt , cmap=RdBu_r, units="°C", seasons=["JJA", "DJF"], 
                    axes=[ax3,ax4], fig=fig, vmax=30, vmin=-15, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.38, 0.02, 0.25],title=True,
                    season_label = ["[C]", "[D]"])

plot_seasonal_mean(variable="Precipitation", data_slt=Alps_0_prec_slt , cmap=Blues, units="mm/month", seasons=["JJA", "DJF"], 
                    axes=[ax5,ax6], fig=fig, vmax=350, vmin=0, levels=22, domain="Europe", level_ticks=6, cbar_pos = [0.90, 0.06, 0.02, 0.25],
                    title=True, data_u10=Alps_0_u10_slt, data_v10=Alps_0_v10_slt, season_label=["[E]", "[F]"])

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "figS4.png"), format= "png", bbox_inches="tight", dpi=600)


#Alps 100%
fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection":
                                                                                                                  projection})
plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=Alps_100_d18op_slt , cmap=RdBu, units="‰", seasons=["JJA", "DJF"], 
                    axes=[ax1,ax2], fig=fig, vmax=5, vmin=-20, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.69, 0.02, 0.25], title=True,
                    season_label= ["[A]        Alps 100% [JJA]", "[B]         Alps 100% [DJF]"])

plot_seasonal_mean(variable="Temperature", data_slt=Alps_100_temp2_slt , cmap=RdBu_r, units="°C", seasons=["JJA", "DJF"], 
                    axes=[ax3,ax4], fig=fig, vmax=30, vmin=-15, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.38, 0.02, 0.25],title=True,
                    season_label = ["[C]", "[D]"])

plot_seasonal_mean(variable="Precipitation", data_slt=Alps_100_prec_slt , cmap=Blues, units="mm/month", seasons=["JJA", "DJF"], 
                    axes=[ax5,ax6], fig=fig, vmax=350, vmin=0, levels=22, domain="Europe", level_ticks=6, cbar_pos = [0.90, 0.06, 0.02, 0.25],
                    title=True, data_u10=Alps_100_u10_slt, data_v10=Alps_100_v10_slt, season_label=["[E]", "[F]"])

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "figS3.png"), format= "png", bbox_inches="tight", dpi=600)
