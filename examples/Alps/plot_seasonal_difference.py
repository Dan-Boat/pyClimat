#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:41:41 2021

@author: dboateng
This is example script for using Climat in visualising ECHAM module output (long-term seasonal differece means)
The script contains directories of module outputs and path to save plot
Note: it is structured solely for the personal needs of the author, therefore, it must be adapted advisably.
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)


import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east_200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"

# for supplementary
exp_name_east_150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"
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
east_50_data, east_50_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_50, years=years,
                                                  period=period)
east_200_data, east_200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_200, years=years,
                                                  period=period)

#  extracting variables 
#control
temp2 = extract_var(Dataset=control_data , varname="temp2", units="°C")
prec = extract_var(Dataset= control_data , varname="prec", units="mm/month")
d18op = extract_var(Dataset=control_data , varname="d18op", units="per mil", Dataset_wiso= control_wiso)
u10 = extract_var(Dataset=control_data , varname="u10")
v10 = extract_var(Dataset=control_data , varname="v10")

#east_150
temp2_east_150 = extract_var(Dataset=east_150_data , varname="temp2", units="°C")
prec_east_150 = extract_var(Dataset= east_150_data , varname="prec", units="mm/month")
d18op_east_150 = extract_var(Dataset=east_150_data , varname="d18op", units="per mil", Dataset_wiso= east_150_wiso)
u10_east_150 = extract_var(Dataset=east_150_data , varname="u10")
v10_east_150 = extract_var(Dataset=east_150_data , varname="v10")

#east_0
temp2_east_0 = extract_var(Dataset=east_0_data , varname="temp2", units="°C")
prec_east_0 = extract_var(Dataset= east_0_data , varname="prec", units="mm/month")
d18op_east_0 = extract_var(Dataset=east_0_data , varname="d18op", units="per mil", Dataset_wiso= east_0_wiso)
u10_east_0 = extract_var(Dataset=east_0_data , varname="u10")
v10_east_0 = extract_var(Dataset=east_0_data , varname="v10")

#east_50
temp2_east_50 = extract_var(Dataset=east_50_data , varname="temp2", units="°C")
prec_east_50 = extract_var(Dataset= east_50_data , varname="prec", units="mm/month")
d18op_east_50 = extract_var(Dataset=east_50_data , varname="d18op", units="per mil", Dataset_wiso= east_50_wiso)
u10_east_50 = extract_var(Dataset=east_50_data , varname="u10")
v10_east_50 = extract_var(Dataset=east_50_data , varname="v10")


#east_200
temp2_east_200 = extract_var(Dataset=east_200_data , varname="temp2", units="°C")
prec_east_200 = extract_var(Dataset= east_200_data , varname="prec", units="mm/month")
d18op_east_200 = extract_var(Dataset=east_200_data , varname="d18op", units="per mil", Dataset_wiso= east_200_wiso)
u10_east_200 = extract_var(Dataset=east_200_data , varname="u10")
v10_east_200 = extract_var(Dataset=east_200_data , varname="v10")

# compute long-term means for control experiment

# seasonal and annual
Alps_100_temp2_slt = compute_lterm_mean(data=temp2, time="season", season_calendar="standard")
Alps_100_prec_slt = compute_lterm_mean(data=prec, time="season", season_calendar="standard")
Alps_100_d18op_slt = compute_lterm_mean(data=d18op , time="season", season_calendar="standard")
Alps_100_u10_slt = compute_lterm_mean(data=u10 , time="season", season_calendar="standard")
Alps_100_v10_slt = compute_lterm_mean(data=v10 , time="season", season_calendar="standard")

Alps_100_temp2_alt = compute_lterm_mean(data=temp2, time="annual")
Alps_100_prec_alt = compute_lterm_mean(data=prec, time="annual")
Alps_100_d18op_alt = compute_lterm_mean(data=d18op , time="annual")
Alps_100_u10_alt = compute_lterm_mean(data=u10 , time="annual")
Alps_100_v10_alt = compute_lterm_mean(data=v10 , time="annual")

# compute seasonal and annual difference means


#east_150
east_150_temp2_slt = compute_lterm_diff(data_control= temp2, data_main=temp2_east_150, time="season", season_calendar="standard")
east_150_prec_slt = compute_lterm_diff(data_control=prec, data_main=prec_east_150, time="season", season_calendar="standard")
east_150_d18op_slt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_150, time="season", season_calendar="standard")
east_150_u10_slt = compute_lterm_mean(data=u10_east_150, time="season", season_calendar="standard")
east_150_v10_slt = compute_lterm_mean(data=v10_east_150, time="season", season_calendar="standard")

east_150_temp2_alt = compute_lterm_diff(data_control= temp2, data_main=temp2_east_150, time="annual")
east_150_prec_alt = compute_lterm_diff(data_control=prec, data_main=prec_east_150, time="annual")
east_150_d18op_alt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_150, time="annual")
east_150_u10_alt = compute_lterm_mean(data=u10_east_150, time="annual")
east_150_v10_alt = compute_lterm_mean(data=v10_east_150, time="annual")

#east_0
east_0_temp2_slt = compute_lterm_diff(data_control= temp2, data_main=temp2_east_0, time="season", season_calendar="standard")
east_0_prec_slt = compute_lterm_diff(data_control=prec, data_main=prec_east_0, time="season", season_calendar="standard")
east_0_d18op_slt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_0, time="season", season_calendar="standard")
east_0_u10_slt = compute_lterm_mean(data=u10_east_0, time="season", season_calendar="standard")
east_0_v10_slt = compute_lterm_mean(data=v10_east_0, time="season", season_calendar="standard")

east_0_temp2_alt = compute_lterm_diff(data_control= temp2, data_main=temp2_east_0, time="annual",)
east_0_prec_alt = compute_lterm_diff(data_control=prec, data_main=prec_east_0, time="annual",)
east_0_d18op_alt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_0, time="annual",)
east_0_u10_alt = compute_lterm_mean(data=u10_east_0, time="annual",)
east_0_v10_alt = compute_lterm_mean(data=v10_east_0, time="annual",)

#east_50
east_50_temp2_slt = compute_lterm_diff(data_control= temp2, data_main=temp2_east_50, time="season", season_calendar="standard")
east_50_prec_slt = compute_lterm_diff(data_control=prec, data_main=prec_east_50, time="season", season_calendar="standard")
east_50_d18op_slt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_50, time="season", season_calendar="standard")
east_50_u10_slt = compute_lterm_mean(data=u10_east_50, time="season", season_calendar="standard")
east_50_v10_slt = compute_lterm_mean(data=v10_east_50, time="season", season_calendar="standard")

east_50_temp2_alt = compute_lterm_diff(data_control= temp2, data_main=temp2_east_50, time="annual")
east_50_prec_alt = compute_lterm_diff(data_control=prec, data_main=prec_east_50, time="annual")
east_50_d18op_alt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_50, time="annual")
east_50_u10_alt = compute_lterm_mean(data=u10_east_50, time="annual")
east_50_v10_alt = compute_lterm_mean(data=v10_east_50, time="annual")

#east_200
east_200_temp2_slt = compute_lterm_diff(data_control= temp2, data_main=temp2_east_200, time="season", season_calendar="standard")
east_200_prec_slt = compute_lterm_diff(data_control=prec, data_main=prec_east_200, time="season", season_calendar="standard")
east_200_d18op_slt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_200, time="season", season_calendar="standard")
east_200_u10_slt = compute_lterm_mean(data=u10_east_200, time="season", season_calendar="standard")
east_200_v10_slt = compute_lterm_mean(data=v10_east_200, time="season", season_calendar="standard")

east_200_temp2_alt = compute_lterm_diff(data_control= temp2, data_main=temp2_east_200, time="annual")
east_200_prec_alt = compute_lterm_diff(data_control=prec, data_main=prec_east_200, time="annual")
east_200_d18op_alt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_200, time="annual")
east_200_u10_alt = compute_lterm_mean(data=u10_east_200, time="annual")
east_200_v10_alt = compute_lterm_mean(data=v10_east_200, time="annual")

# visualising long-term difference

projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")


# d18op
fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})


plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=Alps_100_d18op_slt , cmap=YlGnBu, units="‰", seasons=["JJA"], 
                   axes=[ax1], fig=fig, vmax=2, vmin=-16, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.69, 0.02, 0.25], title=True, 
                   season_label= ["[A]          Alps 100% [JJA]"], bottom_labels=False, left_labels=True)

plot_annual_mean(ax=ax2, variable='$\delta^{18}$Op vs SMOW', data_alt=Alps_100_d18op_alt, cmap=YlGnBu, units="‰", vmax=2, vmin=-16, domain="Europe", 
                 levels=22, level_ticks=11 , title="[B]          Alps 100%", left_labels=False, bottom_labels=False,add_colorbar=False)


plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=east_200_d18op_slt , cmap=RdBu, units="‰", seasons=["JJA"], 
                   axes=[ax3], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25], title=True, 
                   season_label= ["[C]          Alps East 200% - Alps 100% [JJA]"], compare_data1=d18op, compare_data2=d18op_east_200, max_pvalue=0.05, plot_stats=True, 
                   bottom_labels=False)


plot_annual_mean(ax=ax4, variable='$\delta^{18}$Op vs SMOW', data_alt=east_200_d18op_alt, cmap=RdBu, units="‰", vmax=10, vmin=-10, domain="Europe", 
                 levels=22, level_ticks=11 , title="[D]          Alps East 200% - Alps 100%", left_labels=False, bottom_labels=False,add_colorbar=False, compare_data1=d18op,
                 compare_data2=d18op_east_200, max_pvalue=0.05, plot_stats=True)


plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=east_0_d18op_slt , cmap=RdBu, units="‰", seasons=["JJA"], 
                    axes=[ax5], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, title=True,
                    season_label= ["[E]          Alps East 0% - Alps 100% [JJA]"], compare_data1=d18op, compare_data2=d18op_east_0, max_pvalue=0.05, plot_stats=True, 
                    add_colorbar=False, bottom_labels=True, left_labels=True)

plot_annual_mean(ax=ax6, variable='$\delta^{18}$Op vs SMOW', data_alt=east_0_d18op_alt, cmap=RdBu, units="‰", vmax=10, vmin=-10, domain="Europe", 
                 levels=22, level_ticks=11 , title="[F]          Alps East 0% - Alps 100%", left_labels = False, bottom_labels=True, add_colorbar=False, compare_data1=d18op,
                 compare_data2=d18op_east_0, max_pvalue=0.05, plot_stats=True)



fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "fig2.svg"), format= "svg", bbox_inches="tight", dpi=600)



# Temperature
fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})


plot_seasonal_mean(variable="Temperature", data_slt=Alps_100_temp2_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                   axes=[ax1], fig=fig, vmax=30, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.69, 0.02, 0.25], title=True, 
                   season_label= ["[A]          Alps 100% [JJA]"], bottom_labels=False, left_labels=True)

plot_annual_mean(ax=ax2, variable="Temperature", data_alt=Alps_100_temp2_alt, cmap=RdBu_r, units="°C", vmax=30, vmin=-10, domain="Europe", 
                 levels=22, level_ticks=11 , title="[B]          Alps 100%", left_labels=False, bottom_labels=False,add_colorbar=False)


plot_seasonal_mean(variable="Temperature", data_slt=east_200_temp2_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                   axes=[ax3], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25], title=True, 
                   season_label= ["[C]          Alps East 200% - Alps 100% [JJA]"], compare_data1=temp2, compare_data2=temp2_east_200, max_pvalue=0.05, plot_stats=True, 
                   bottom_labels=False)


plot_annual_mean(ax=ax4, variable="Temperature", data_alt=east_200_temp2_alt, cmap=RdBu_r, units="°C", vmax=10, vmin=-10, domain="Europe", 
                 levels=22, level_ticks=11 , title="[D]          Alps East 200% - Alps 100%", left_labels=False, bottom_labels=False,add_colorbar=False, compare_data1=temp2,
                 compare_data2=temp2_east_200, max_pvalue=0.05, plot_stats=True)


plot_seasonal_mean(variable="Temperature", data_slt=east_0_temp2_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                    axes=[ax5], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, title=True,
                    season_label= ["[E]          Alps East 0% - Alps 100% [JJA]"], compare_data1=temp2, compare_data2=temp2_east_0, max_pvalue=0.05, plot_stats=True, 
                    add_colorbar=False, bottom_labels=True, left_labels=True)

plot_annual_mean(ax=ax6, variable="Temperature", data_alt=east_0_temp2_alt, cmap=RdBu_r, units="°C", vmax=10, vmin=-10, domain="Europe", 
                 levels=22, level_ticks=11 , title="[F]          Alps East 0% - Alps 100%", left_labels = False, bottom_labels=True, add_colorbar=False, compare_data1=temp2,
                 compare_data2=temp2_east_0, max_pvalue=0.05, plot_stats=True)



fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "fig3.svg"), format= "svg", bbox_inches="tight", dpi=600)



# Precipitation 
fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})


plot_seasonal_mean(variable="Precipitation", data_slt=Alps_100_prec_slt , cmap=Blues, units="mm/month", seasons=["JJA"], 
                   axes=[ax1], fig=fig, vmax=300, vmin=0, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.69, 0.02, 0.25], title=True, 
                   season_label= ["[A]          Alps 100% [JJA]"], bottom_labels=False, left_labels=True, plot_winds_pattern=True, data_v=Alps_100_v10_slt,
                   data_u= Alps_100_u10_slt)

plot_annual_mean(ax=ax2, variable="Precipitation", data_alt=Alps_100_prec_alt, cmap=Blues, units="mm/month", vmax=300, vmin=0, domain="Europe", 
                 levels=22, level_ticks=11 , title="[B]          Alps 100%", left_labels=False, bottom_labels=False,add_colorbar=False, data_v10=Alps_100_v10_alt,
                 data_u10=Alps_100_u10_alt,)


plot_seasonal_mean(variable="Precipitation", data_slt=east_200_prec_slt , cmap=BrBG, units="mm/month", seasons=["JJA"], 
                   axes=[ax3], fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25], title=True, 
                   season_label= ["[C]          Alps East 200% - Alps 100% [JJA]"], compare_data1=prec, compare_data2=prec_east_200, max_pvalue=0.05, plot_stats=True, 
                   bottom_labels=False)


plot_annual_mean(ax=ax4, variable="Precipitation", data_alt=east_200_prec_alt, cmap=BrBG, units="mm/month", vmax=100, vmin=-100, domain="Europe", 
                 levels=22, level_ticks=11 , title="[D]          Alps East 200% - Alps 100%", left_labels=False, bottom_labels=False,add_colorbar=False, compare_data1=prec,
                 compare_data2=prec_east_200, max_pvalue=0.05, plot_stats=True)


plot_seasonal_mean(variable="Precipitation", data_slt=east_0_prec_slt , cmap=BrBG, units="mm/month", seasons=["JJA"], 
                    axes=[ax5], fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=11, title=True,
                    season_label= ["[E]          Alps East 0% - Alps 100% [JJA]"], compare_data1=prec, compare_data2=prec_east_0, max_pvalue=0.05, plot_stats=True, 
                    add_colorbar=False, bottom_labels=True, left_labels=True)

plot_annual_mean(ax=ax6, variable="Precipitation", data_alt=east_0_prec_alt, cmap=BrBG, units="mm/month", vmax=100, vmin=-100, domain="Europe", 
                 levels=22, level_ticks=11 , title="[F]          Alps East 0% - Alps 100%", left_labels = False, bottom_labels=True, add_colorbar=False, compare_data1=prec,
                 compare_data2=prec_east_0, max_pvalue=0.05, plot_stats=True)



fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "fig4.svg"), format= "svg", bbox_inches="tight", dpi=600)

# implementation for supplementary figures
