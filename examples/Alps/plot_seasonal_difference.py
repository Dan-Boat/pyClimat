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

temp2_aw100e100_alt = compute_lterm_mean(data=temp2_aw100e100, time="annual")
prec_aw100e100_alt = compute_lterm_mean(data=prec_aw100e100,  time="annual")
d18op_aw100e100_alt = compute_lterm_mean(data=d18op_aw100e100, time="annual")
u10_aw100e100_alt = compute_lterm_mean(data=u10_aw100e100, time="annual")
v10_aw100e100_alt = compute_lterm_mean(data=v10_aw100e100, time="annual")

#aw100e0
u10_aw100e0_slt = compute_lterm_mean(data=u10_aw100e0, time="season", season_calendar="standard")
v10_aw100e0_slt = compute_lterm_mean(data=v10_aw100e0, time="season", season_calendar="standard")

u10_aw100e0_alt = compute_lterm_mean(data=u10_aw100e0, time="annual")
v10_aw100e0_alt = compute_lterm_mean(data=v10_aw100e0, time="annual")

#aw100e200
u10_aw100e200_slt = compute_lterm_mean(data=u10_aw100e200, time="season", season_calendar="standard")
v10_aw100e200_slt = compute_lterm_mean(data=v10_aw100e200, time="season", season_calendar="standard")

u10_aw100e200_alt = compute_lterm_mean(data=u10_aw100e200, time="annual")
v10_aw100e200_alt = compute_lterm_mean(data=v10_aw100e200, time="annual")

#aw200e100
u10_aw200e100_slt = compute_lterm_mean(data=u10_aw200e100, time="season", season_calendar="standard")
v10_aw200e100_slt = compute_lterm_mean(data=v10_aw200e100, time="season", season_calendar="standard")

u10_aw200e100_alt = compute_lterm_mean(data=u10_aw200e100, time="annual")
v10_aw200e100_alt = compute_lterm_mean(data=v10_aw200e100, time="annual")

#aw200e0
u10_aw200e0_slt = compute_lterm_mean(data=u10_aw200e0, time="season", season_calendar="standard")
v10_aw200e0_slt = compute_lterm_mean(data=v10_aw200e0, time="season", season_calendar="standard")

u10_aw200e0_alt = compute_lterm_mean(data=u10_aw200e0, time="annual")
v10_aw200e0_alt = compute_lterm_mean(data=v10_aw200e0, time="annual")

#aw200e200
u10_aw200e200_slt = compute_lterm_mean(data=u10_aw200e200, time="season", season_calendar="standard")
v10_aw200e200_slt = compute_lterm_mean(data=v10_aw200e200, time="season", season_calendar="standard")

u10_aw200e200_alt = compute_lterm_mean(data=u10_aw200e200, time="annual")
v10_aw200e200_alt = compute_lterm_mean(data=v10_aw200e200, time="annual")


# compute seasonal and annual difference lterm means

#aw100e0
temp2_aw100e0_slt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw100e0, time="season", season_calendar="standard")
prec_aw100e0_slt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw100e0, time="season", season_calendar="standard")
d18op_aw100e0_slt = compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw100e0, time="season", season_calendar="standard")

temp2_aw100e0_alt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw100e0, time="annual")
prec_aw100e0_alt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw100e0, time="annual")
d18op_aw100e0_alt = compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw100e0, time="annual")

#aw100e200
temp2_aw100e200_slt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw100e200, time="season", season_calendar="standard")
prec_aw100e200_slt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw100e200, time="season", season_calendar="standard")
d18op_aw100e200_slt= compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw100e200, time="season", season_calendar="standard")

temp2_aw100e200_alt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw100e200, time="annual")
prec_aw100e200_alt =  compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw100e200, time="annual")
d18op_aw100e200_alt =  compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw100e200, time="annual")

#aw200e100
temp2_aw200e100_slt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw200e100, time="season", season_calendar="standard")
prec_aw200e100_slt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw200e100, time="season", season_calendar="standard")
d18op_aw200e100_slt = compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw200e100, time="season", season_calendar="standard")

temp2_aw200e100_alt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw200e100, time="annual")
prec_aw200e100_alt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw200e100, time="annual")
d18op_aw200e100_alt = compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw200e100, time="annual")

#aw200e0
temp2_aw200e0_slt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw200e0, time="season", season_calendar="standard")
prec_aw200e0_slt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw200e0, time="season", season_calendar="standard")
d18op_aw200e0_slt = compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw200e0, time="season", season_calendar="standard")

temp2_aw200e0_alt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw200e0, time="annual")
prec_aw200e0_alt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw200e0, time="annual")
d18op_aw200e0_alt = compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw200e0, time="annual")

#aw200e200
temp2_aw200e200_slt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw200e200, time="season", season_calendar="standard")
prec_aw200e200_slt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw200e200, time="season", season_calendar="standard")
d18op_aw200e200_slt = compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw200e200, time="season", season_calendar="standard")

temp2_aw200e200_alt = compute_lterm_diff(data_control= temp2_aw100e100, data_main=temp2_aw200e200, time="annual")
prec_aw200e200_alt = compute_lterm_diff(data_control= prec_aw100e100, data_main=prec_aw200e200, time="annual")
d18op_aw200e200_alt = compute_lterm_diff(data_control= d18op_aw100e100, data_main=d18op_aw200e200, time="annual")



# visualising long-term difference

projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

#apply figure font and style
apply_style(fontsize=22, style=None, linewidth=2)

def plot_summer_diff(varname =None):
    if varname == "d18op":
        # d18op
        fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})
        
        
        plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=d18op_aw100e100_slt , cmap=YlGnBu, units="‰", seasons=["JJA"], 
                           axes=[ax1], fig=fig, vmax=2, vmin=-16, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.68, 0.02, 0.25], title=True, 
                           season_label= ["[A]  CTL"], bottom_labels=False, left_labels=True)
        
        plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=d18op_aw200e100_slt , cmap=RdBu, units="‰", seasons=["JJA"], 
                           axes=[ax2], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25], title=True, 
                           season_label= ["[B]  W2E1 - CTL"], bottom_labels=False, left_labels=False, compare_data1=d18op_aw100e100,
                           compare_data2=d18op_aw200e100, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=d18op_aw100e0_slt , cmap=RdBu, units="‰", seasons=["JJA"], 
                           axes=[ax3], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[C]  W1E0 - CTL"], bottom_labels=False, left_labels=True, compare_data1=d18op_aw100e100,
                           compare_data2=d18op_aw100e0, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=d18op_aw200e0_slt , cmap=RdBu, units="‰", seasons=["JJA"], 
                           axes=[ax4], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[D]  W2E0 - CTL"], bottom_labels=False, left_labels=False, compare_data1=d18op_aw100e100,
                           compare_data2=d18op_aw200e0, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=d18op_aw100e200_slt , cmap=RdBu, units="‰", seasons=["JJA"], 
                           axes=[ax5], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[E]  W1E2 - CTL"], bottom_labels=True, left_labels=True, compare_data1=d18op_aw100e100,
                           compare_data2=d18op_aw100e200, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=d18op_aw200e200_slt , cmap=RdBu, units="‰", seasons=["JJA"], 
                           axes=[ax6], fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[F]  W2E2 - CTL"], bottom_labels=True, left_labels=False, compare_data1=d18op_aw100e100,
                           compare_data2=d18op_aw200e200, max_pvalue=0.05, plot_stats=True)
        
        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout()
        plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
        plt.savefig(os.path.join(path_to_store, "fig2.svg"), format= "svg", bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(path_to_store, "fig2.png"), format= "png", bbox_inches="tight", dpi=300)
    
    
    if varname == "Temperature":
    # Temperature 
        fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})
        
        
        plot_seasonal_mean(variable="Temperature", data_slt=temp2_aw100e100_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                           axes=[ax1], fig=fig, vmax=30, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.68, 0.02, 0.25], title=True, 
                           season_label= ["[A]  CTL"], bottom_labels=False, left_labels=True)
        
        plot_seasonal_mean(variable="Temperature", data_slt=temp2_aw200e100_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                           axes=[ax2], fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25], title=True, 
                           season_label= ["[B]  W2E1 - CTL"], bottom_labels=False, left_labels=False, compare_data1=temp2_aw100e100,
                           compare_data2=temp2_aw200e100, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable="Temperature", data_slt=temp2_aw100e0_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                           axes=[ax3], fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[C]  W1E0 - CTL"], bottom_labels=False, left_labels=True, compare_data1=temp2_aw100e100,
                           compare_data2=temp2_aw100e0, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable="Temperature", data_slt=temp2_aw200e0_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                           axes=[ax4], fig=fig, vmax=13, vmin=-13, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[D]  W2E0 - CTL"], bottom_labels=False, left_labels=False, compare_data1=temp2_aw100e100,
                           compare_data2=temp2_aw200e0, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable="Temperature", data_slt=temp2_aw100e200_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                           axes=[ax5], fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[E]  W1E2 - CTL"], bottom_labels=True, left_labels=True, compare_data1=temp2_aw100e100,
                           compare_data2=temp2_aw100e200, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable="Temperature", data_slt=temp2_aw200e200_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                           axes=[ax6], fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[F]  W2E2 - CTL"], bottom_labels=True, left_labels=False, compare_data1=temp2_aw100e100,
                           compare_data2=temp2_aw200e200, max_pvalue=0.05, plot_stats=True)
    
    
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout()
        plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
        plt.savefig(os.path.join(path_to_store, "fig3.svg"), format= "svg", bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(path_to_store, "fig3.png"), format= "png", bbox_inches="tight", dpi=300)
    
    if varname == "Precipitation":
        #prec
        fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})
        
        
        plot_seasonal_mean(variable="Precipitation", data_slt=prec_aw100e100_slt , cmap=Blues, units="mm/month", seasons=["JJA"], 
                           axes=[ax1], fig=fig, vmax=250, vmin=0, levels=22, domain="Europe", level_ticks=6, cbar_pos = [0.90, 0.68, 0.02, 0.25], title=True, 
                           season_label= ["[A]  CTL"], bottom_labels=False, left_labels=True, plot_winds_pattern=True, data_v= v10_aw100e100_slt,
                           data_u=u10_aw100e100_slt)
        
        plot_seasonal_mean(variable="Precipitation", data_slt=prec_aw200e100_slt , cmap=BrBG, units="mm/month", seasons=["JJA"], 
                           axes=[ax2], fig=fig, vmax=125, vmin=-125, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25], title=True, 
                           season_label= ["[B]  W2E1 - CTL"], bottom_labels=False, left_labels=False, compare_data1=prec_aw100e100,
                           compare_data2=prec_aw200e100, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable="Precipitation", data_slt=prec_aw100e0_slt , cmap=BrBG, units="mm/month", seasons=["JJA"], 
                           axes=[ax3], fig=fig, vmax=125, vmin=-125, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[C]  W1E0 - CTL"], bottom_labels=False, left_labels=True, compare_data1=prec_aw100e100,
                           compare_data2=prec_aw100e0, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable="Precipitation", data_slt=prec_aw200e0_slt , cmap=BrBG, units="mm/month", seasons=["JJA"], 
                           axes=[ax4], fig=fig, vmax=125, vmin=-125, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[D]  W2E0 - CTL"], bottom_labels=False, left_labels=False, compare_data1=prec_aw100e100,
                           compare_data2=prec_aw200e0, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable="Precipitation", data_slt=prec_aw100e200_slt , cmap=BrBG, units="mm/month", seasons=["JJA"], 
                           axes=[ax5], fig=fig, vmax=125, vmin=-125, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[E]  W1E2 - CTL"], bottom_labels=True, left_labels=True, compare_data1=prec_aw100e100,
                           compare_data2=prec_aw100e200, max_pvalue=0.05, plot_stats=True)
        
        plot_seasonal_mean(variable="Precipitation", data_slt=prec_aw200e200_slt , cmap=BrBG, units="mm/month", seasons=["JJA"], 
                           axes=[ax6], fig=fig, vmax=125, vmin=-125, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                           season_label= ["[F]  W2E2 - CTL"], bottom_labels=True, left_labels=False, compare_data1=prec_aw100e100,
                           compare_data2=prec_aw200e200, max_pvalue=0.05, plot_stats=True)
        
        
        fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
        plt.tight_layout()
        plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
        plt.savefig(os.path.join(path_to_store, "fig4.svg"), format= "svg", bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(path_to_store, "fig4.png"), format= "png", bbox_inches="tight", dpi=300)


# supplementary for annual difference

def plot_annual_diff():

    # d18op
    fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})
    
    plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_aw100e100_alt, cmap=YlGnBu, units="‰", 
                       ax=ax1, fig=fig, vmax=2, vmin=-16, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.68, 0.02, 0.25],  
                       title= "[A]  CTL", bottom_labels=False, left_labels=True)
    
    plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_aw200e100_alt , cmap=RdBu, units="‰",  
                       ax=ax2, fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25],  
                       title= "[B]  W2E1 - CTL", bottom_labels=False, left_labels=False, compare_data1=d18op_aw100e100,
                       compare_data2=d18op_aw200e100, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_aw100e0_alt , cmap=RdBu, units="‰", 
                       ax=ax3, fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, add_colorbar=False,
                       title= "[C]  W1E0 - CTL", bottom_labels=False, left_labels=True, compare_data1=d18op_aw100e100,
                       compare_data2=d18op_aw100e0, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_aw200e0_alt , cmap=RdBu, units="‰",  
                       ax=ax4, fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, add_colorbar=False,
                       title= "[D]  W2E0 - CTL", bottom_labels=False, left_labels=False, compare_data1=d18op_aw100e100,
                       compare_data2=d18op_aw200e0, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_aw100e200_alt , cmap=RdBu, units="‰",
                       ax=ax5, fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, add_colorbar=False,
                       title= "[E]  W1E2 - CTL", bottom_labels=True, left_labels=True, compare_data1=d18op_aw100e100,
                       compare_data2=d18op_aw100e200, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable='$\delta^{18}$Op vs SMOW', data_alt=d18op_aw200e200_alt , cmap=RdBu, units="‰",
                       ax=ax6, fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, add_colorbar=False,
                       title= "[F]  W2E2 - CTL", bottom_labels=True, left_labels=False, compare_data1=d18op_aw100e100,
                       compare_data2=d18op_aw200e200, max_pvalue=0.05, plot_stats=True)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
    plt.savefig(os.path.join(path_to_store, "figS5.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "figS5.png"), format= "png", bbox_inches="tight", dpi=300)
    
    
    
    # Temperature 
    fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})
    
    
    plot_annual_mean(variable="Temperature", data_alt=temp2_aw100e100_alt , cmap=RdBu_r, units="°C",  
                       ax=ax1, fig=fig, vmax=30, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.68, 0.02, 0.25], 
                       title= "[A]  CTL", bottom_labels=False, left_labels=True)
    
    plot_annual_mean(variable="Temperature", data_alt=temp2_aw200e100_alt , cmap=RdBu_r, units="°C",  
                       ax=ax2, fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25],  
                       title= "[B]  W2E1 - CTL", bottom_labels=False, left_labels=False, compare_data1=temp2_aw100e100,
                       compare_data2=temp2_aw200e100, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable="Temperature", data_alt=temp2_aw100e0_alt , cmap=RdBu_r, units="°C",  
                       ax=ax3, fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11,  add_colorbar=False,
                       title= "[C]  W1E0 - CTL", bottom_labels=False, left_labels=True, compare_data1=temp2_aw100e100,
                       compare_data2=temp2_aw100e0, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable="Temperature", data_alt=temp2_aw200e0_alt , cmap=RdBu_r, units="°C",  
                       ax=ax4, fig=fig, vmax=14, vmin=-14, levels=22, domain="Europe", level_ticks=11,  add_colorbar=False,
                       title= "[D]  W2E0 - CTL", bottom_labels=False, left_labels=False, compare_data1=temp2_aw100e100,
                       compare_data2=temp2_aw200e0, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable="Temperature", data_alt=temp2_aw100e200_alt , cmap=RdBu_r, units="°C", 
                       ax=ax5, fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11,  add_colorbar=False,
                       title= "[E]  W1E2 - CTL", bottom_labels=True, left_labels=True, compare_data1=temp2_aw100e100,
                       compare_data2=temp2_aw100e200, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable="Temperature", data_alt=temp2_aw200e200_alt , cmap=RdBu_r, units="°C", 
                       ax=ax6, fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11,  add_colorbar=False,
                       title= "[F]  W2E2 - CTL", bottom_labels=True, left_labels=False, compare_data1=temp2_aw100e100,
                       compare_data2=temp2_aw200e200, max_pvalue=0.05, plot_stats=True)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
    plt.savefig(os.path.join(path_to_store, "figS6.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "figS6.png"), format= "png", bbox_inches="tight", dpi=300)
    
    
    #prec
    fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})
    
    
    plot_annual_mean(variable="Precipitation", data_alt=prec_aw100e100_alt , cmap=Blues, units="mm/month",
                       ax=ax1, fig=fig, vmax=300, vmin=0, levels=22, domain="Europe", level_ticks=6, cbar_pos = [0.90, 0.68, 0.02, 0.25], 
                       title= "[A]  CTL", bottom_labels=False, left_labels=True, data_v10= v10_aw100e100_alt,
                       data_u10=u10_aw100e100_alt)
    
    plot_annual_mean(variable="Precipitation", data_alt=prec_aw200e100_alt , cmap=BrBG, units="mm/month",  
                       ax=ax2, fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25], 
                       title= "[B]  W2E1 - CTL", bottom_labels=False, left_labels=False, compare_data1=prec_aw100e100,
                       compare_data2=prec_aw200e100, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable="Precipitation", data_alt=prec_aw100e0_alt , cmap=BrBG, units="mm/month",  
                       ax=ax3, fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=11,  add_colorbar=False,
                       title= "[C]  W1E0 - CTL", bottom_labels=False, left_labels=True, compare_data1=prec_aw100e100,
                       compare_data2=prec_aw100e0, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable="Precipitation", data_alt=prec_aw200e0_alt , cmap=BrBG, units="mm/month",
                       ax=ax4, fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=11,  add_colorbar=False,
                       title= "[D]  W2E0 - CTL", bottom_labels=False, left_labels=False, compare_data1=prec_aw100e100,
                       compare_data2=prec_aw200e0, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable="Precipitation", data_alt=prec_aw100e200_alt , cmap=BrBG, units="mm/month", 
                       ax=ax5, fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=11,  add_colorbar=False,
                       title= "[E]  W1E2 - CTL", bottom_labels=True, left_labels=True, compare_data1=prec_aw100e100,
                       compare_data2=prec_aw100e200, max_pvalue=0.05, plot_stats=True)
    
    plot_annual_mean(variable="Precipitation", data_alt=prec_aw200e200_alt , cmap=BrBG, units="mm/month",  
                       ax=ax6, fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=11,  add_colorbar=False,
                       title= "[F]  W2E2 - CTL", bottom_labels=True, left_labels=False, compare_data1=prec_aw100e100,
                       compare_data2=prec_aw200e200, max_pvalue=0.05, plot_stats=True)
    
    
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
    plt.savefig(os.path.join(path_to_store, "figS7.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "figS7.png"), format= "png", bbox_inches="tight", dpi=300)

    
if __name__ == '__main__':
    plot_summer_diff(varname="Temperature")
    plot_summer_diff(varname="d18op")
    plot_summer_diff(varname="Precipitation")
    
    # # annual plots 
    #plot_annual_diff()
