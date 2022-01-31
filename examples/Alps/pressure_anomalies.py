#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 14:09:16 2022

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
aw100e100_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e100, years=years,
                                                  period=period, read_wiso=False, add_name="slp_JJA")
aw100e0_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e0, years=years,
                                                  period=period, read_wiso=False, add_name="slp_JJA")
aw100e200_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e200, years=years,
                                                  period=period, read_wiso=False, add_name="slp_JJA")
aw200e100_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e100, years=years,
                                                  period=period, read_wiso=False, add_name="slp_JJA")
aw200e0_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e0, years=years,
                                                  period=period, read_wiso=False, add_name="slp_JJA")
aw200e200_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e200, years=years,
                                                  period=period, read_wiso=False, add_name="slp_JJA")

#  extracting variables
 
#aw100e100
aps_aw100e100 = extract_var(Dataset=aw100e100_data , varname="slp", units="hPa")
aps_aw100e0 = extract_var(Dataset=aw100e0_data , varname="slp", units="hPa")
aps_aw100e200 = extract_var(Dataset=aw100e200_data , varname="slp", units="hPa")
aps_aw200e100 = extract_var(Dataset=aw200e100_data , varname="slp", units="hPa")
aps_aw200e0 = extract_var(Dataset=aw200e0_data , varname="slp", units="hPa")
aps_aw200e200 = extract_var(Dataset=aw200e200_data , varname="slp", units="hPa")

aps_aw100e100_slt = compute_lterm_mean(data=aps_aw100e100 )
aps_aw100e0_slt = compute_lterm_mean(data=aps_aw100e0)
aps_aw100e200_slt = compute_lterm_diff(data_control= aps_aw100e100, data_main=aps_aw100e200,)
aps_aw200e100_slt = compute_lterm_diff(data_control= aps_aw100e100, data_main=aps_aw200e100,)
aps_aw200e0_slt = compute_lterm_diff(data_control= aps_aw100e100, data_main=aps_aw200e0, )
aps_aw200e200_slt = compute_lterm_diff(data_control= aps_aw100e100, data_main=aps_aw200e200,)


projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")


fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})


plot_annual_mean(variable="aps", data_alt=aps_aw100e100_slt , cmap=RdBu_r, units="hPa", 
                   ax=ax1, fig=fig, vmax=1030, vmin=850, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.68, 0.02, 0.25],  
                   title= ["[A]  AW100E100"], bottom_labels=False, left_labels=True)

plot_annual_mean(variable="aps", data_alt=aps_aw200e100_slt , cmap=RdBu_r, units="hPa",
                   ax=ax2, fig=fig, vmax=10, vmin=-10, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25],  
                   title= ["[B]  AW200E100 - AW100E100"], bottom_labels=False, left_labels=False, compare_data1=aps_aw100e100,
                   compare_data2=aps_aw200e100, max_pvalue=0.05, plot_stats=True)

plot_annual_mean(variable="aps", data_alt=aps_aw100e0_slt , cmap=RdBu_r, units="hPa",
                   ax=ax3, fig=fig, vmax=1030, vmin=850, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.25, 0.02, 0.25],  
                   title= ["[B]  AW100E0 - AW100E100"], bottom_labels=False, left_labels=False, compare_data1=aps_aw100e100,
                   compare_data2=aps_aw200e100, max_pvalue=0.05, plot_stats=True)


plot_seasonal_mean(variable="aps", data_slt=aps_aw100e0_slt , cmap=RdBu_r, units="hPa", seasons=["JJA"], 
                   axes=[ax3], fig=fig, vmax=5, vmin=-5, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                   season_label= ["[C]  AW100E0 - AW100E100"], bottom_labels=False, left_labels=True, compare_data1=aps_aw100e100,
                   compare_data2=aps_aw100e0, max_pvalue=0.05, plot_stats=True)

plot_seasonal_mean(variable="aps", data_slt=aps_aw200e0_slt , cmap=RdBu_r, units="hPa", seasons=["JJA"], 
                   axes=[ax4], fig=fig, vmax=5, vmin=-5, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                   season_label= ["[D]  AW200E0 - AW100E100"], bottom_labels=False, left_labels=False, compare_data1=aps_aw100e100,
                   compare_data2=aps_aw200e0, max_pvalue=0.05, plot_stats=True)

plot_seasonal_mean(variable="aps", data_slt=aps_aw100e200_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                   axes=[ax5], fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                   season_label= ["[E]  AW100E200 - AW100E100"], bottom_labels=True, left_labels=True, compare_data1=aps_aw100e100,
                   compare_data2=aps_aw100e200, max_pvalue=0.05, plot_stats=True)

plot_seasonal_mean(variable="aps", data_slt=aps_aw200e200_slt , cmap=RdBu_r, units="°C", seasons=["JJA"], 
                   axes=[ax6], fig=fig, vmax=12, vmin=-12, levels=22, domain="Europe", level_ticks=11, title=True, add_colorbar=False,
                   season_label= ["[F]  AW200E200 - AW100E100"], bottom_labels=True, left_labels=False, compare_data1=aps_aw100e100,
                   compare_data2=aps_aw200e200, max_pvalue=0.05, plot_stats=True)


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "fig3.svg"), format= "svg", bbox_inches="tight", dpi=600)
