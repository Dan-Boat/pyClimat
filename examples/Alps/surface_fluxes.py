#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 19:13:46 2022

@author: dboateng
"""

#importing modules 
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *


# defining module paths
module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east_200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"

years= "1003_1017"
period = "1m"

#dataset
control_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period, read_wiso=False)
east_0_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_0, years=years,
                                                  period=period, read_wiso=False)
east_200_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_200, years=years,
                                                  period=period, read_wiso=False)

#latent heat
ahfl = extract_var(Dataset=control_data , varname="latent heat")
ahfl_east_0 = extract_var(Dataset=east_0_data , varname="latent heat")
ahfl_east_200 = extract_var(Dataset=east_200_data , varname="latent heat")

#lterm means
ahfl_slt = compute_lterm_mean(data=ahfl, time="season", season_calendar="standard")
# ahfl_east_0_slt = compute_lterm_diff(data_control=ahfl, data_main=ahfl_east_0, time="season", season_calendar="standard")
# ahfl_east_200_slt = compute_lterm_diff(data_control=ahfl, data_main=ahfl_east_200, time="season", season_calendar="standard")


ahfl_east_0_slt = compute_lterm_mean(data=ahfl_east_0, time="season", season_calendar="standard")
ahfl_east_200_slt = compute_lterm_mean(data=ahfl_east_200, time="season", season_calendar="standard")


# E-P history 
e_p = extract_var(Dataset=control_data , varname="E-P", units="mm/month")
e_p_east_0 = extract_var(Dataset=east_0_data , varname="E-P", units="mm/month")
e_p_east_200 = extract_var(Dataset=east_200_data , varname="E-P", units="mm/month")

#lterm means
e_p_slt = compute_lterm_mean(data=e_p, time="season", season_calendar="standard")
e_p_east_0_slt = compute_lterm_mean(data=e_p_east_0, time="season", season_calendar="standard")
e_p_east_200_slt = compute_lterm_mean(data=e_p_east_200, time="season", season_calendar="standard")


projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(nrows = 2, ncols = 3, figsize=(17, 12), subplot_kw={"projection": projection})

plot_seasonal_mean(variable="Latent heat flux", data_slt=ahfl_slt , cmap=RdYlBu_r, units="W/m²", seasons=["JJA"], 
                    axes=[ax1], fig=fig, vmax=120, vmin=0, levels=22, domain="Europe", level_ticks=6, cbar_pos = [0.90, 0.60, 0.02, 0.25], title=True, 
                    bottom_labels=False, left_labels=True, season_label=["[A]          AW100E100"])

plot_seasonal_mean(variable="Latent heat flux", data_slt=ahfl_east_0_slt , cmap=RdYlBu_r, units="W/m²", seasons=["JJA"], 
                    axes=[ax2], vmax=120, vmin=0, levels=22, domain="Europe", level_ticks=6, title=True, 
                    bottom_labels=False, left_labels=False, season_label=["[B]          AW100E0"], add_colorbar=False)

plot_seasonal_mean(variable="Latent heat flux", data_slt=ahfl_east_200_slt , cmap=RdYlBu_r, units="W/m²", seasons=["JJA"], 
                    axes=[ax3], vmax=120, vmin=0, levels=22, domain="Europe", level_ticks=6, title=True, 
                    bottom_labels=False, left_labels=False, season_label=["[C]          AW100E200"], add_colorbar=False)

plot_seasonal_mean(variable="E-P", data_slt=e_p_slt , cmap=RdBu_r, units="mm/month", seasons=["JJA"], 
                    axes=[ax4], fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=6, cbar_pos = [0.90, 0.05, 0.02, 0.25], title=True, 
                    bottom_labels=False, left_labels=True, season_label=["[D]"])

plot_seasonal_mean(variable="E-P", data_slt=e_p_east_0_slt , cmap=RdBu_r, units="mm/month", seasons=["JJA"], 
                    axes=[ax5], fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=6, title=True, 
                    bottom_labels=False, left_labels=False, season_label=["[E]"], add_colorbar=False)

plot_seasonal_mean(variable="E-P", data_slt=e_p_east_200_slt , cmap=RdBu_r, units="mm/month", seasons=["JJA"], 
                    axes=[ax6], fig=fig, vmax=100, vmin=-100, levels=22, domain="Europe", level_ticks=6, title=True, 
                    bottom_labels=False, left_labels=False, season_label=["[F]"], add_colorbar=False)

