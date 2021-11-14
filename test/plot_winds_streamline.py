#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 12:33:41 2021

@author: dboateng

This script visualise the winds at 850 hPa for all experiments with elevation as background plot
"""

# importing package

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name_Alps100 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east150 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
years= "1003_1017"
period = "1m"



# reading winds and elevation data. 

# Alps 100%
data_v = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps100, years=years,
                                                  period=period, add_name="v", read_wiso=False)
data_u = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps100, years=years,
                                                  period=period, add_name="u", read_wiso=False)
data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps100, years=years,
                                                  period=period, read_wiso=False)
elev = extract_var(Dataset=data , varname="elev", units="m")
v = extract_var(Dataset=data_v, varname="v", lev_units="hPa", lev=850)
u = extract_var(Dataset=data_u, varname="u", lev_units="hPa", lev=850)
u10 = extract_var(Dataset=data , varname="u10")
v10 = extract_var(Dataset=data , varname="v10")
ep = extract_var(Dataset=data , varname="e/p")



# Alps east 0%
data_v_east0 = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east0, years=years,
                                                  period=period, add_name="v", read_wiso=False)
data_u_east0 = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east0, years=years,
                                                  period=period, add_name="u", read_wiso=False)
data_east0 = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east0, years=years,
                                                  period=period, read_wiso=False)
elev_east0 = extract_var(Dataset=data_east0 , varname="elev", units="m")
v_east0 = extract_var(Dataset=data_v_east0, varname="v", lev_units="hPa", lev=850)
u_east0 = extract_var(Dataset=data_u_east0, varname="u", lev_units="hPa", lev=850)
u10_east0 = extract_var(Dataset=data_east0 , varname="u10")
v10_east0 = extract_var(Dataset=data_east0 , varname="v10")


# Alps east 150%
data_v_east150 = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east150, years=years,
                                                  period=period, add_name="v", read_wiso=False)
data_u_east150 = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east150, years=years,
                                                  period=period, add_name="u", read_wiso=False)
data_east150 = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east150, years=years,
                                                  period=period, read_wiso=False)
elev_east150 = extract_var(Dataset=data_east150, varname="elev", units="m")
v_east150 = extract_var(Dataset=data_v_east150, varname="v", lev_units="hPa", lev=850)
u_east150 = extract_var(Dataset=data_u_east150, varname="u", lev_units="hPa", lev=850)
u10_east150 = extract_var(Dataset=data_east150 , varname="u10")
v10_east150 = extract_var(Dataset=data_east150 , varname="v10")


# computing seasonal long-term means
elev_slt = compute_lterm_mean(data=elev, time="season", season_calendar="standard")
v_slt = compute_lterm_mean(data=v, time="season", season_calendar="standard")
u_slt = compute_lterm_mean(data=u, time="season", season_calendar="standard")
u10_slt = compute_lterm_mean(data=u10, time="season", season_calendar="standard")
v10_slt = compute_lterm_mean(data=v10, time="season", season_calendar="standard")
ep_slt = compute_lterm_mean(data=ep, time="season", season_calendar="standard")


elev_east0_slt = compute_lterm_mean(data=elev_east0, time="season", season_calendar="standard")
v_east0_slt = compute_lterm_mean(data=v_east0, time="season", season_calendar="standard")
u_east0_slt = compute_lterm_mean(data=u_east0, time="season", season_calendar="standard")
u10_east0_slt = compute_lterm_mean(data=u10_east0, time="season", season_calendar="standard")
v10_east0_slt = compute_lterm_mean(data=v10_east0, time="season", season_calendar="standard")

elev_east150_slt = compute_lterm_mean(data=elev_east150, time="season", season_calendar="standard")
v_east150_slt = compute_lterm_mean(data=v_east150, time="season", season_calendar="standard")
u_east150_slt = compute_lterm_mean(data=u_east150, time="season", season_calendar="standard")
u10_east150_slt = compute_lterm_mean(data=u10_east150, time="season", season_calendar="standard")
v10_east150_slt = compute_lterm_mean(data=v10_east150, time="season", season_calendar="standard")


# visualising wind patterns

projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection": projection})

plot_seasonal_mean(variable="Elevation", data_slt=elev_slt, cmap=Greys, units="m", seasons=["JJA", "DJF"], 
                   axes=[ax1,ax2], fig=fig, vmax=3500, vmin=0, domain="Europe", cbar_pos = [0.90, 0.38, 0.02, 0.25], title=True,
                   season_label= ["[A] Alps 100%", "[B]"], levels=22, data_u=u10_slt, data_v=v10_slt,
                   plot_winds_streamline=True, level_ticks=6, plot_winds_pattern=False)

plot_seasonal_mean(variable="Elevation", data_slt=elev_east0_slt, cmap=Greys, units="m", seasons=["JJA", "DJF"], 
                   axes=[ax3,ax4], fig=fig, vmax=3500, vmin=0, domain="Europe", cbar_pos = [0.90, 0.38, 0.02, 0.25], title=True,
                   season_label= ["[C] Alps east 0%", "[D]"], levels=22, data_u=u10_east0_slt, data_v=v10_east0_slt,
                   plot_winds_streamline=True, level_ticks=6, plot_winds_pattern=False)

plot_seasonal_mean(variable="Elevation", data_slt=elev_east150_slt, cmap=Greys, units="m", seasons=["JJA", "DJF"], 
                   axes=[ax5,ax6], fig=fig, vmax=3500, vmin=0, domain="Europe", cbar_pos = [0.90, 0.38, 0.02, 0.25], title=True,
                   season_label= ["[E] Alps east 150%", "[F]"], levels=22, data_u=u10_east150_slt, data_v=v10_east150_slt,
                   plot_winds_streamline=True, level_ticks=6, plot_winds_pattern=False)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas first 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "fig9.png"), format= "png", bbox_inches="tight", dpi=600)









