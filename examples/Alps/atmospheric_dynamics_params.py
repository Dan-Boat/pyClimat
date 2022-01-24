#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 13:22:28 2022

@author: dboateng
"""

#importing modules 
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *


# defining module paths
module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name_AW100E100 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_AW100E0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_AW100E200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"

exp_name_AW200E100 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
exp_name_AW200E0 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
#exp_name_AW200E200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"


years= "1003_1017"
period = "1m"

# reading datasets
#omega
AW100E100_omega = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E100, years=years,
                                                  period=period, add_name="omega", read_wiso=False)
AW100E0_omega = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E0, years=years,
                                                  period=period, add_name="omega", read_wiso=False)
AW100E200_omega = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E200, years=years,
                                                  period=period, add_name="omega", read_wiso=False)

# AW200E200_omega = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW200E200, years=years,
#                                                   period=period, add_name="omega", read_wiso=False)
AW200E0_omega = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW200E0, years=years,
                                                  period=period, add_name="omega", read_wiso=False)
AW200E100_omega = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW200E100, years=years,
                                                  period=period, add_name="omega", read_wiso=False)


#pressure variables 
AW100E100_plev = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E100, years=years,
                                                  period=period, add_name="plev", read_wiso=False)
AW100E0_plev = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E0, years=years,
                                                  period=period, add_name="plev", read_wiso=False)
AW100E200_plev = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E200, years=years,
                                                  period=period, add_name="plev", read_wiso=False)

#dataset
AW100E100_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E100, years=years,
                                                  period=period, read_wiso=False)
AW100E0_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E0, years=years,
                                                  period=period, read_wiso=False)
AW100E200_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E200, years=years,
                                                  period=period, read_wiso=False)
#coordinates
maxlon = 21
minlon = -6
maxlat = 48
minlat = 45
#extracting variables 

#vertical velocity
omega = extract_var(Dataset=AW100E100_omega, varname="omega", units="Pa/s", lev_units="hPa")
omega_AW100E0 = extract_var(Dataset=AW100E0_omega, varname="omega", units="Pa/s", lev_units="hPa")
omega_AW100E200 = extract_var(Dataset=AW100E200_omega, varname="omega", units="Pa/s", lev_units="hPa")

#computing seasonal long-term means 
omega_slt = compute_lterm_mean(data=omega, time="season", season_calendar="standard")
omega_AW100E0_slt = compute_lterm_mean(data=omega_AW100E0, time="season", season_calendar="standard")
omega_AW100E200_slt = compute_lterm_mean(data=omega_AW100E200, time="season", season_calendar="standard")

#extract profile 
df_omega = extract_vertical_section(data=omega_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
df_omega_AW100E0 = extract_vertical_section(data=omega_AW100E0_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
df_omega_AW100E200 = extract_vertical_section(data=omega_AW100E200_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")


#elevation
elev = extract_var(Dataset=AW100E100_data , varname="elev", units="m")
elev_AW100E0 = extract_var(Dataset=AW100E0_data , varname="elev", units="m")
elev_AW100E200 = extract_var(Dataset=AW100E200_data , varname="elev", units="m")

#lterm means
elev_slt = compute_lterm_mean(data=elev, time="season", season_calendar="standard")
elev_AW100E0_slt = compute_lterm_mean(data=elev_AW100E0, time="season", season_calendar="standard")
elev_AW100E200_slt = compute_lterm_mean(data=elev_AW100E200, time="season", season_calendar="standard")

# extract profile 
df_elev = extract_profile(data = elev_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", to_pandas=True)
df_elev_AW100E0 = extract_profile(data = elev_AW100E0_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", to_pandas=True)
df_elev_AW100E200 = extract_profile(data = elev_AW100E200_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", to_pandas=True)

#cloud cover 
aclcac = extract_var(Dataset=AW100E100_plev, varname="aclcac", lev_units="hPa")
aclcac_AW100E0 = extract_var(Dataset=AW100E0_plev, varname="aclcac", lev_units="hPa")
aclcac_AW100E200 = extract_var(Dataset=AW100E200_plev, varname="aclcac", lev_units="hPa")

#lterm means
aclcac_slt = compute_lterm_mean(data=aclcac, time="season", season_calendar="standard")
aclcac_AW100E0_slt = compute_lterm_mean(data=aclcac_AW100E0, time="season", season_calendar="standard")
aclcac_AW100E200_slt = compute_lterm_mean(data=aclcac_AW100E200, time="season", season_calendar="standard")

#extract_profile 
df_aclcac = extract_vertical_section(data=aclcac_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
df_aclcac_AW100E0 = extract_vertical_section(data=aclcac_AW100E0_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
df_aclcac_AW100E200 = extract_vertical_section(data=aclcac_AW100E200_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")

#specific humidity
q = extract_var(Dataset=AW100E100_plev, varname="q", lev_units="hPa", units="g/kg")
q_AW100E0 = extract_var(Dataset=AW100E0_plev, varname="q", lev_units="hPa", units="g/kg")
q_AW100E200 = extract_var(Dataset=AW100E200_plev, varname="q", lev_units="hPa", units="g/kg")

#lterm means
q_slt = compute_lterm_mean(data=q, time="season", season_calendar="standard")
q_AW100E0_slt = compute_lterm_mean(data=q_AW100E0, time="season", season_calendar="standard")
q_AW100E200_slt = compute_lterm_mean(data=q_AW100E200, time="season", season_calendar="standard")

#extract_profile 
df_q = extract_vertical_section(data=q_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
df_q_AW100E0 = extract_vertical_section(data=q_AW100E0_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
df_q_AW100E200 = extract_vertical_section(data=q_AW100E200_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")



path_to_store = os.path.join(module_output_main_path, "plots")

#visualization 
plt.rcParams.update({'font.size': 22, "font.weight":"bold"})

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9))= plt.subplots(nrows = 3, ncols = 3, sharey=False, sharex=False, figsize=(25, 20))

plot_vertical_section(variable="Omega", data=df_omega , cmap=YlGnBu, units="Pa/s", season="JJA", vmax=0.08, vmin=-0.08, levels=22,
                      level_ticks=6, plot_colorbar=True, cbar_pos=[0.90, 0.72, 0.02, 0.23], geosp_data=df_elev, dim="lon", ax=ax1, fig=fig, 
                      bottom_labels=False, right_labels=False, left_labels=True, title= "[A]            AW100E100")
plot_vertical_section(variable="Omega", data=df_omega_AW100E0 , cmap=YlGnBu, units="Pa/s", season="JJA", vmax=0.08, vmin=-0.08, levels=22,
                      level_ticks=6, plot_colorbar=False, geosp_data=df_elev_AW100E0 , dim="lon", ax=ax2, fig=fig,
                      bottom_labels=False, right_labels=False, left_labels=False, title="[B]             AW100E0")
plot_vertical_section(variable="Omega", data=df_omega_AW100E200 , cmap=YlGnBu, units="Pa/s", season="JJA", vmax=0.08, vmin=-0.08, levels=22,
                      level_ticks=6, plot_colorbar=False, geosp_data=df_elev_AW100E200, dim="lon", ax=ax3, fig=fig,
                      bottom_labels=False, right_labels=True, left_labels=False, title="[C]             AW100E200")

plot_vertical_section(variable="Clouds", data=df_aclcac , cmap=PuBu, units="-", season="JJA", vmax=0.3, vmin=0, levels=22,
                      level_ticks=6, plot_colorbar=True, cbar_pos=[0.90, 0.39, 0.02, 0.23], geosp_data=df_elev, dim="lon", ax=ax4, fig=fig, 
                      bottom_labels=False, right_labels=False, left_labels=True, title="[D]")
plot_vertical_section(variable="Clouds", data=df_aclcac_AW100E0, cmap=PuBu, units="-", season="JJA", vmax=0.3, vmin=0, levels=22,
                      level_ticks=6, plot_colorbar=False, geosp_data=df_elev_AW100E0 , dim="lon", ax=ax5, fig=fig,
                      bottom_labels=False, right_labels=False, left_labels=False, title="[E]")
plot_vertical_section(variable="Clouds", data=df_aclcac_AW100E200 , cmap=PuBu, units="-", season="JJA", vmax=0.3, vmin=0, levels=22,
                      level_ticks=6, plot_colorbar=False, geosp_data=df_elev_AW100E200, dim="lon", ax=ax6, fig=fig,
                      bottom_labels=False, right_labels=True, left_labels=False, title="[F]")

plot_vertical_section(variable="q", data=df_q , cmap=RdYlBu, units="g/kg", season="JJA", vmax=10, vmin=0, levels=22,
                      level_ticks=6, plot_colorbar=True, cbar_pos=[0.90, 0.04, 0.02, 0.23], geosp_data=df_elev, dim="lon", ax=ax7, fig=fig, 
                      bottom_labels=True, right_labels=False, left_labels=True, title="[G]")
plot_vertical_section(variable="q", data=df_q_AW100E0, cmap=RdYlBu, units="g/kg", season="JJA", vmax=10, vmin=0, levels=22,
                      level_ticks=6, plot_colorbar=False, geosp_data=df_elev_AW100E0 , dim="lon", ax=ax8, fig=fig,
                      bottom_labels=True, right_labels=False, left_labels=False, title="[H]")
plot_vertical_section(variable="q", data=df_q_AW100E200 , cmap=RdYlBu, units="g/kg", season="JJA", vmax=10, vmin=0, levels=22,
                      level_ticks=6, plot_colorbar=False, geosp_data=df_elev_AW100E200, dim="lon", ax=ax9, fig=fig,
                      bottom_labels=True, right_labels=True, left_labels=False, title="[I]")



plt.tight_layout() 
plt.subplots_adjust(left=0.02, right=0.86, top=0.98, bottom=0.03)
plt.savefig(os.path.join(path_to_store, "fig9.svg"), format= "svg", bbox_inches="tight", dpi=600)




