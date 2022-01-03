#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 14:09:15 2021

@author: dboateng
"""

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"

exp_name_control = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"

years= "1003_1017"
period = "1m"


# reading dataset
control_omega = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period, add_name="omega", read_wiso=False)

control_plev = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period, add_name="plev", read_wiso=False)

control_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period, read_wiso=False)
control_v = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period, add_name="v", read_wiso=False)
control_u = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period, add_name="u", read_wiso=False)


omega = extract_var(Dataset=control_omega, varname="omega", units="Pa/s", lev_units="hPa")
elev = extract_var(Dataset=control_data , varname="elev", units="m")
aclcac = extract_var(Dataset=control_plev, varname="aclcac", lev_units="hPa")
q = extract_var(Dataset=control_plev, varname="q", lev_units="hPa")
u = extract_var(Dataset=control_u, varname="u", lev_units="hPa")
v = extract_var(Dataset=control_v, varname="v", lev_units="hPa")

#omega_slt = compute_lterm_mean(data=omega, time="season", season_calendar="standard")
elev_slt = compute_lterm_mean(data=elev, time="season", season_calendar="standard")
aclcac_slt = compute_lterm_mean(data=aclcac, time="season", season_calendar="standard")
#q_slt = compute_lterm_mean(data=q, time="season", season_calendar="standard")
#u_slt = compute_lterm_mean(data=u, time="season", season_calendar="standard")
#v_slt = compute_lterm_mean(data=v, time="season", season_calendar="standard")

maxlon = 21
minlon = -6
maxlat = 47
minlat = 46

#data_om = extract_vertical_section(data=omega_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
h_data = extract_profile(data = elev_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", to_pandas=True)
data_cc = extract_vertical_section(data=aclcac_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
#data_q = extract_vertical_section(data=q_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
#data_u = extract_vertical_section(data=u_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")
#data_v = extract_vertical_section(data=v_slt, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, dim="lon", season="JJA")

# plotting 
# import matplotlib
# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 22}

# matplotlib.rc('font', **font)

plt.rcParams.update({'font.size': 22})
fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2,sharey=False, sharex=False, figsize=(20, 8))

plot_vertical_section(variable="Cloud Cover", data=data_cc, cmap=plt.cm.YlGnBu, units="-", season="JJA", vmax=0.25, vmin=0, levels=22, output_name=None, 
                     output_format=None, level_ticks=6, title=None, path_to_store=None, plot_colorbar=True,
                     cbar_pos=None, fig_title=None, season_label=None, geosp_data=h_data, dim="lon", ax=ax1, fig=fig)

plt.tight_layout() 
plt.subplots_adjust(left=0.02, right=0.86, top=0.98, bottom=0.03)


# import matplotlib.pyplot as plt 

# # extracting coordinates from data 
# X = data_u.index.values
# Y = data_u.columns.values 
# U = data_u.values.T
# V = data_v.values.T

# skip = (slice(None, None, 2), slice(None, None, 3))

# x,y = np.meshgrid(X,Y)

# fig,ax = plt.subplots()
# norm = MidpointNormalize(midpoint = 0)
# a = ax.quiver(x[skip],y[skip], U[skip], V[skip], pivot= "mid", scale= 50,)
# #a =ax.contourf(x,y,Z, cmap=plt.cm.YlGnBu, )#levels=21, vmax=0.3, vmin=0,)
# cb = plt.colorbar(a)
# ax.invert_yaxis()

# ax2 = ax.twinx()
# ax2.grid(False)

# df_geosp = h_data["JJA"]
# ax2.fill_between(df_geosp.index, df_geosp/1000,  0, color="black", edgecolor = "black", linestyle="-", linewidth=2.5)
# ax2.set_ylim(0, 16)