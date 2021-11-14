#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 14:00:23 2021

@author: dboateng

This script uses the EOF analyis implemented in the Climat_analysis and its visualisation in Climat_plots to extract the three most dominant EOFs
from slp or geopoth model ouptuts and reanalyis

!!! EOF sometimes forgot the axis of rotation, therefore knowledge of the dipole axis is required to validate the results (multiply by -ve to reverse the pattern)
"""

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east_150 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_east_50 = "a004_hpc-bw_e5w2.3_t159_PI_Alps_east_50_t159l31.6h"
years= "1003_1017"
period = "1m"

#ERA-interim

path_to_msl = "/home/dboateng/Datasets/Reanalysis/era_interim/msl_monthly.nc"

slp_ERA = read_ERA_processed(path=path_to_msl, varname="msl") / 100 #Pa--> hpa
slp_ERA = slp_ERA.rename({"longitude": "lon", "latitude":"lat"})

# reading dataset
control_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period, add_name="slp", read_wiso=False)
east_0_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_0, years=years,
                                                  period=period, add_name="slp", read_wiso=False)
east_150_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_150, years=years,
                                                  period=period, add_name="slp", read_wiso=False)

# reading sea level presure

slp = extract_var(Dataset=control_data , varname="slp", units="hPa")
slp_east_0 = extract_var(Dataset= east_0_data , varname="slp", units="hPa")
slp_east_150 = extract_var(Dataset=east_150_data , varname="slp", units="hPa")

# performing EOF analysis (using winter data extracted for North Atlantic region)

eofs, pcs, var_frac = EOF_analysis(data=slp, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="DJF",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=3,
                    return_pcs=True, npcs=4)
eofs_east_0, pcs_east_0, var_frac_east_0 = EOF_analysis(data=slp_east_0, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="DJF",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=3,
                    return_pcs=True, npcs=4)
eofs_east_150, pcs_east_150, var_frac_east_150 = EOF_analysis(data=slp_east_150, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="DJF",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=3,
                    return_pcs=True, npcs=4)

eofs_ERA, pcs_east_ERA, var_frac_ERA = EOF_analysis(data=slp_ERA, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="DJF",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=3,
                    return_pcs=True, npcs=4)



projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

#projection = ccrs.AlbersEqualArea(central_latitude=40, central_longitude=-35, standard_parallels=(0, 80))

# NAO
fig, ((ax1,ax2), (ax3,ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(18,10), subplot_kw={"projection":projection})

plot_eofsAsCovariance(variable= "slp", data=(eofs_ERA[0]*1), mode_var=var_frac_ERA[0], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal",use_AlberEqualArea=False,
                      ax=ax1, fig=fig, title="[A] ERA-Interim")
plot_eofsAsCovariance(variable= "slp", data=(eofs[0]*-1), mode_var=var_frac[0], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] Alps 100")
plot_eofsAsCovariance(variable= "mslp", data=(eofs_east_0[0]*-1), mode_var=var_frac_east_0[0], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C] Alps east 0")
plot_eofsAsCovariance(variable= "mslp", data=(eofs_east_150[0]*1), mode_var=var_frac_east_150[0], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D] Alps east 150")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.07)
plt.savefig(os.path.join(path_to_store, "figS10.png"), format= "png", bbox_inches="tight", dpi=600)


# EA
fig, ((ax1,ax2), (ax3,ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(18,10), subplot_kw={"projection":projection})

plot_eofsAsCovariance(variable= "slp", data=(eofs_ERA[1]*-1), mode_var=var_frac_ERA[1], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal",use_AlberEqualArea=False,
                      ax=ax1, fig=fig, title="[A] ERA-Interim")
plot_eofsAsCovariance(variable= "slp", data=(eofs[1]*-1), mode_var=var_frac[1], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] Alps 100")
plot_eofsAsCovariance(variable= "slp", data=(eofs_east_0[1]*-1), mode_var=var_frac_east_0[1], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C] Alps east 0")
plot_eofsAsCovariance(variable= "slp", data=(eofs_east_150[1]*-1), mode_var=var_frac_east_150[1], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D] Alps east 150")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.07)
plt.savefig(os.path.join(path_to_store, "fig10.png"), format= "png", bbox_inches="tight", dpi=600)


# SCAN

fig, ((ax1,ax2), (ax3,ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(18,10), subplot_kw={"projection":projection})

plot_eofsAsCovariance(variable= "slp", data=(eofs_ERA[2]*-1), mode_var=var_frac_ERA[2], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal",use_AlberEqualArea=False,
                      ax=ax1, fig=fig, title="[A] ERA-Interim")
plot_eofsAsCovariance(variable= "slp", data=(eofs[2]*-1), mode_var=var_frac[2], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] Alps 100")
plot_eofsAsCovariance(variable= "slp", data=(eofs_east_0[2]*-1), mode_var=var_frac_east_0[2], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C] Alps east 0")
plot_eofsAsCovariance(variable= "slp", data=(eofs_east_150[2]*-1), mode_var=var_frac_east_150[2], units="hPa", vmax=10, vmin=-10, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D] Alps east 150")

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.94, bottom=0.07)
plt.savefig(os.path.join(path_to_store, "figS11.png"), format= "png", bbox_inches="tight", dpi=600)



