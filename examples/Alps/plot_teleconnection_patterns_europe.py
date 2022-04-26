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
exp_name_AW100E100 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_AW100E0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_AW100E200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_AW200E100 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
exp_name_AW200E0 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
exp_name_AW200E200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"

years= "1003_1017"
period = "1m"

#ERA-interim

path_to_msl = "/home/dboateng/Datasets/ERA5/monthly_1950_2021/msl_monthly.nc"

slp_ERA = read_ERA_processed(path=path_to_msl, varname="msl") / 100 #Pa--> hpa
slp_ERA = slp_ERA.rename({"longitude": "lon", "latitude":"lat"})

# reading dataset
AW100E100_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E100, years=years,
                                                  period=period, add_name="slp", read_wiso=False)
AW100E0_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E0, years=years,
                                                  period=period, add_name="slp", read_wiso=False)
AW100E200_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW100E200, years=years,
                                                  period=period, add_name="slp", read_wiso=False)
AW200E100_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW200E100, years=years,
                                                  period=period, add_name="slp", read_wiso=False)
AW200E0_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW200E0, years=years,
                                                  period=period, add_name="slp", read_wiso=False)
AW200E200_data = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_AW200E200, years=years,
                                                  period=period, add_name="slp", read_wiso=False)

# reading sea level presure

slp_AW100E100 = extract_var(Dataset=AW100E100_data , varname="slp", units="hPa")
slp_AW100E0 = extract_var(Dataset= AW100E0_data , varname="slp", units="hPa")
slp_AW100E200 = extract_var(Dataset=AW100E200_data , varname="slp", units="hPa")
slp_AW200E100 = extract_var(Dataset=AW200E100_data , varname="slp", units="hPa")
slp_AW200E0 = extract_var(Dataset=AW200E0_data , varname="slp", units="hPa")
slp_AW200E200 = extract_var(Dataset=AW200E200_data , varname="slp", units="hPa")


#performing EOF analysis (using winter data extracted for North Atlantic region)
eofs_AW100E100, pcs_AW100E100, var_frac_AW100E100 = EOF_analysis(data=slp_AW100E100, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="JJA",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=4,
                    return_pcs=True, npcs=4)
eofs_AW100E0, pcs_AW100E0, var_frac_AW100E0 = EOF_analysis(data=slp_AW100E0, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="JJA",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=4,
                    return_pcs=True, npcs=4)
eofs_AW100E200, pcs_AW100E200, var_frac_AW100E200 = EOF_analysis(data=slp_AW100E200, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="JJA",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=4,
                    return_pcs=True, npcs=4)

eofs_AW200E100, pcs_AW200E100, var_frac_AW200E100 = EOF_analysis(data=slp_AW200E100, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="JJA",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=4,
                    return_pcs=True, npcs=4)

eofs_AW200E0, pcs_AW200E0, var_frac_AW200E0 = EOF_analysis(data=slp_AW200E0, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="JJA",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=4,
                    return_pcs=True, npcs=4)

eofs_AW200E200, pcs_AW200E200, var_frac_AW200E200 = EOF_analysis(data=slp_AW200E200, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="JJA",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=4,
                    return_pcs=True, npcs=4)

eofs_ERA, pcs_ERA, var_frac_ERA = EOF_analysis(data=slp_ERA, maxlon=60, minlon=-80, maxlat=80, minlat=20, season="JJA",
                    neofs=4, pcscaling=1, apply_coslat_weights=True, return_variance=True, neigs=4,
                    return_pcs=True, npcs=4)



projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

#apply font and style 

apply_style(fontsize=22, style=None, linewidth=2)
#plt.rcParams.update({'font.family':'Helvitica'})

#plotting teleconnections for ERA5 in Summer
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols=2, figsize=(15, 8), subplot_kw={"projection": projection})

plot_eofsAsCovariance(variable= "slp", data=(eofs_ERA[0]*-1), mode_var=var_frac_ERA[0], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", use_AlberEqualArea=False,
                      ax=ax1, fig=fig, title="[A]", bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_ERA[1]*-1), mode_var=var_frac_ERA[1], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=False,
                      ax=ax2, fig=fig, title="[B]", bottom_labels=False, left_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_ERA[2]), mode_var=var_frac_ERA[2], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=False,
                      ax=ax3, fig=fig, title="[C]", )

plot_eofsAsCovariance(variable= "slp", data=(eofs_ERA[3]), mode_var=var_frac_ERA[3], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, use_AlberEqualArea=False,
                      ax=ax4, fig=fig, title="[D]", left_labels=False)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.94, top=0.94, bottom=0.10)
plt.savefig(os.path.join(path_to_store, "figS13.svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(path_to_store, "figS13.png"), format= "png", bbox_inches="tight", dpi=300)




# NAO
fig, ((ax1,ax2), (ax3,ax4), (ax5, ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20,16), subplot_kw={"projection":projection})

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E100[0]*-1), mode_var=var_frac_AW100E100[0], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", 
                      ax=ax1, fig=fig, title="[A] CTL", bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E0[0]), mode_var=var_frac_AW100E0[0], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] W1E0", left_labels=False, bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E200[0]*-1), mode_var=var_frac_AW100E200[0], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C] W1E2", bottom_labels=False, )

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E0[0]), mode_var=var_frac_AW200E0[0], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D] W2E0", left_labels=False, bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E100[0]), mode_var=var_frac_AW200E100[0], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax5, fig=fig, title="[E] W2E1",)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E200[0]), mode_var=var_frac_AW200E200[0], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax6, fig=fig, title="[F] W2E2", left_labels=False)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.94, top=0.94, bottom=0.10)
plt.savefig(os.path.join(path_to_store, "fig12.svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(path_to_store, "fig12.png"), format= "png", bbox_inches="tight", dpi=300)


# SCAN

fig, ((ax1,ax2), (ax3,ax4), (ax5, ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20,16), subplot_kw={"projection":projection})

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E100[1]*-1), mode_var=var_frac_AW100E100[1], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", 
                      ax=ax1, fig=fig, title="[A] CTL", bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E0[1]*-1), mode_var=var_frac_AW100E0[1], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] W1E0", left_labels=False, bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E200[1]), mode_var=var_frac_AW100E200[1], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C] W1E2", bottom_labels=False, )

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E0[1]), mode_var=var_frac_AW200E0[1], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D] W2E0", left_labels=False, bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E100[2]*-1), mode_var=var_frac_AW200E100[2], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax5, fig=fig, title="[E] W2E1",)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E200[2]*-1), mode_var=var_frac_AW200E200[2], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax6, fig=fig, title="[F] W2E2", left_labels=False)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.94, top=0.94, bottom=0.10)
plt.savefig(os.path.join(path_to_store, "fig13.svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(path_to_store, "fig13.png"), format= "png", bbox_inches="tight", dpi=300)



# EAWR

fig, ((ax1,ax2), (ax3,ax4), (ax5, ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20,16), subplot_kw={"projection":projection})

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E100[2]), mode_var=var_frac_AW100E100[2], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", 
                      ax=ax1, fig=fig, title="[A] CTL", bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E0[2]), mode_var=var_frac_AW100E0[2], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] W1E0", left_labels=False, bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E200[2]*-1), mode_var=var_frac_AW100E200[2], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C] W1E2", bottom_labels=False, )

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E0[3]*-1), mode_var=var_frac_AW200E0[3], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D] W2E0", left_labels=False, bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E100[3]), mode_var=var_frac_AW200E100[3], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax5, fig=fig, title="[E] W2E1",)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E200[1]), mode_var=var_frac_AW200E200[1], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax6, fig=fig, title="[F] W2E2", left_labels=False)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.94, top=0.94, bottom=0.10)
plt.savefig(os.path.join(path_to_store, "figS14.svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(path_to_store, "figS14.png"), format= "png", bbox_inches="tight", dpi=300)




#EA

fig, ((ax1,ax2), (ax3,ax4), (ax5, ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20,16), subplot_kw={"projection":projection})

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E100[3]*-1), mode_var=var_frac_AW100E100[3], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=True, cbar_position= [0.30, 0.07, 0.30, 0.02], cbar_orientation="horizontal", 
                      ax=ax1, fig=fig, title="[A] CTL", bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E0[3]), mode_var=var_frac_AW100E0[3], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax2, fig=fig, title="[B] W1E0", left_labels=False, bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW100E200[3]*-1), mode_var=var_frac_AW100E200[3], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax3, fig=fig, title="[C] W1E2", bottom_labels=False, )

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E0[2]), mode_var=var_frac_AW200E0[2], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax4, fig=fig, title="[D] W2E0", left_labels=False, bottom_labels=False)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E100[1]*-1), mode_var=var_frac_AW200E100[1], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax5, fig=fig, title="[E] W2E1",)

plot_eofsAsCovariance(variable= "slp", data=(eofs_AW200E200[3]*-1), mode_var=var_frac_AW200E200[3], units="hPa", vmax=5, vmin=-5, cmap=RdBu_r, domain="NH", levels=22,
                      level_ticks=11, cbar=False, ax=ax6, fig=fig, title="[F] W2E2", left_labels=False)

fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.94, top=0.94, bottom=0.10)
plt.savefig(os.path.join(path_to_store, "figS15.svg"), format= "svg", bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(path_to_store, "figS15.png"), format= "png", bbox_inches="tight", dpi=300)

plt.show()