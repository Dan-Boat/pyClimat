#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 18:43:39 2021

@author: dboateng
This script displays all the topo input configuarations used for sensitivity experiments of Alps project. 
"""
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

#setting paths 
module_output_main_path = "/home/dboateng/Model_output_pst"

a001 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
a002 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
a003 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
#a004 = "a004_hpc-bw_e5w2.3_t159_PI_Alps_east_50_t159l31.6h"
a005 = "a005_hpc-bw_e5w2.3_t159_PI_Alps_150_t159l31.6h"
a006 = "a006_hpc-bw_e5w2.3_t159_PI_Alps_0_t159l31.6h"
t017 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"
a008 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"
a009 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
a010 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"

filename = "T159_jan_surf.nc"

#reading data

a001_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = a001, filename=filename, read_var=True, varname="OROMEA")
a002_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = a002, filename=filename, read_var=True, varname="OROMEA")
a003_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = a003, filename=filename, read_var=True, varname="OROMEA")
a008_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = a008, filename=filename, read_var=True, varname="OROMEA")
a005_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = a005, filename=filename, read_var=True, varname="OROMEA")
a006_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = a006, filename=filename, read_var=True, varname="OROMEA")
t017_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = t017, filename=filename, read_var=True, varname="OROMEA")
a009_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = a009, filename=filename, read_var=True, varname="OROMEA")
a010_data = read_ECHAM_input(main_path=module_output_main_path , exp_name = a010, filename=filename, read_var=True, varname="OROMEA")


#visualisation
projection = ccrs.PlateCarree()
path_to_store = os.path.join(module_output_main_path, "plots")

fig, ((ax1,ax2),(ax3, ax4), (ax5,ax6)) = plt.subplots(nrows = 3, ncols = 2, figsize=(20, 15), subplot_kw={"projection":projection})

plot_echam_topo(variable="Elevation", data=a001_data, cmap=Greys, units="m", ax=ax1, vmax=3500, vmin=0, levels=12, level_ticks=6,
                domain="Europe", title="AW100E200", cbar=False)

plot_echam_topo(variable="Elevation", data=t017_data, cmap=Greys, units="m", ax=ax5, vmax=3500, vmin=0, levels=12, level_ticks=6, 
                domain="Europe", title="A200", cbar=False)

plot_echam_topo(variable="Elevation", data=a003_data, cmap=Greys, units="m", ax=ax2, vmax=3500, vmin=0, levels=12, level_ticks=6, 
                domain="Europe", title="AW100E0", cbar=False)
plot_echam_topo(variable="Elevation", data=a006_data, cmap=Greys, units="m", ax=ax6, vmax=3500, vmin=0, levels=12, level_ticks=6, 
                domain="Europe", title="A0", cbar_position= [0.90, 0.30, 0.02, 0.40], cbar_orientation="vertical", cbar=True, fig=fig)
plot_echam_topo(variable="Elevation", data=a009_data, cmap=Greys, units="m", ax=ax3, vmax=3500, vmin=0, levels=12, level_ticks=6, 
                domain="Europe", title="AW200E100", cbar=False)
plot_echam_topo(variable="Elevation", data=a010_data, cmap=Greys, units="m", ax=ax4, vmax=3500, vmin=0, levels=12, level_ticks=6, 
                domain="Europe", title = "AW200E0",cbar=False)


fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
plt.tight_layout() 
plt.subplots_adjust(left=0.05, right=0.89, top=0.95, bottom=0.06)
plt.savefig(os.path.join(path_to_store, "figS1_new.svg"), format= "svg", bbox_inches="tight", dpi=600)
plt.show()