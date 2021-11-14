#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 10:18:26 2021

@author: dboateng
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"

exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east_150 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_Alps_150 = "a005_hpc-bw_e5w2.3_t159_PI_Alps_150_t159l31.6h"
exp_name_Alps_0 = "a006_hpc-bw_e5w2.3_t159_PI_Alps_0_t159l31.6h"
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
Alps_150_data, Alps_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_150, years=years,
                                                  period=period)
Alps_0_data, Alps_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_0, years=years,
                                                  period=period)
east_50_data, east_50_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_50, years=years,
                                                  period=period)

#extracting variables and computing long-term means

# control
d18op = extract_var(Dataset=control_data , varname="d18op", units="per mil", Dataset_wiso= control_wiso)
elev = extract_var(Dataset=control_data , varname="elev", units="m")
d18op_alt = compute_lterm_mean(data=d18op, time="annual")
elev_alt = compute_lterm_mean(data=elev, time="annual")

d18op_east_150 = extract_var(Dataset=east_150_data , varname="d18op", units="per mil", Dataset_wiso= east_150_wiso)
elev_east_150 = extract_var(Dataset=east_150_data , varname="elev", units="m")
east_150_d18op_alt = compute_lterm_mean(data=d18op_east_150, time="annual")
east_150_elev_alt = compute_lterm_mean(data=elev_east_150, time="annual")

d18op_east_0 = extract_var(Dataset=east_0_data , varname="d18op", units="per mil", Dataset_wiso= east_0_wiso)
elev_east_0 = extract_var(Dataset=east_0_data , varname="elev", units="m")
east_0_d18op_alt = compute_lterm_mean(data=d18op_east_0, time="annual")
east_0_elev_alt = compute_lterm_mean(data=elev_east_0, time="annual")

d18op_east_50 = extract_var(Dataset=east_50_data , varname="d18op", units="per mil", Dataset_wiso= east_50_wiso)
elev_east_50 = extract_var(Dataset=east_50_data , varname="elev", units="m")
east_50_d18op_alt = compute_lterm_mean(data=d18op_east_50, time="annual")
east_50_elev_alt = compute_lterm_mean(data=elev_east_50, time="annual")

d18op_Alps_150 = extract_var(Dataset=Alps_150_data , varname="d18op", units="per mil", Dataset_wiso= Alps_150_wiso)
elev_Alps_150 = extract_var(Dataset=Alps_150_data , varname="elev", units="m")
Alps_150_d18op_alt = compute_lterm_mean(data=d18op_Alps_150, time="annual")
Alps_150_elev_alt = compute_lterm_mean(data=elev_Alps_150, time="annual")

d18op_Alps_0 = extract_var(Dataset=Alps_0_data , varname="d18op", units="per mil", Dataset_wiso= Alps_0_wiso)
elev_Alps_0 = extract_var(Dataset=Alps_0_data , varname="elev", units="m")
Alps_0_d18op_alt = compute_lterm_mean(data=d18op_Alps_0, time="annual")
Alps_0_elev_alt = compute_lterm_mean(data=elev_Alps_0, time="annual")



control_d18op_lon_alt = extract_profile(data = d18op_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
control_geosp_lon_alt = extract_profile(data = elev_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
control_d18op_lat_alt = extract_profile(data = d18op_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")
control_geosp_lat_alt = extract_profile(data = elev_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")

east_150_d18op_lon_alt = extract_profile(data = east_150_d18op_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
east_150_geosp_lon_alt = extract_profile(data = east_150_elev_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
east_150_d18op_lat_alt = extract_profile(data = east_150_d18op_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")
east_150_geosp_lat_alt = extract_profile(data = east_150_elev_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")

east_0_d18op_lon_alt = extract_profile(data = east_0_d18op_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
east_0_geosp_lon_alt = extract_profile(data = east_0_elev_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
east_0_d18op_lat_alt = extract_profile(data = east_0_d18op_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")
east_0_geosp_lat_alt = extract_profile(data = east_0_elev_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")

east_50_d18op_lon_alt = extract_profile(data = east_50_d18op_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
east_50_geosp_lon_alt = extract_profile(data = east_50_elev_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
east_50_d18op_lat_alt = extract_profile(data = east_50_d18op_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")
east_50_geosp_lat_alt = extract_profile(data = east_50_elev_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")

Alps_150_d18op_lon_alt = extract_profile(data = Alps_150_d18op_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
Alps_150_geosp_lon_alt = extract_profile(data = Alps_150_elev_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
Alps_150_d18op_lat_alt = extract_profile(data = Alps_150_d18op_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")
Alps_150_geosp_lat_alt = extract_profile(data = Alps_150_elev_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")

Alps_0_d18op_lon_alt = extract_profile(data = Alps_0_d18op_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
Alps_0_geosp_lon_alt = extract_profile(data = Alps_0_elev_alt, maxlon=20, minlon=-5, maxlat=47, minlat=46, dim="lon", to_pandas="yes")
Alps_0_d18op_lat_alt = extract_profile(data = Alps_0_d18op_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")
Alps_0_geosp_lat_alt = extract_profile(data = Alps_0_elev_alt, maxlon=14, minlon=12, maxlat=54, minlat=40, dim="lat", to_pandas="yes")

# visualisation

path_to_store = os.path.join(module_output_main_path, "plots")
fig, (ax1,ax2) = plt.subplots(nrows = 2, ncols = 1, figsize=(17, 13),sharey=False, sharex=False)
xmax = 20
xmin = -5
ymax = 3000
ymin = 0
isomax = -2
isomin = -18

#apply_style_2()
plot_iso_profiles(df_iso=control_d18op_lon_alt , df_geosp=control_geosp_lon_alt , dim="lon", iso_color=black, iso_label="Alps 100%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax1, title="[A]")
plot_iso_profiles(df_iso=east_150_d18op_lon_alt , df_geosp=east_150_geosp_lon_alt , dim="lon", iso_color=blue, iso_label="Alps east 150%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax1)
plot_iso_profiles(df_iso=east_0_d18op_lon_alt , df_geosp=east_0_geosp_lon_alt , dim="lon", iso_color=red, iso_label="Alps east 0%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax1)
plot_iso_profiles(df_iso=east_50_d18op_lon_alt , df_geosp=east_50_geosp_lon_alt , dim="lon", iso_color=green, iso_label="Alps east 50%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax1)

xmax = 54
xmin = 40
ymax = 3000
ymin = 0
isomax = -2
isomin = -18

plot_iso_profiles(df_iso=control_d18op_lat_alt , df_geosp=control_geosp_lat_alt , dim="lat", iso_color=black, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax2, title="[B]")
plot_iso_profiles(df_iso=east_150_d18op_lat_alt , df_geosp=east_150_geosp_lat_alt , dim="lat", iso_color=blue, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax2)
plot_iso_profiles(df_iso=east_0_d18op_lat_alt , df_geosp=east_0_geosp_lat_alt , dim="lat", iso_color=red, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax2)
plot_iso_profiles(df_iso=east_50_d18op_lat_alt , df_geosp=east_50_geosp_lat_alt , dim="lat", iso_color=green, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax2)

fig.legend(frameon=True, fontsize=18, loc="upper right",)
plt.tight_layout() 
plt.subplots_adjust(left=0.06, right=0.80, top=0.98, bottom=0.05)
plt.savefig(os.path.join(path_to_store, "figS6.svg"), format= "svg", bbox_inches="tight", dpi=600)

fig, ((ax1,ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(20, 13),sharey=False, sharex=False)

xmax = 20
xmin = -5
ymax = 3000
ymin = 0
isomax = -2
isomin = -18

plot_iso_profiles(df_iso=control_d18op_lon_alt , df_geosp=control_geosp_lon_alt , dim="lon", iso_color=black, iso_label="Alps 100%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax1, title="[A]")
plot_iso_profiles(df_iso=Alps_150_d18op_lon_alt , df_geosp=Alps_150_geosp_lon_alt , dim="lon", iso_color=red, iso_label="Alps 150%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax1)
plot_iso_profiles(df_iso=east_150_d18op_lon_alt , df_geosp=east_150_geosp_lon_alt, dim="lon", iso_color=blue, iso_label="Alps east 150%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax1)

plot_iso_profiles(df_iso=control_d18op_lon_alt , df_geosp=control_geosp_lon_alt , dim="lon", iso_color=black, iso_label="Alps 100%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax2, title="[B]")
plot_iso_profiles(df_iso=Alps_0_d18op_lon_alt , df_geosp=Alps_0_geosp_lon_alt , dim="lon", iso_color=purple, iso_label="Alps 0%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax2)
plot_iso_profiles(df_iso=east_0_d18op_lon_alt , df_geosp=east_0_geosp_lon_alt, dim="lon", iso_color=green, iso_label="Alps east 0%", 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax2)
xmax = 55
xmin = 40
ymax = 3000
ymin = 0
isomax = -2
isomin = -18

plot_iso_profiles(df_iso=control_d18op_lat_alt , df_geosp=control_geosp_lat_alt , dim="lat", iso_color=black, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax3, title="[C]")
plot_iso_profiles(df_iso=Alps_150_d18op_lat_alt , df_geosp=Alps_150_geosp_lat_alt , dim="lat", iso_color=red, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax3)
plot_iso_profiles(df_iso=east_150_d18op_lat_alt , df_geosp=east_150_geosp_lat_alt, dim="lat", iso_color=blue, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax3)

plot_iso_profiles(df_iso=control_d18op_lat_alt , df_geosp=control_geosp_lat_alt , dim="lat", iso_color=black, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax4, title="[D]")
plot_iso_profiles(df_iso=Alps_0_d18op_lat_alt , df_geosp=Alps_0_geosp_lat_alt , dim="lat", iso_color=purple, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax4)
plot_iso_profiles(df_iso=east_0_d18op_lat_alt , df_geosp=east_0_geosp_lat_alt, dim="lat", iso_color=green, iso_label=None, 
                  xmax=xmax, xmin=xmin, ymin=ymin, ymax=ymax, isomax=isomax, isomin=isomin, ax=ax4)

fig.legend(frameon=True, fontsize=18, loc="upper right",)
plt.tight_layout() 
plt.subplots_adjust(left=0.15, right=0.95, top=0.98, bottom=0.05)
plt.savefig(os.path.join(path_to_store, "figS7.png"), format= "png", bbox_inches="tight", dpi=600)

