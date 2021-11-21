#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 10:48:42 2021

@author: dboateng
This is example script for using Climat in analysing and visualising ECHAM module output (long-term seasonal differece means)
The script contains directories of module outputs and path to save plot
Note: it is structured solely for the personal needs of the author, therefore, it must be adapted advisably.
"""
#Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"

exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east_200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_Alps_150 = "a005_hpc-bw_e5w2.3_t159_PI_Alps_150_t159l31.6h"
exp_name_east_150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"
#exp_name_Alps_0 = "a006_hpc-bw_e5w2.3_t159_PI_Alps_0_t159l31.6h"
exp_name_Alps_0 = "t015_dkrz-mistral_e5w2.3_PI_Alps0_t159l31.6h"
exp_name_east_50 = "a004_hpc-bw_e5w2.3_t159_PI_Alps_east_50_t159l31.6h"
exp_name_Alps_200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"

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
Alps_200_data, Alps_200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_Alps_200, years=years,
                                                  period=period)
east_50_data, east_50_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_50, years=years,
                                                  period=period)
east_200_data, east_200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_200, years=years,
                                                  period=period)

#extracting variables and computing long-term means

#control
d18op = extract_var(Dataset=control_data , varname="d18op", units="per mil", Dataset_wiso= control_wiso)
elev = extract_var(Dataset=control_data , varname="elev", units="m")
d18op_slt = compute_lterm_mean(data=d18op, time="season", season_calendar="standard")
elev_slt = compute_lterm_mean(data=elev, time="season", season_calendar="standard")

d18op_alt = compute_lterm_mean(data=d18op, time="annual",)
elev_alt = compute_lterm_mean(data=elev, time="annual",)

# east_150
d18op_east_150 = extract_var(Dataset=east_150_data , varname="d18op", units="per mil", Dataset_wiso= east_150_wiso)
elev_east_150 = extract_var(Dataset=east_150_data , varname="elev", units="m")

east_150_d18op_slt = compute_lterm_mean(data=d18op_east_150, time="season", season_calendar="standard")
east_150_elev_slt = compute_lterm_mean(data=elev_east_150, time="season", season_calendar="standard")

east_150_d18op_alt = compute_lterm_mean(data=d18op_east_150, time="annual",)
east_150_elev_alt = compute_lterm_mean(data=elev_east_150, time="annual",)

#east_0
d18op_east_0 = extract_var(Dataset=east_0_data , varname="d18op", units="per mil", Dataset_wiso= east_0_wiso)
elev_east_0 = extract_var(Dataset=east_0_data , varname="elev", units="m")
east_0_d18op_slt = compute_lterm_mean(data=d18op_east_0, time="season", season_calendar="standard")
east_0_elev_slt = compute_lterm_mean(data=elev_east_0, time="season", season_calendar="standard")

east_0_d18op_alt = compute_lterm_mean(data=d18op_east_0, time="annual",)
east_0_elev_alt = compute_lterm_mean(data=elev_east_0, time="annual",)

#east_50
d18op_east_50 = extract_var(Dataset=east_50_data , varname="d18op", units="per mil", Dataset_wiso= east_50_wiso)
elev_east_50 = extract_var(Dataset=east_50_data , varname="elev", units="m")
east_50_d18op_slt = compute_lterm_mean(data=d18op_east_50, time="season", season_calendar="standard")
east_50_elev_slt = compute_lterm_mean(data=elev_east_50, time="season", season_calendar="standard")

east_50_d18op_alt = compute_lterm_mean(data=d18op_east_50, time="annual",)
east_50_elev_alt = compute_lterm_mean(data=elev_east_50, time="annual",)

#east_200
d18op_east_200 = extract_var(Dataset=east_200_data , varname="d18op", units="per mil", Dataset_wiso= east_200_wiso)
elev_east_200 = extract_var(Dataset=east_200_data , varname="elev", units="m")
east_200_d18op_slt = compute_lterm_mean(data=d18op_east_200, time="season", season_calendar="standard")
east_200_elev_slt = compute_lterm_mean(data=elev_east_200, time="season", season_calendar="standard")

east_200_d18op_alt = compute_lterm_mean(data=d18op_east_200, time="annual",)
east_200_elev_alt = compute_lterm_mean(data=elev_east_200, time="annual",)

#Alps_150
d18op_Alps_150 = extract_var(Dataset=Alps_150_data , varname="d18op", units="per mil", Dataset_wiso= Alps_150_wiso)
elev_Alps_150 = extract_var(Dataset=Alps_150_data , varname="elev", units="m")
Alps_150_d18op_slt = compute_lterm_mean(data=d18op_Alps_150, time="season", season_calendar="standard")
Alps_150_elev_slt = compute_lterm_mean(data=elev_Alps_150, time="season", season_calendar="standard")

Alps_150_d18op_alt = compute_lterm_mean(data=d18op_Alps_150, time="annual",)
Alps_150_elev_alt = compute_lterm_mean(data=elev_Alps_150, time="annual",)

#Alps_0
d18op_Alps_0 = extract_var(Dataset=Alps_0_data , varname="d18op", units="per mil", Dataset_wiso= Alps_0_wiso)
elev_Alps_0 = extract_var(Dataset=Alps_0_data , varname="elev", units="m")
Alps_0_d18op_slt = compute_lterm_mean(data=d18op_Alps_0, time="season", season_calendar="standard")
Alps_0_elev_slt = compute_lterm_mean(data=elev_Alps_0, time="season", season_calendar="standard")

Alps_0_d18op_alt = compute_lterm_mean(data=d18op_Alps_0, time="annual",)
Alps_0_elev_alt = compute_lterm_mean(data=elev_Alps_0, time="annual",)


#Alps_200
d18op_Alps_200 = extract_var(Dataset=Alps_200_data , varname="d18op", units="per mil", Dataset_wiso= Alps_200_wiso)
elev_Alps_200 = extract_var(Dataset=Alps_200_data , varname="elev", units="m")
Alps_200_d18op_slt = compute_lterm_mean(data=d18op_Alps_200, time="season", season_calendar="standard")
Alps_200_elev_slt = compute_lterm_mean(data=elev_Alps_200, time="season", season_calendar="standard")

Alps_200_d18op_alt = compute_lterm_mean(data=d18op_Alps_200, time="annual",)
Alps_200_elev_alt = compute_lterm_mean(data=elev_Alps_200, time="annual",)

# extracting profiles 
#coordinates 

maxlon_A = 20
minlon_A = -5
maxlat_A = 47
minlat_A = 46

maxlon_B = 15
minlon_B = 13
maxlat_B = 54
minlat_B = 40


#control
control_d18op_lon = extract_profile(data = d18op_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
control_geosp_lon = extract_profile(data = elev_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
control_d18op_lat = extract_profile(data = d18op_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
control_geosp_lat = extract_profile(data = elev_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


control_d18op_lon_a = extract_profile(data = d18op_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
control_geosp_lon_a = extract_profile(data = elev_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
control_d18op_lat_a = extract_profile(data = d18op_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
control_geosp_lat_a = extract_profile(data = elev_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


#east_150
east_150_d18op_lon = extract_profile(data = east_150_d18op_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_150_geosp_lon = extract_profile(data = east_150_elev_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_150_d18op_lat = extract_profile(data = east_150_d18op_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
east_150_geosp_lat = extract_profile(data = east_150_elev_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


east_150_d18op_lon_a = extract_profile(data = east_150_d18op_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_150_geosp_lon_a = extract_profile(data = east_150_elev_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_150_d18op_lat_a = extract_profile(data = east_150_d18op_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
east_150_geosp_lat_a = extract_profile(data = east_150_elev_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#east_0
east_0_d18op_lon = extract_profile(data = east_0_d18op_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_0_geosp_lon = extract_profile(data = east_0_elev_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_0_d18op_lat = extract_profile(data = east_0_d18op_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
east_0_geosp_lat = extract_profile(data = east_0_elev_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


east_0_d18op_lon_a = extract_profile(data = east_0_d18op_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_0_geosp_lon_a = extract_profile(data = east_0_elev_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_0_d18op_lat_a = extract_profile(data = east_0_d18op_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
east_0_geosp_lat_a = extract_profile(data = east_0_elev_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#east_50
east_50_d18op_lon = extract_profile(data = east_50_d18op_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_50_geosp_lon = extract_profile(data = east_50_elev_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_50_d18op_lat = extract_profile(data = east_50_d18op_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
east_50_geosp_lat = extract_profile(data = east_50_elev_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


east_50_d18op_lon_a = extract_profile(data = east_50_d18op_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_50_geosp_lon_a = extract_profile(data = east_50_elev_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_50_d18op_lat_a = extract_profile(data = east_50_d18op_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
east_50_geosp_lat_a = extract_profile(data = east_50_elev_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#east_200
east_200_d18op_lon = extract_profile(data = east_200_d18op_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_200_geosp_lon = extract_profile(data = east_200_elev_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_200_d18op_lat = extract_profile(data = east_200_d18op_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
east_200_geosp_lat = extract_profile(data = east_200_elev_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


east_200_d18op_lon_a = extract_profile(data = east_200_d18op_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_200_geosp_lon_a = extract_profile(data = east_200_elev_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
east_200_d18op_lat_a = extract_profile(data = east_200_d18op_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
east_200_geosp_lat_a = extract_profile(data = east_200_elev_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


#Alps_150
Alps_150_d18op_lon = extract_profile(data = Alps_150_d18op_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_150_geosp_lon = extract_profile(data = Alps_150_elev_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_150_d18op_lat = extract_profile(data = Alps_150_d18op_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
Alps_150_geosp_lat = extract_profile(data = Alps_150_elev_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


Alps_150_d18op_lon_a = extract_profile(data = Alps_150_d18op_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_150_geosp_lon_a = extract_profile(data = Alps_150_elev_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_150_d18op_lat_a = extract_profile(data = Alps_150_d18op_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
Alps_150_geosp_lat_a = extract_profile(data = Alps_150_elev_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#Alps_0
Alps_0_d18op_lon = extract_profile(data = Alps_0_d18op_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_0_geosp_lon = extract_profile(data = Alps_0_elev_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_0_d18op_lat = extract_profile(data = Alps_0_d18op_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
Alps_0_geosp_lat = extract_profile(data = Alps_0_elev_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


Alps_0_d18op_lon_a = extract_profile(data = Alps_0_d18op_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_0_geosp_lon_a = extract_profile(data = Alps_0_elev_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_0_d18op_lat_a = extract_profile(data = Alps_0_d18op_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
Alps_0_geosp_lat_a = extract_profile(data = Alps_0_elev_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


#Alps_200
Alps_200_d18op_lon = extract_profile(data = Alps_200_d18op_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_200_geosp_lon = extract_profile(data = Alps_200_elev_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_200_d18op_lat = extract_profile(data = Alps_200_d18op_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
Alps_200_geosp_lat = extract_profile(data = Alps_200_elev_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


Alps_200_d18op_lon_a = extract_profile(data = Alps_200_d18op_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_200_geosp_lon_a = extract_profile(data = Alps_200_elev_alt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
Alps_200_d18op_lat_a = extract_profile(data = Alps_200_d18op_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
Alps_200_geosp_lat_a = extract_profile(data = Alps_200_elev_alt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

# visualisation

path_to_store = os.path.join(module_output_main_path, "plots")
#fig, ((ax1,ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(20, 13),sharey=False, sharex=False)

apply_style_2(fontsize=22)

fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2,sharey=False, sharex=False, figsize=(20, 8))

plot_iso_profiles(df_iso=control_d18op_lon , df_geosp=control_geosp_lon , dim="lon", iso_color=black, iso_label="Alps 100%",
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, title="[A]", 
                  right_labels =False)
plot_iso_profiles(df_iso=east_150_d18op_lon , df_geosp=east_150_geosp_lon , dim="lon", iso_color=blue, iso_label="Alps east 150%",
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, 
                  right_labels =False)
plot_iso_profiles(df_iso=east_0_d18op_lon , df_geosp=east_0_geosp_lon , dim="lon", iso_color=red, iso_label="Alps east 0%",
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, 
                  right_labels =False)
plot_iso_profiles(df_iso=east_200_d18op_lon , df_geosp=east_200_geosp_lon , dim="lon", iso_color=green, iso_label="Alps east 200%",
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, right_labels =False)


plot_iso_profiles(df_iso=control_d18op_lat , df_geosp=control_geosp_lat , dim="lat", iso_color=black, iso_label=None,
                  season="JJA", xmax=54, xmin=40, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, title="[B]", left_labels=False,
                  )
plot_iso_profiles(df_iso=east_150_d18op_lat , df_geosp=east_150_geosp_lat , dim="lat", iso_color=blue, iso_label=None,
                  season="JJA", xmax=54, xmin=40, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, left_labels=False,
                  )
plot_iso_profiles(df_iso=east_0_d18op_lat , df_geosp=east_0_geosp_lat , dim="lat", iso_color=red, iso_label=None,
                  season="JJA", xmax=54, xmin=40, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, left_labels=False,
                  )
plot_iso_profiles(df_iso=east_200_d18op_lat , df_geosp=east_200_geosp_lat , dim="lat", iso_color=green, iso_label=None,
                  season="JJA", xmax=54, xmin=40, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, left_labels=False,
                 )
# xmax = 20
# xmin = -5
# ymax = 3000
# ymin = 0
# isomax = -3
# isomin = -18

# plot_iso_profiles(df_iso=control_d18op_lon_a, df_geosp= control_geosp_lon_a, dim="lon", iso_color=black, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax3,
#                   right_labels=False)
# plot_iso_profiles(df_iso=east_150_d18op_lon_a, df_geosp= east_150_geosp_lon_a, dim="lon", iso_color=blue, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax3,
#                   right_labels=False)
# plot_iso_profiles(df_iso=east_0_d18op_lon_a, df_geosp= east_0_geosp_lon_a, dim="lon", iso_color=red, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax3,
#                   right_labels=False)
# plot_iso_profiles(df_iso=east_200_d18op_lon_a, df_geosp= east_200_geosp_lon_a, dim="lon", iso_color=green, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax3,
#                   right_labels=False)

# xmax = 54
# xmin = 40

# plot_iso_profiles(df_iso=control_d18op_lat_a, df_geosp= control_geosp_lat_a, dim="lat", iso_color=black, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax4, left_labels =False)
# plot_iso_profiles(df_iso=east_150_d18op_lat_a, df_geosp= east_150_geosp_lat_a, dim="lat", iso_color=blue, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax4, left_labels=False)
# plot_iso_profiles(df_iso=east_0_d18op_lat_a, df_geosp= east_0_geosp_lat_a, dim="lat", iso_color=red, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax4, left_labels=False)
# plot_iso_profiles(df_iso=east_200_d18op_lat_a, df_geosp= east_200_geosp_lat_a, dim="lat", iso_color=green, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax4, left_labels=False)

fig.legend(frameon=True, fontsize=18, loc="upper right",)
plt.tight_layout() 
plt.subplots_adjust(left=0.02, right=0.85, top=0.98, bottom=0.03)
plt.savefig(os.path.join(path_to_store, "fig5.svg"), format= "svg", bbox_inches="tight", dpi=600)

#fig, ((ax1,ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2,sharey=False, sharex=False)

fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2,sharey=False, sharex=False, figsize=(20,8))

apply_style_2()
plot_iso_profiles(df_iso=control_d18op_lon , df_geosp=control_geosp_lon , dim="lon", iso_color=black, iso_label="Alps 100%",
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, title="[A]",
                  right_labels=False)
plot_iso_profiles(df_iso=east_200_d18op_lon , df_geosp=east_200_geosp_lon , dim="lon", iso_color=blue, iso_label="Alps east 200%",
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, 
                  right_labels=False)
plot_iso_profiles(df_iso=Alps_200_d18op_lon , df_geosp=Alps_200_geosp_lon, dim="lon", iso_color=red, iso_label="Alps 200%",
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1,
                  right_labels=False)


plot_iso_profiles(df_iso=control_d18op_lon , df_geosp=control_geosp_lon , dim="lon", iso_color=black, iso_label=None,
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, title="[B]", left_labels=False,
                 )
plot_iso_profiles(df_iso=east_0_d18op_lon , df_geosp=east_0_geosp_lon , dim="lon", iso_color=green, iso_label="Alps east 0%", 
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, left_labels=False,
                  )
plot_iso_profiles(df_iso=Alps_0_d18op_lon , df_geosp=Alps_0_geosp_lon, dim="lon", iso_color=purple, iso_label="Alps 0%",
                  season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, left_labels=False,)


#fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2,sharey=False, sharex=False, figsize=(20,8))
# xmax = 20
# xmin = -5
# ymax = 3000
# ymin = 0
# isomax = -3
# isomin = -18

# plot_iso_profiles(df_iso=control_d18op_lon_a, df_geosp= control_geosp_lon_a, dim="lon", iso_color=black, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax3, title="[C]",
#                   right_labels=False)
# plot_iso_profiles(df_iso=east_200_d18op_lon_a, df_geosp= east_200_geosp_lon_a, dim="lon", iso_color=blue, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax3,
#                   right_labels=False)
# plot_iso_profiles(df_iso=Alps_200_d18op_lon_a, df_geosp= Alps_200_geosp_lon_a, dim="lon", iso_color=red, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax3,
#                   right_labels=False)


# plot_iso_profiles(df_iso=control_d18op_lon_a, df_geosp= control_geosp_lon_a, dim="lon", iso_color=black, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax4, left_labels =False, title="[D]")
# plot_iso_profiles(df_iso=east_0_d18op_lon_a, df_geosp= east_0_geosp_lon_a, dim="lon", iso_color=green, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax4, left_labels=False)
# plot_iso_profiles(df_iso=Alps_0_d18op_lon_a, df_geosp= Alps_0_geosp_lon_a, dim="lon", iso_color=purple, iso_label=None,
#                   xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, isomax=isomax, isomin=isomin, ax=ax4, left_labels=False)


fig.legend(frameon=True, fontsize=18, loc="upper right",)
plt.tight_layout() 
plt.subplots_adjust(left=0.02, right=0.85, top=0.98, bottom=0.03)
plt.savefig(os.path.join(path_to_store, "fig6.svg"), format= "svg", bbox_inches="tight", dpi=600)