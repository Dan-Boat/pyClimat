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
exp_name_aw100e100 = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_aw100e0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_aw100e200 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
exp_name_aw100e150 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"

# west set up 
exp_name_aw200e100 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
exp_name_aw200e0 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
exp_name_aw200e200 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"
exp_name_aw0e0 = "t015_dkrz-mistral_e5w2.3_PI_Alps0_t159l31.6h"


# for supplementary (same but for annual)
years= "1003_1017"
period = "1m"


# reading dataset
aw100e100_data, aw100e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e100, years=years,
                                                  period=period)
aw100e0_data, aw100e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e0, years=years,
                                                  period=period)
aw100e200_data, aw100e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e200, years=years,
                                                  period=period)
aw100e150_data, aw100e150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw100e150, years=years,
                                                  period=period)

aw200e100_data, aw200e100_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e100, years=years,
                                                  period=period)
aw200e0_data, aw200e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e0, years=years,
                                                  period=period)
aw200e200_data, aw200e200_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw200e200, years=years,
                                                  period=period)

aw0e0_data, aw0e0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_aw0e0, years=years,
                                                  period=period)

#extracting variables and computing long-term means

#aw100e100
d18op_aw100e100 = extract_var(Dataset=aw100e100_data , varname="d18op", units="per mil", Dataset_wiso= aw100e100_wiso)
elev_aw100e100 = extract_var(Dataset=aw100e100_data , varname="elev", units="m")

d18op_aw100e100_slt = compute_lterm_mean(data=d18op_aw100e100, time="season", season_calendar="standard")
elev_aw100e100_slt = compute_lterm_mean(data=elev_aw100e100 , time="season", season_calendar="standard")

#aw100e0
d18op_aw100e0 = extract_var(Dataset=aw100e0_data , varname="d18op", units="per mil", Dataset_wiso= aw100e0_wiso)
elev_aw100e0 = extract_var(Dataset=aw100e0_data , varname="elev", units="m")

d18op_aw100e0_slt = compute_lterm_mean(data=d18op_aw100e0, time="season", season_calendar="standard")
elev_aw100e0_slt = compute_lterm_mean(data=elev_aw100e0 , time="season", season_calendar="standard")

#aw100e200
d18op_aw100e200 = extract_var(Dataset=aw100e200_data , varname="d18op", units="per mil", Dataset_wiso= aw100e200_wiso)
elev_aw100e200 = extract_var(Dataset=aw100e200_data , varname="elev", units="m")

d18op_aw100e200_slt = compute_lterm_mean(data=d18op_aw100e200, time="season", season_calendar="standard")
elev_aw100e200_slt = compute_lterm_mean(data=elev_aw100e200 , time="season", season_calendar="standard")

#aw100e150
d18op_aw100e150 = extract_var(Dataset=aw100e150_data , varname="d18op", units="per mil", Dataset_wiso= aw100e150_wiso)
elev_aw100e150 = extract_var(Dataset=aw100e150_data , varname="elev", units="m")

d18op_aw100e150_slt = compute_lterm_mean(data=d18op_aw100e150, time="season", season_calendar="standard")
elev_aw100e150_slt = compute_lterm_mean(data=elev_aw100e150 , time="season", season_calendar="standard")


#aw200e100
d18op_aw200e100 = extract_var(Dataset=aw200e100_data , varname="d18op", units="per mil", Dataset_wiso= aw200e100_wiso)
elev_aw200e100 = extract_var(Dataset=aw200e100_data , varname="elev", units="m")

d18op_aw200e100_slt = compute_lterm_mean(data=d18op_aw200e100, time="season", season_calendar="standard")
elev_aw200e100_slt = compute_lterm_mean(data=elev_aw200e100 , time="season", season_calendar="standard")

#aw200e0
d18op_aw200e0 = extract_var(Dataset=aw200e0_data , varname="d18op", units="per mil", Dataset_wiso= aw200e0_wiso)
elev_aw200e0 = extract_var(Dataset=aw200e0_data , varname="elev", units="m")

d18op_aw200e0_slt = compute_lterm_mean(data=d18op_aw200e0, time="season", season_calendar="standard")
elev_aw200e0_slt = compute_lterm_mean(data=elev_aw200e0 , time="season", season_calendar="standard")

#aw200e200
d18op_aw200e200 = extract_var(Dataset=aw200e200_data , varname="d18op", units="per mil", Dataset_wiso= aw200e200_wiso)
elev_aw200e200 = extract_var(Dataset=aw200e200_data , varname="elev", units="m")

d18op_aw200e200_slt = compute_lterm_mean(data=d18op_aw200e200, time="season", season_calendar="standard")
elev_aw200e200_slt = compute_lterm_mean(data=elev_aw200e200 , time="season", season_calendar="standard")

#aw0e0
d18op_aw0e0 = extract_var(Dataset=aw0e0_data , varname="d18op", units="per mil", Dataset_wiso= aw0e0_wiso)
elev_aw0e0 = extract_var(Dataset=aw0e0_data , varname="elev", units="m")

d18op_aw0e0_slt = compute_lterm_mean(data=d18op_aw0e0, time="season", season_calendar="standard")
elev_aw0e0_slt = compute_lterm_mean(data=elev_aw0e0 , time="season", season_calendar="standard")


# extracting profiles 
#coordinates 

maxlon_A = 20
minlon_A = -5
maxlat_A = 47
minlat_A = 46

maxlon_B = 13
minlon_B = 10
maxlat_B = 54
minlat_B = 40


#aw100e100
aw100e100_d18op_lon = extract_profile(data = d18op_aw100e100_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw100e100_geosp_lon = extract_profile(data = elev_aw100e100_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw100e100_d18op_lat = extract_profile(data = d18op_aw100e100_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
aw100e100_geosp_lat = extract_profile(data = elev_aw100e100_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#aw100e0
aw100e0_d18op_lon = extract_profile(data = d18op_aw100e0_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw100e0_geosp_lon = extract_profile(data = elev_aw100e0_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw100e0_d18op_lat = extract_profile(data = d18op_aw100e0_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
aw100e0_geosp_lat = extract_profile(data = elev_aw100e0_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#aw100e200
aw100e200_d18op_lon = extract_profile(data = d18op_aw100e200_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw100e200_geosp_lon = extract_profile(data = elev_aw100e200_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw100e200_d18op_lat = extract_profile(data = d18op_aw100e200_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
aw100e200_geosp_lat = extract_profile(data = elev_aw100e200_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#aw100e150
aw100e150_d18op_lon = extract_profile(data = d18op_aw100e150_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw100e150_geosp_lon = extract_profile(data = elev_aw100e150_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw100e150_d18op_lat = extract_profile(data = d18op_aw100e150_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
aw100e150_geosp_lat = extract_profile(data = elev_aw100e150_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#aw200e200
aw200e200_d18op_lon = extract_profile(data = d18op_aw200e200_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw200e200_geosp_lon = extract_profile(data = elev_aw200e200_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw200e200_d18op_lat = extract_profile(data = d18op_aw200e200_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
aw200e200_geosp_lat = extract_profile(data = elev_aw200e200_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#aw200e100
aw200e100_d18op_lon = extract_profile(data = d18op_aw200e100_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw200e100_geosp_lon = extract_profile(data = elev_aw200e100_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw200e100_d18op_lat = extract_profile(data = d18op_aw200e100_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
aw200e100_geosp_lat = extract_profile(data = elev_aw200e100_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

#aw200e0
aw200e0_d18op_lon = extract_profile(data = d18op_aw200e0_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw200e0_geosp_lon = extract_profile(data = elev_aw200e0_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw200e0_d18op_lat = extract_profile(data = d18op_aw200e0_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
aw200e0_geosp_lat = extract_profile(data = elev_aw200e0_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)


#aw0e0
aw0e0_d18op_lon = extract_profile(data = d18op_aw0e0_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw0e0_geosp_lon = extract_profile(data = elev_aw0e0_slt, maxlon=maxlon_A, minlon=minlon_A, maxlat=maxlat_A, minlat=minlat_A, dim="lon", to_pandas=True)
aw0e0_d18op_lat = extract_profile(data = d18op_aw0e0_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)
aw0e0_geosp_lat = extract_profile(data = elev_aw0e0_slt, maxlon=maxlon_B, minlon=minlon_B, maxlat=maxlat_B, minlat=minlat_B, dim="lat", to_pandas=True)

# visualisation
path_to_store = os.path.join(module_output_main_path, "plots")

#apply font and style 

apply_style(fontsize=22, style=None, linewidth=2)

def main_text():

    fig, ((ax1,ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(20, 15),sharey=False, sharex=False)
    
    plot_iso_profiles(df_iso=aw100e100_d18op_lon , df_geosp=aw100e100_geosp_lon , dim="lon", iso_color=black, iso_label="CTL",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, title="[A]", 
                      right_labels =False, bottom_labels=False, shade_color="dimgrey", shade_alpha=0.15)
    
    plot_iso_profiles(df_iso=aw100e0_d18op_lon , df_geosp=aw100e0_geosp_lon , dim="lon", iso_color=blue, iso_label="W1E0",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1,
                      right_labels =False, bottom_labels=False, shade_color="dimgrey", shade_alpha=0.2)
    
    plot_iso_profiles(df_iso=aw100e200_d18op_lon , df_geosp=aw100e200_geosp_lon , dim="lon", iso_color=red, iso_label="W1E2",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1,
                      right_labels =False, bottom_labels=False, shade_color="dimgrey", shade_alpha=0.3)
    
    plot_iso_profiles(df_iso=aw100e150_d18op_lon , df_geosp=aw100e150_geosp_lon , dim="lon", iso_color=green, iso_label="W1E1.5",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1,
                      right_labels =False, bottom_labels=False, shade_color="dimgrey", shade_alpha=0.25)
    
    ax1.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    
    plot_iso_profiles(df_iso=aw100e100_d18op_lon , df_geosp=aw100e100_geosp_lon , dim="lon", iso_color=black, iso_label=None,
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax3, title="[C]", 
                      right_labels =False, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.2)
    
    plot_iso_profiles(df_iso=aw200e100_d18op_lon , df_geosp=aw200e100_geosp_lon , dim="lon", iso_color=golden, iso_label="W2E1",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax3, 
                      right_labels =False, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.3)
    
    plot_iso_profiles(df_iso=aw200e0_d18op_lon , df_geosp=aw200e0_geosp_lon , dim="lon", iso_color=purple, iso_label="W2E0",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax3,  
                      right_labels =False, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.1)
    
    ax3.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    ax3.text(-5, -350, "W", fontsize=20)
    ax3.text(20, -350, "E", fontsize=20)
    
    plot_iso_profiles(df_iso=aw100e100_d18op_lat , df_geosp=aw100e100_geosp_lat , dim="lat", iso_color=black, iso_label=None,
                      season="JJA", xmax=54, xmin=40, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, title="[B]", left_labels=False,
                      bottom_labels=False, shade_color="dimgrey", shade_alpha=0.15)
    
    plot_iso_profiles(df_iso=aw100e0_d18op_lat , df_geosp=aw100e0_geosp_lat , dim="lat", iso_color=blue, iso_label=None,
                      season="JJA", xmax=54, xmin=40, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, left_labels=False,
                      bottom_labels=False, shade_color="dimgrey", shade_alpha=0.2)
    
    plot_iso_profiles(df_iso=aw100e200_d18op_lat , df_geosp=aw100e200_geosp_lat , dim="lat", iso_color=red, iso_label=None,
                      season="JJA", xmax=54, xmin=40, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2,  left_labels=False,
                      bottom_labels=False, shade_color="dimgrey", shade_alpha=0.3)
    
    plot_iso_profiles(df_iso=aw100e150_d18op_lat , df_geosp=aw100e150_geosp_lat , dim="lat", iso_color=green, iso_label=None,
                      season="JJA", xmax=54, xmin=40, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2,  left_labels=False,
                      bottom_labels=False, shade_color="dimgrey", shade_alpha=0.25)
    
    
    ax2.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    
    plot_iso_profiles(df_iso=aw100e100_d18op_lat , df_geosp=aw100e100_geosp_lat , dim="lat", iso_color=black, iso_label=None,
                      season="JJA", xmax=54, xmin=40, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax4, title="[D]", left_labels=False,
                      shade_color="dimgrey", shade_alpha=0.2)
    plot_iso_profiles(df_iso=aw200e100_d18op_lat , df_geosp=aw200e100_geosp_lat , dim="lat", iso_color=golden, iso_label=None,
                      season="JJA", xmax=54, xmin=40, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax4,  left_labels=False,
                      shade_color="dimgrey", shade_alpha=0.3)
    plot_iso_profiles(df_iso=aw200e0_d18op_lat , df_geosp=aw200e0_geosp_lat , dim="lat", iso_color=purple, iso_label=None,
                      season="JJA", xmax=54, xmin=40, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax4, left_labels=False,
                      shade_color="dimgrey", shade_alpha=0.1)
    
    
    ax4.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    ax4.text(40, -350, "S", fontsize=20)
    ax4.text(54, -350, "N", fontsize=20)
    
    fig.legend(frameon=True, fontsize=22, loc="upper right",)
    plt.tight_layout() 
    plt.subplots_adjust(left=0.04, right=0.84, top=0.94, bottom=0.04)
    plt.savefig(os.path.join(path_to_store, "fig5.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "fig5.png"), format= "png", bbox_inches="tight", dpi=300)
    
    
def supplementary():
    
    fig, (ax1,ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(18, 13),sharey=False, sharex=False)

    plot_iso_profiles(df_iso=aw100e100_d18op_lon , df_geosp=aw100e100_geosp_lon , dim="lon", iso_color=black, iso_label="CTL",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1, title="[A]", 
                      right_labels =False, shade_color=black, shade_alpha=0.2)

    plot_iso_profiles(df_iso=aw100e200_d18op_lon , df_geosp=aw100e200_geosp_lon , dim="lon", iso_color=blue, iso_label="W1E2",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1,
                      right_labels =False,  shade_color=blue, shade_alpha=0.15)

    plot_iso_profiles(df_iso=aw200e200_d18op_lon , df_geosp=aw200e200_geosp_lon , dim="lon", iso_color=red, iso_label="W2E2",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax1,
                      right_labels =False, shade_color=red, shade_alpha=0.1)
    
    ax1.grid(visible=True, linestyle="--", linewidth=1.0, color=grey)
    ax1.text(-5, -150, "W", fontsize=20)
    ax1.text(20, -150, "E", fontsize=20)
    
    plot_iso_profiles(df_iso=aw100e100_d18op_lon , df_geosp=aw100e100_geosp_lon , dim="lon", iso_color=black, iso_label=None,
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2, title="[B]", 
                      left_labels =False, shade_color=black, shade_alpha=0.15)

    plot_iso_profiles(df_iso=aw100e0_d18op_lon , df_geosp=aw100e0_geosp_lon , dim="lon", iso_color=green, iso_label="W1E0",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2,
                      left_labels =False, shade_color=green, shade_alpha=0.2)

    plot_iso_profiles(df_iso=aw0e0_d18op_lon , df_geosp=aw0e0_geosp_lon , dim="lon", iso_color=purple, iso_label="W0E0",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3000, isomax=0, isomin=-12, ax=ax2,
                      left_labels =False, shade_color=purple, shade_alpha=0.3)
    
    ax2.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    ax2.text(-5, -150, "W", fontsize=20)
    ax2.text(20, -150, "E", fontsize=20)
    
    
    fig.legend(frameon=True, fontsize=22, loc="upper right",)
    plt.tight_layout() 
    plt.subplots_adjust(left=0.04, right=0.86, top=0.94, bottom=0.04)
    plt.savefig(os.path.join(path_to_store, "figS8.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "figS8.png"), format= "png", bbox_inches="tight", dpi=300)
    

def fig_for_egu_talk():
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(19, 13), sharey=False, sharex=False)
    
    plot_iso_profiles(df_iso=aw100e100_d18op_lon , df_geosp=aw100e100_geosp_lon , dim="lon", iso_color=black, iso_label="CTL",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=0, isomin=-16, ax=ax2, title=None, 
                      right_labels =True, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.15)
    
    plot_iso_profiles(df_iso=aw100e0_d18op_lon , df_geosp=aw100e0_geosp_lon , dim="lon", iso_color=blue, iso_label="W1E0",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=0, isomin=-16, ax=ax2,
                      right_labels =True, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.2)
    
    plot_iso_profiles(df_iso=aw100e200_d18op_lon , df_geosp=aw100e200_geosp_lon , dim="lon", iso_color=red, iso_label="W1E2",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=0, isomin=-16, ax=ax2,
                      right_labels =True, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.3)
    
    plot_iso_profiles(df_iso=aw100e150_d18op_lon , df_geosp=aw100e150_geosp_lon , dim="lon", iso_color=green, iso_label="W1E1.5",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=0, isomin=-16, ax=ax2,
                      right_labels =True, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.25)
    
    ax1.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    ax1.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    ax1.text(-5, -250, "W", fontsize=24)
    ax1.text(20, -250, "E", fontsize=24)
    
    
    plot_iso_profiles(df_iso=aw100e100_d18op_lon , df_geosp=aw100e100_geosp_lon , dim="lon", iso_color=black, iso_label=None,
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax1, title=None, 
                      right_labels =True, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.2)
    
    plot_iso_profiles(df_iso=aw200e100_d18op_lon , df_geosp=aw200e100_geosp_lon , dim="lon", iso_color=golden, iso_label="W2E1",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax1, 
                      right_labels =True, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.3)
    
    plot_iso_profiles(df_iso=aw200e0_d18op_lon , df_geosp=aw200e0_geosp_lon , dim="lon", iso_color=purple, iso_label="W2E0",
                      season="JJA", xmax=20, xmin=-5, ymin=0, ymax=3500, isomax=-2, isomin=-16, ax=ax1,  
                      right_labels =True, bottom_labels=True, shade_color="dimgrey", shade_alpha=0.1)
    
    ax2.grid(visible=True, linestyle="--", linewidth=0.8, color=grey)
    ax2.text(-5, -350, "W", fontsize=24)
    ax2.text(20, -350, "E", fontsize=24)
    
    fig.legend(frameon=True, fontsize=22, loc="upper right",)
    plt.tight_layout() 
    plt.subplots_adjust(left=0.04, right=0.86, top=0.94, bottom=0.04)
    plt.savefig(os.path.join(path_to_store, "egu1.svg"), format= "svg", bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(path_to_store, "egu1.png"), format= "png", bbox_inches="tight", dpi=300)
    
if __name__ == '__main__':
    #supplementary()
    #main_text()
    fig_for_egu_talk()
    
    
    

