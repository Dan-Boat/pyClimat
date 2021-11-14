#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 11:34:19 2021

@author: dboateng
"""
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *

# using scipy to perfrom student t-test analysis and spearmanr to perform monotonicity analysis between two datasets

# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"
exp_name_control = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_east_0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_east_150 = "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"
years= "1003_1017"
period = "1m"


# reading dataset
control_data, control_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_control, years=years,
                                                  period=period)
east_0_data, east_0_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_0, years=years,
                                                  period=period)
east_150_data, east_150_wiso = read_ECHAM_processed(main_path=module_output_main_path , exp_name= exp_name_east_150, years=years,
                                                  period=period)

#  extracting variables 
#control
temp2 = extract_var(Dataset=control_data , varname="temp2", units="°C")
prec = extract_var(Dataset= control_data , varname="prec", units="mm/month")
d18op = extract_var(Dataset=control_data , varname="d18op", units="per mil", Dataset_wiso= control_wiso)
u10 = extract_var(Dataset=control_data , varname="u10")
v10 = extract_var(Dataset=control_data , varname="v10")

#east_150
temp2_east_150 = extract_var(Dataset=east_150_data , varname="temp2", units="°C")
prec_east_150 = extract_var(Dataset= east_150_data , varname="prec", units="mm/month")
d18op_east_150 = extract_var(Dataset=east_150_data , varname="d18op", units="per mil", Dataset_wiso= east_150_wiso)
u10_east_150 = extract_var(Dataset=east_150_data , varname="u10")
v10_east_150 = extract_var(Dataset=east_150_data , varname="v10")

#east_0
temp2_east_0 = extract_var(Dataset=east_0_data , varname="temp2", units="°C")
prec_east_0 = extract_var(Dataset= east_0_data , varname="prec", units="mm/month")
d18op_east_0 = extract_var(Dataset=east_0_data , varname="d18op", units="per mil", Dataset_wiso= east_0_wiso)
u10_east_0 = extract_var(Dataset=east_0_data , varname="u10")
v10_east_0 = extract_var(Dataset=east_0_data , varname="v10")


east_150_d18op_slt = compute_lterm_diff(data_control=d18op, data_main=d18op_east_150, time="season", season_calendar="standard")

# minlon = -15
# maxlon = 40
# minlat = 35
# maxlat = 65

# # testing
# d18op_w = d18op.groupby("time.season")["DJF"]
# d18op_east_150_w = d18op_east_150.groupby("time.season")["DJF"]

# #stats_results = correlation(dataA=prec, dataB=d18op, use_pearsonr=True, return_pvalue=True)

# stats_results = student_t_test_btn_datasets(dataA=d18op_w, dataB=d18op_east_150_w, return_pvalue=True, minlat=minlat, minlon=minlon, maxlon=maxlon, maxlat=maxlat,
#                                             max_pvalue=0.2)
projection = ccrs.PlateCarree()
fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize=(20, 15), subplot_kw={"projection":projection})

plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=east_150_d18op_slt , cmap=RdBu, units="‰", seasons=["JJA", "DJF"], 
                   axes=[ax1,ax2], fig=fig, vmax=5, vmin=-5, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.69, 0.02, 0.25], title=True, 
                   season_label= ["[A]       Alps east 150% - Alps 100% [JJA]", "[B]        Alps east 150% - Alps 100% [DJF]"],  compare_data1=d18op,
                   compare_data2=d18op_east_150, max_pvalue=0.2, plot_stats=True)

# plot_seasonal_mean(variable='$\delta^{18}$Op vs SMOW', data_slt=east_150_d18op_slt , cmap=RdBu, units="‰", seasons=["DJF"], 
#                    axes=[ax], fig=fig, vmax=5, vmin=-5, levels=22, domain="Europe", level_ticks=11, cbar_pos = [0.90, 0.69, 0.02, 0.25], title=True, 
#                    season_label= ["[A]       Alps east 150% - Alps 100% [JJA]", "[B]        Alps east 150% - Alps 100% [DJF]"], compare_data1=d18op,
#                    compare_data2=d18op_east_150, max_pvalue=0.2, plot_stats=True)

#ax.contourf(stats_results.lon.values, stats_results.lat.values, stats_results.t_statistic.values, colors="none", hatches=["//"],)


print(stats_results)

# corr_spearmanr = np.zeros((temp2.shape[1], temp2.shape[2]), dtype=float)

# #looping over lat and lon
# for i in range(0, temp2.shape[1]):
#     for j in range(0, temp2.shape[2]):
#        #r,p= stats.spearmanr(temp2.isel(lat=[i], lon=[j]).values.ravel(), d18op.isel(lat=[i], lon=[j]).values.ravel())
#        r,p= stats.ttest_ind(temp2.isel(lat=[i], lon=[j]).values, temp2_east_150.isel(lat=[i], lon=[j]).values)

#        if p >= 0.1:
#            corr_spearmanr[i,j] = np.nan
#        else:
#             corr_spearmanr[i,j] = r
# print(corr_spearmanr)

# # storing in dataset
# stats_dataset = xr.Dataset(coords={"lon": (["lon"], temp2.lon), "lat": (["lat"], temp2.lat)})
# stats_dataset["spearmanr"]= (["lat", "lon"], corr_spearmanr)