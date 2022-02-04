#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 15:52:38 2022

@author: dboateng
"""

# Importing ClimatPackages (this will later be compile into a package wheere __init__.py with import all the functions)
import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit")

from Package import *




# Path to experiments
module_output_main_path = "/home/dboateng/Model_output_pst"
dir_ERAtp = "/home/dboateng/Datasets/era_interim/tp_monthly.nc"
dir_ERAt2m = "/home/dboateng/Datasets/era_interim/t2m_monthly.nc"


exp_name_t040 = "t040_dkrz-mistral_e5w2.3_MIOCS_450_ctl_t159l31.6h"
exp_name_t041 = "t041_dkrz-mistral_e5w2.3_MIOCS_278_ctl_t159l31.6h"
exp_name_t018 = "t018_hpc-bw_e5w2.3_PI_Alps100_t159l31.6h"

years = "1003_1009"
period = "1m"

t040_data, t040_wiso = read_ECHAM_processed(main_path=module_output_main_path, exp_name=exp_name_t040, years=years,
                                            period=period)
t041_data, t041_wiso = read_ECHAM_processed(main_path=module_output_main_path, exp_name=exp_name_t041, years=years,
                                            period=period)
years = "1003_1017"

t018_data, t018_wiso = read_ECHAM_processed(main_path=module_output_main_path, exp_name=exp_name_t018, years=years,
                                            period=period)

tp = read_ERA_processed(path=dir_ERAtp, varname="tp")
t2m = read_ERA_processed(path=dir_ERAt2m, varname="t2m")

# extracting variables
d18op_t040 = extract_var(Dataset=t040_data, varname="d18op",
                         units="per mil", Dataset_wiso=t040_wiso)

d18ov_t040 = extract_var(Dataset=t040_data, varname="d18ov",
                         units="per mil", Dataset_wiso=t040_wiso, lev=30)
temp_t040 = extract_var(Dataset=t040_data, varname="temp2", units="°C")
prec_t040 = extract_var(Dataset=t040_data, varname="prec", units="mm/month")
relhum_t040 = extract_var(Dataset=t040_data, varname="relhum", units="%")
evap_t040 = extract_var(Dataset=t040_data, varname="evap", units="mm/month")


d18op_t041 = extract_var(Dataset=t041_data, varname="d18op",
                         units="per mil", Dataset_wiso=t041_wiso)
d18ov_t041 = extract_var(Dataset=t041_data, varname="d18ov",
                         units="per mil", Dataset_wiso=t041_wiso, lev=30)
temp_t041 = extract_var(Dataset=t041_data, varname="temp2", units="°C")
prec_t041 = extract_var(Dataset=t041_data, varname="prec", units="mm/month")
relhum_t041 = extract_var(Dataset=t041_data, varname="relhum", units="%")
evap_t041 = extract_var(Dataset=t041_data, varname="evap", units="mm/month")


d18op_t018 = extract_var(Dataset=t018_data, varname="d18op",
                         units="per mil", Dataset_wiso=t018_wiso)
d18ov_t018 = extract_var(Dataset=t018_data, varname="d18ov",
                         units="per mil", Dataset_wiso=t018_wiso, lev=30)
temp_t018 = extract_var(Dataset=t018_data, varname="temp2", units="°C")
prec_t018 = extract_var(Dataset=t018_data, varname="prec", units="mm/month")
relhum_t018 = extract_var(Dataset=t018_data, varname="relhum", units="%")
evap_t018 = extract_var(Dataset=t018_data, varname="evap", units="mm/month")




max_lon, max_lat, min_lon, min_lat = 25, 55, -2, 40

def weighted_mean(data):
    weights = np.cos(np.deg2rad(data.lat))
    weights.name = "weights"
    
    data_weighted = data.weighted(weights)
    
    data_mean = data_weighted.mean(dim=("lat", "lon"), skipna=True)
    return data_mean 

def weighted_std(data):
    weights = np.cos(np.deg2rad(data.lat))
    weights.name = "weights"
    
    data_weighted = data.weighted(weights)
    
    data_std = data_weighted.std(dim=("lat", "lon"), skipna=True)
    return data_std



#central europe 
def monthly_means_central_europe(d18op_data, d18ov_data, temp_data, prec_data, evap_data, relhum_data, Dataset):
    max_lon, max_lat, min_lon, min_lat = 25, 55, -2, 40
    
    d18op_europe = extract_transect(data=d18op_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=None, Dataset=Dataset)
    temp_europe = extract_transect(data=temp_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=None, Dataset=Dataset)
    prec_europe = extract_transect(data=prec_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=None, Dataset=Dataset)
    relhum_europe = extract_transect(data=relhum_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=None, Dataset=Dataset)
    evap_europe = extract_transect(data=evap_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=None, Dataset=Dataset)
    d18ov_europe = extract_transect(data=d18ov_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=None, Dataset=Dataset)
    
    return d18op_europe, d18ov_europe, temp_europe, prec_europe, evap_europe, relhum_europe 

def monthly_means_high_elevation(d18op_data, d18ov_data, temp_data, prec_data, evap_data, relhum_data, Dataset):
    max_lon, max_lat, min_lon, min_lat = 16, 48, 3, 42 
    
    d18op_high = extract_transect(data=d18op_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=1000, maxelev=None, Dataset=Dataset)
    temp_high = extract_transect(data=temp_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=1000, maxelev=None, Dataset=Dataset)
    prec_high = extract_transect(data=prec_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=1000, maxelev=None, Dataset=Dataset)
    relhum_high = extract_transect(data=relhum_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=1000, maxelev=None, Dataset=Dataset)
    evap_high = extract_transect(data=evap_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=1000, maxelev=None, Dataset=Dataset)
    d18ov_high = extract_transect(data=d18ov_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=1000, maxelev=None, Dataset=Dataset)
    
    return d18op_high, d18ov_high, temp_high, prec_high, evap_high, relhum_high


def monthly_means_low_elevation(d18op_data, d18ov_data, temp_data, prec_data, evap_data, relhum_data, Dataset):
    max_lon, max_lat, min_lon, min_lat = 16, 51, 2, 48
    
    d18op_low = extract_transect(data=d18op_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=500, Dataset=Dataset)
    temp_low = extract_transect(data=temp_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=500, Dataset=Dataset)
    prec_low = extract_transect(data=prec_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=500, Dataset=Dataset)
    relhum_low = extract_transect(data=relhum_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=500, Dataset=Dataset)
    evap_low = extract_transect(data=evap_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=500, Dataset=Dataset)
    d18ov_low = extract_transect(data=d18ov_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                         sea_land_mask=True, minelev=None, maxelev=500, Dataset=Dataset)
    
    return d18op_low, d18ov_low, temp_low, prec_low, evap_low, relhum_low

def monthly_means_era(temp_data, prec_data, area, Dataset):
    if area == "europe":
        max_lon, max_lat, min_lon, min_lat = 25, 55, -2, 40
        temp = extract_transect(data=temp_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                             sea_land_mask=True, minelev=None, maxelev=None, Dataset=Dataset)
        prec = extract_transect(data=prec_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                             sea_land_mask=True, minelev=None, maxelev=None, Dataset=Dataset)
        
    elif area == "high":
        max_lon, max_lat, min_lon, min_lat = 16, 48, 3, 42
        temp = extract_transect(data=temp_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                                 sea_land_mask=True, minelev=1000, maxelev=None, Dataset=Dataset)
        prec = extract_transect(data=prec_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                                 sea_land_mask=True, minelev=1000, maxelev=None, Dataset=Dataset)
    
    elif area == "low":
        max_lon, max_lat, min_lon, min_lat = 16, 51, 2, 48
        temp = extract_transect(data=temp_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                             sea_land_mask=True, minelev=None, maxelev=500, Dataset=Dataset)
        prec = extract_transect(data=prec_data, maxlon=max_lon, minlon=min_lon, maxlat=max_lat, minlat=min_lat,
                                             sea_land_mask=True, minelev=None, maxelev=500, Dataset=Dataset)
    return temp, prec



def compute_mean_std(data, save=False, path_to_save=None):
    data_mean = weighted_mean(data)
    data_std = weighted_std(data)
    
    import calendar
    mnames = [calendar.month_abbr[im+1] for im in np.arange(12)]
    df = pd.DataFrame(index = mnames, columns = ["mean","std"])
    df["mean"] = data_mean
    df["std"] = data_std
    if save == True:
        df.to_csv(path_to_save)
    return df 
    
# t040
d18op_europe_t040, d18ov_europe_t040, temp_europe_t040, prec_europe_t040, evap_europe_t040, relhum_europe_t040 = monthly_means_central_europe(
    d18op_t040, d18ov_t040, temp_t040, prec_t040, evap_t040, relhum_t040, t040_data)

d18op_high_t040, d18ov_high_t040, temp_high_t040, prec_high_t040, evap_high_t040, relhum_high_t040 = monthly_means_high_elevation(
    d18op_t040, d18ov_t040, temp_t040, prec_t040, evap_t040, relhum_t040, t040_data)

d18op_low_t040, d18ov_low_t040, temp_low_t040, prec_low_t040, evap_low_t040, relhum_low_t040 = monthly_means_low_elevation(
    d18op_t040, d18ov_t040, temp_t040, prec_t040, evap_t040, relhum_t040, t040_data)

#t041
d18op_europe_t041, d18ov_europe_t041, temp_europe_t041, prec_europe_t041, evap_europe_t041, relhum_europe_t041 = monthly_means_central_europe(
    d18op_t041, d18ov_t041, temp_t041, prec_t041, evap_t041, relhum_t041, t041_data)

d18op_high_t041, d18ov_high_t041, temp_high_t041, prec_high_t041, evap_high_t041, relhum_high_t041 = monthly_means_high_elevation(
    d18op_t041, d18ov_t041, temp_t041, prec_t041, evap_t041, relhum_t041, t041_data)

d18op_low_t041, d18ov_low_t041, temp_low_t041, prec_low_t041, evap_low_t041, relhum_low_t041 = monthly_means_low_elevation(
    d18op_t041, d18ov_t041, temp_t041, prec_t041, evap_t041, relhum_t041, t041_data)

#t018
d18op_europe_t018, d18ov_europe_t018, temp_europe_t018, prec_europe_t018, evap_europe_t018, relhum_europe_t018 = monthly_means_central_europe(
    d18op_t018, d18ov_t018, temp_t018, prec_t018, evap_t018, relhum_t018, t018_data)

d18op_high_t018, d18ov_high_t018, temp_high_t018, prec_high_t018, evap_high_t018, relhum_high_t018 = monthly_means_high_elevation(
    d18op_t018, d18ov_t018, temp_t018, prec_t018, evap_t018, relhum_t018, t018_data)

d18op_low_t018, d18ov_low_t018, temp_low_t018, prec_low_t018, evap_low_t018, relhum_low_t018 = monthly_means_low_elevation(
    d18op_t018, d18ov_t018, temp_t018, prec_t018, evap_t018, relhum_t018, t018_data)

# ERA
tp = tp[:,1:,:] * 1000 #>>mm/month
t2m = t2m[:,1:,:] -273.15 

tp_mm = tp.groupby("time.month").mean("time")
t2m_mm = t2m.groupby("time.month").mean("time")

temp_ERA_eu, prec_ERA_eu = monthly_means_era(t2m_mm, tp_mm, area="europe", Dataset=t018_data) #use PI land mask to extract era
temp_ERA_high, prec_ERA_high = monthly_means_era(t2m_mm, tp_mm, area="high", Dataset=t018_data)
temp_ERA_low, prec_ERA_low = monthly_means_era(t2m_mm, tp_mm, area="low", Dataset=t018_data)

#europe

d18op_europe_t040 = compute_mean_std(d18op_europe_t040,save=True, path_to_save="mio450_d18op_eu.csv")
d18op_europe_t041 = compute_mean_std(d18op_europe_t041,save=True, path_to_save="mio278_d18op_eu.csv")
d18op_europe_t018 = compute_mean_std(d18op_europe_t018,save=True, path_to_save="PI_d18op_eu.csv")

d18ov_europe_t040 = compute_mean_std(d18ov_europe_t040,save=True, path_to_save="mio450_d18ov_eu.csv")
d18ov_europe_t041 = compute_mean_std(d18ov_europe_t041,save=True, path_to_save="mio278_d18ov_eu.csv")
d18ov_europe_t018 = compute_mean_std(d18ov_europe_t018,save=True, path_to_save="PI_d18ov_eu.csv")

temp_europe_t040 = compute_mean_std(temp_europe_t040,save=True, path_to_save="mio450_temp_eu.csv")
temp_europe_t041 = compute_mean_std(temp_europe_t041,save=True, path_to_save="mio278_temp_eu.csv")
temp_europe_t018 = compute_mean_std(temp_europe_t018,save=True, path_to_save="PI_temp_eu.csv")
temp_europe_era  = compute_mean_std(temp_ERA_eu, save=True, path_to_save="era_temp_eu.csv")

prec_europe_t040 = compute_mean_std(prec_europe_t040,save=True, path_to_save="mio450_prec_eu.csv")
prec_europe_t041 = compute_mean_std(prec_europe_t041,save=True, path_to_save="mio278_prec_eu.csv")
prec_europe_t018 = compute_mean_std(prec_europe_t018,save=True, path_to_save="PI_prec_eu.csv")
prec_europe_era  = compute_mean_std(prec_ERA_eu, save=True, path_to_save="era_prec_eu.csv")

evap_europe_t040 = compute_mean_std(evap_europe_t040,save=True, path_to_save="mio450_evap_eu.csv")
evap_europe_t041 = compute_mean_std(evap_europe_t041,save=True, path_to_save="mio278_evap_eu.csv")
evap_europe_t018 = compute_mean_std(evap_europe_t018,save=True, path_to_save="PI_evap_eu.csv")

relhum_europe_t040 = compute_mean_std(relhum_europe_t040,save=True, path_to_save="mio450_relhum_eu.csv")
relhum_europe_t041 = compute_mean_std(relhum_europe_t041,save=True, path_to_save="mio278_relhum_eu.csv")
relhum_europe_t018 = compute_mean_std(relhum_europe_t018,save=True, path_to_save="PI_relhum_eu.csv")


#high
d18op_high_t040 = compute_mean_std(d18op_high_t040,save=True, path_to_save="mio450_d18op_high.csv")
d18op_high_t041 = compute_mean_std(d18op_high_t041,save=True, path_to_save="mio278_d18op_high.csv")
d18op_high_t018 = compute_mean_std(d18op_high_t018,save=True, path_to_save="PI_d18op_high.csv")

d18ov_high_t040 = compute_mean_std(d18ov_high_t040,save=True, path_to_save="mio450_d18ov_high.csv")
d18ov_high_t041 = compute_mean_std(d18ov_high_t041,save=True, path_to_save="mio278_d18ov_high.csv")
d18ov_high_t018 = compute_mean_std(d18ov_high_t018,save=True, path_to_save="PI_d18ov_high.csv")

temp_high_t040 = compute_mean_std(temp_high_t040,save=True, path_to_save="mio450_temp_high.csv")
temp_high_t041 = compute_mean_std(temp_high_t041,save=True, path_to_save="mio278_temp_high.csv")
temp_high_t018 = compute_mean_std(temp_high_t018,save=True, path_to_save="PI_temp_high.csv")
temp_high_era  = compute_mean_std(temp_ERA_high, save=True, path_to_save="era_temp_high.csv")

prec_high_t040 = compute_mean_std(prec_high_t040,save=True, path_to_save="mio450_prec_high.csv")
prec_high_t041 = compute_mean_std(prec_high_t041,save=True, path_to_save="mio278_prec_high.csv")
prec_high_t018 = compute_mean_std(prec_high_t018,save=True, path_to_save="PI_prec_high.csv")
prec_high_era  = compute_mean_std(prec_ERA_high, save=True, path_to_save="era_prec_high.csv")

evap_high_t040 = compute_mean_std(evap_high_t040,save=True, path_to_save="mio450_evap_high.csv")
evap_high_t041 = compute_mean_std(evap_high_t041,save=True, path_to_save="mio278_evap_high.csv")
evap_high_t018 = compute_mean_std(evap_high_t018,save=True, path_to_save="PI_evap_high.csv")

relhum_high_t040 = compute_mean_std(relhum_high_t040,save=True, path_to_save="mio450_relhum_high.csv")
relhum_high_t041 = compute_mean_std(relhum_high_t041,save=True, path_to_save="mio278_relhum_high.csv")
relhum_high_t018 = compute_mean_std(relhum_high_t018,save=True, path_to_save="PI_relhum_high.csv")

#low
d18op_low_t040 = compute_mean_std(d18op_low_t040,save=True, path_to_save="mio450_d18op_low.csv")
d18op_low_t041 = compute_mean_std(d18op_low_t041,save=True, path_to_save="mio278_d18op_low.csv")
d18op_low_t018 = compute_mean_std(d18op_low_t018,save=True, path_to_save="PI_d18op_low.csv")

d18ov_low_t040 = compute_mean_std(d18ov_low_t040,save=True, path_to_save="mio450_d18ov_low.csv")
d18ov_low_t041 = compute_mean_std(d18ov_low_t041,save=True, path_to_save="mio278_d18ov_low.csv")
d18ov_low_t018 = compute_mean_std(d18ov_low_t018,save=True, path_to_save="PI_d18ov_low.csv")

temp_low_t040 = compute_mean_std(temp_low_t040,save=True, path_to_save="mio450_temp_low.csv")
temp_low_t041 = compute_mean_std(temp_low_t041,save=True, path_to_save="mio278_temp_low.csv")
temp_low_t018 = compute_mean_std(temp_low_t018,save=True, path_to_save="PI_temp_low.csv")
temp_low_era  = compute_mean_std(temp_ERA_low, save=True, path_to_save="era_temp_low.csv")

prec_low_t040 = compute_mean_std(prec_low_t040,save=True, path_to_save="mio450_prec_low.csv")
prec_low_t041 = compute_mean_std(prec_low_t041,save=True, path_to_save="mio278_prec_low.csv")
prec_low_t018 = compute_mean_std(prec_low_t018,save=True, path_to_save="PI_prec_low.csv")
prec_low_era  = compute_mean_std(prec_ERA_low, save=True, path_to_save="era_prec_low.csv")

evap_low_t040 = compute_mean_std(evap_low_t040,save=True, path_to_save="mio450_evap_low.csv")
evap_low_t041 = compute_mean_std(evap_low_t041,save=True, path_to_save="mio278_evap_low.csv")
evap_low_t018 = compute_mean_std(evap_low_t018,save=True, path_to_save="PI_evap_low.csv")

relhum_low_t040 = compute_mean_std(relhum_low_t040,save=True, path_to_save="mio450_relhum_low.csv")
relhum_low_t041 = compute_mean_std(relhum_low_t041,save=True, path_to_save="mio278_relhum_low.csv")
relhum_low_t018 = compute_mean_std(relhum_low_t018,save=True, path_to_save="PI_relhum_low.csv")



def apply_style(fontsize=22):
    small_size=fontsize-2
    #plt.style.use('seaborn')
#    plt.style.use('dark_background')
    #plt.style.use("fivethirtyeight")
    plt.style.use('bmh')    
    plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
    mpl.rc('text', usetex=True)
    mpl.rc('font', size=22, family='serif')
    mpl.rc('xtick', labelsize=small_size)
    mpl.rc('ytick', labelsize=small_size)
    mpl.rc('legend', fontsize=small_size)
    mpl.rc('axes', labelsize=20)
    mpl.rc('lines', linewidth=2.5)
    mpl.rc("font", weight="bold")


def plot_all_data(data_Mio278, data_Mio450, data_PI, axes, text, ymin, ymax, label, ERA_data=None):
    axes.plot(data_Mio278["mean"], "--", color=green, label="Mio278", linewidth=2)
    axes.plot(data_Mio450["mean"], "--", color=red, label="Mio450", linewidth=2)
    axes.plot(data_PI["mean"], "--", color=black, label="PI", linewidth=2)
    if ERA_data is not None:
        axes.plot(ERA_data["mean"], "--", color=blue, label="ERA", linewidth=2)
    else:
        None
    axes.set_ylabel(text, fontsize= 20, fontweight="bold")
    axes.set_ylim(bottom=ymin, top=ymax)
    axes.set_title(label, fontdict= {"fontsize": 15, "fontweight":"bold"}, loc= "left")
    Mio278_low, Mio278_high = data_Mio278["mean"]- data_Mio278["std"], data_Mio278["mean"] + data_Mio278["std"]
    Mio450_low, Mio450_high = data_Mio450["mean"] - data_Mio450["std"], data_Mio450["mean"] + data_Mio450["std"]
    PI_low, PI_high = data_PI["mean"] - data_PI["std"], data_PI["mean"] + data_PI["std"]
    
    
    axes.fill_between(data_Mio278.index, Mio278_low, Mio278_high, color=green, alpha=0.1)
    axes.fill_between(data_Mio450.index, Mio450_low, Mio450_high, color=red, alpha=0.1)
    axes.fill_between(data_PI.index, PI_low, PI_high, color=black, alpha=0.1)


apply_style()
fig, ((ax1,ax2),(ax3,ax4), (ax5,ax6)) = plt.subplots(3, 2, sharex = False, sharey = False, figsize=(20,15))
plt.subplots_adjust(left=0.12, right=1-0.01, top=0.95, bottom=0.06, hspace=0.01)

plot_all_data(data_Mio278=temp_europe_t041, data_Mio450=temp_europe_t040, data_PI=temp_europe_t018, ERA_data=temp_europe_era, axes=ax1, text="Temperature [°C]", ymin=-5, ymax=35, label="(A)")
plot_all_data(data_Mio278=prec_europe_t041, data_Mio450=prec_europe_t040, data_PI=prec_europe_t018, ERA_data=prec_europe_era, axes=ax2, text="Precipitation [mm/month]", ymin=0, ymax=125, label="(B)")
plot_all_data(data_Mio278=relhum_europe_t041, data_Mio450=relhum_europe_t040, data_PI=relhum_europe_t018, ERA_data=None, axes=ax3, text="Relative Humidity [\%]", ymin=30, ymax=100, label="(C)")
plot_all_data(data_Mio278=evap_europe_t041, data_Mio450=evap_europe_t040, data_PI=evap_europe_t018, ERA_data=None, axes=ax4, text="Evaporation [mm/month]", ymin=0, ymax=140, label="(D)")
plot_all_data(data_Mio278=d18op_europe_t041, data_Mio450=d18op_europe_t040, data_PI=d18op_europe_t018, ERA_data=None, axes=ax6, text=u'$\delta^{18}$Op ‰ vs SMOW', ymin=-13, ymax=-2, label="(E)")
plot_all_data(data_Mio278=d18ov_europe_t041, data_Mio450=d18ov_europe_t040, data_PI=d18ov_europe_t018, ERA_data=None, axes=ax5, text=u'$\delta^{18}$Ov ‰ vs SMOW', ymin=-21, ymax=-12, label="(F)")

ax2.legend(bbox_to_anchor=(1.04,1),frameon=True, fontsize=11, loc="upper left")
plt.tight_layout()
plt.savefig("Europe_final.pdf", format="pdf", bbox_inches="tight")


fig, ((ax1,ax2),(ax3,ax4), (ax5,ax6)) = plt.subplots(3, 2, sharex = False, sharey = False, figsize=(20,15))
plt.subplots_adjust(left=0.12, right=1-0.01, top=0.95, bottom=0.06, hspace=0.01)

plot_all_data(data_Mio278=temp_high_t041, data_Mio450=temp_high_t040, data_PI=temp_high_t018, ERA_data=temp_high_era, axes=ax1, text="Temperature [°C]", ymin=-10, ymax=25, label="(A)")
plot_all_data(data_Mio278=prec_high_t041, data_Mio450=prec_high_t040, data_PI=prec_high_t018, ERA_data=prec_high_era, axes=ax2, text="Precipitation [mm/month]", ymin=0, ymax=200, label="(B)")
plot_all_data(data_Mio278=relhum_high_t041, data_Mio450=relhum_high_t040, data_PI=relhum_high_t018, ERA_data=None, axes=ax3, text="Relative Humidity [\%]", ymin=40, ymax=100, label="(C)")
plot_all_data(data_Mio278=evap_high_t041, data_Mio450=evap_high_t040, data_PI=evap_high_t018, ERA_data=None, axes=ax4, text="Evaporation [mm/month]", ymin=0, ymax=140, label="(D)")
plot_all_data(data_Mio278=d18op_high_t041, data_Mio450=d18op_high_t040, data_PI=d18op_high_t018, ERA_data=None, axes=ax6, text=u'$\delta^{18}$Op ‰ vs SMOW', ymin=-16, ymax=-2, label="(E)")
plot_all_data(data_Mio278=d18ov_high_t041, data_Mio450=d18ov_high_t040, data_PI=d18ov_high_t018, ERA_data=None, axes=ax5, text=u'$\delta^{18}$Ov ‰ vs SMOW', ymin=-25, ymax=-13, label="(F)")

ax2.legend(bbox_to_anchor=(1.04,1),frameon=True, fontsize=11, loc="upper left")
plt.tight_layout()
plt.savefig("High_elevation_final.pdf", format="pdf", bbox_inches="tight")


fig, ((ax1,ax2),(ax3,ax4), (ax5,ax6)) = plt.subplots(3, 2, sharex = False, sharey = False, figsize=(20,15))
plt.subplots_adjust(left=0.12, right=1-0.01, top=0.95, bottom=0.06, hspace=0.01)

plot_all_data(data_Mio278=temp_low_t041, data_Mio450=temp_low_t040, data_PI=temp_low_t018, ERA_data=temp_low_era, axes=ax1, text="Temperature [°C]", ymin=-10, ymax=25, label="(A)")
plot_all_data(data_Mio278=prec_low_t041, data_Mio450=prec_low_t040, data_PI=prec_low_t018, ERA_data=prec_low_era, axes=ax2, text="Precipitation [mm/month]", ymin=0, ymax=150, label="(B)")
plot_all_data(data_Mio278=relhum_low_t041, data_Mio450=relhum_low_t040, data_PI=relhum_low_t018, ERA_data=None, axes=ax3, text="Relative Humidity [\%]", ymin=60, ymax=100, label="(C)")
plot_all_data(data_Mio278=evap_low_t041, data_Mio450=evap_low_t040, data_PI=evap_low_t018, ERA_data=None, axes=ax4, text="Evaporation [mm/month]", ymin=0, ymax=140, label="(D)")
plot_all_data(data_Mio278=d18op_low_t041, data_Mio450=d18op_low_t040, data_PI=d18op_low_t018, ERA_data=None, axes=ax6, text=u'$\delta^{18}$Op ‰ vs SMOW', ymin=-12, ymax=-2, label="(E)")
plot_all_data(data_Mio278=d18ov_low_t041, data_Mio450=d18ov_low_t040, data_PI=d18ov_low_t018, ERA_data=None, axes=ax5, text=u'$\delta^{18}$Ov ‰ vs SMOW', ymin=-20, ymax=-13, label="(F)")

ax2.legend(bbox_to_anchor=(1.04,1),frameon=True, fontsize=11, loc="upper left")
plt.tight_layout()
plt.savefig("low_elevation_final.pdf", format="pdf", bbox_inches="tight")
