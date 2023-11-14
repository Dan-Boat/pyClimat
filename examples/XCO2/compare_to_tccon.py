# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 13:20:05 2023

@author: dboateng
"""

import os
import numpy as np
import pandas as pd 
import xarray as xr
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from math import radians, sin, cos, sqrt, atan2

path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/scratch/OCO2_OCO3/OCO2_OCO3_grid_1"

ds = xr.open_dataset("D:/Datasets/OCO3/monthly_gridded/raw_OCO2_OCO3_gridded_monthly_1.0_grid.nc")

main_path_to_data_tccon = "D:/Datasets/OCO2/TCCON/"



#path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/scratch/OCO2"

#ds = xr.open_dataset("D:/Datasets/OCO3/monthly_gridded/raw_OCO2_gridded_monthly_1.0_grid.nc")



# def haversine(lon1, lat1, lon2, lat2):
#     """
#     Calculate the great circle distance in kilometers between two points 
#     on the earth (specified in decimal degrees)
#     """
#     # Convert decimal degrees to radians
#     lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

#     # Haversine formula
#     dlon = lon2 - lon1
#     dlat = lat2 - lat1
#     a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
#     c = 2 * atan2(sqrt(a), sqrt(1 - a))
#     radius = 6371  # Radius of Earth in kilometers
#     distance = radius * c

#     return distance

def haversine(lon1, lat1, lon2, lat2):
# convert decimal degrees to radians 
    lon1 = np.deg2rad(lon1)
    lon2 = np.deg2rad(lon2)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6371
    return c * r    #km


def extract_indices_around(dataset, lat, lon, radius):
    close_grids = lambda lat_, lon_: haversine(lat, lon, lat_, lon_) <= radius
    
    
    if hasattr(dataset, "longitude"):
        LON, LAT = np.meshgrid(dataset.longitude, dataset.latitude)
    else:
        LON, LAT = np.meshgrid(dataset.lon, dataset.lat)
        
    grids_index = np.where(close_grids(LAT, LON))
    
    return grids_index

def plot_data(df, title, fig_path):
    
    from pyClimat.plot_utils import apply_style
    from matplotlib.dates import YearLocator
    
    apply_style(fontsize=25, style="seaborn-talk", linewidth=4,)
    
    fig, ax = plt.subplots(1, 1, figsize= (20, 15), sharex=True)
    ax.plot(df["TCCON"], linestyle="-", color="black", label="TCCON")
    ax.plot(df["OCO2"], linestyle="-", color="red", label="OCO2-OCO3 L2")
    ax.legend(frameon = True, fontsize=20)
    ax.set_ylabel("XCO2 [ppm]", fontweight="bold", fontsize=20)
    ax.set_ylim([400, 430])
    ax.set_title(title, fontsize=20, fontweight="bold", loc="left") 
    plt.tight_layout()
    
    fig_name = title.replace(",","")
    fig_name = fig_name.replace(" ","_") + ".png"
    
    
    plt.savefig(os.path.join(fig_path, fig_name), bbox_inches="tight", format= "png")
    
    
    pass


def calculate_regional_means(ds, lon_target, lat_target, radius_deg, method="old"):
    """
    Calculate regional means around a specific longitude and latitude location
    with a given radius for a NetCDF dataset using xarray.
    """
    # Find indices of the nearest grid point to the target location
    
    if hasattr(ds, "longitude"):
        ds = ds.rename({"longitude":"lon", "latitude":"lat"})
        
    ds = ds.assign_coords({"lon": (((ds.lon + 180) % 360) - 180)})
    
    if method=="old":
    
        indices = extract_indices_around(ds, lat_target, lon_target, radius_deg)
        
        regional_mean = ds.isel(lat=indices[0], lon=indices[1]).mean(dim=("lon", "lat")).data
        
    else:
        
        
        lon_diff = np.abs(ds.lon - lon_target)
        lat_diff = np.abs(ds.lat - lat_target)
        idx_lon = lon_diff.argmin().item()
        idx_lat = lat_diff.argmin().item()
    
        # Create a mask based on the haversine distance within the specified radius
        mask = np.zeros_like(ds.lon, dtype=bool)
        for i in range(ds.lon.size):
            for j in range(ds.lat.size):
                lon_i = ds.lon[i].item()
                if lon_i > 180:
                    lon_i -= 360
                distance = haversine(lon_target, lat_target, lon_i, ds.lat[j].item())
                if distance <= radius_deg:
                    mask[i] = True
    
        # Apply the mask to the dataset
        ds_region = ds.sel(lon=ds.lon[mask])
    
        # Calculate the regional mean
        regional_mean = ds_region.mean(dim='lon')

    return np.float64(regional_mean)


def count_files_in_directory(path_to_data, glob_pattern="*.nc"):
    files = glob.glob(path_to_data + "/" + glob_pattern)
    return files   


files = count_files_in_directory(main_path_to_data_tccon)


dfs_tccon = []
dfs_oco2 = []


for file in files:
    tccon_data  = xr.open_dataset(os.path.join(main_path_to_data_tccon, file))
    
    lat = tccon_data["lat"][0].data
    lon = tccon_data["long"][0].data
    
    
    df_xco2 = pd.DataFrame({
                            'Xco2': tccon_data['xco2'],
                            "time": pd.to_datetime(tccon_data["time"]),
                            "year": tccon_data["year"],
                            })
    
    df = df_xco2.set_index("time")
    
    #compute monthly
    df_month = df.resample("MS").mean()
    
    df = pd.DataFrame(columns=["OCO2", "TCCON"], index=df_month.index.values)
    
    df_all = pd.DataFrame(columns=["TCCON", "OCO2"]) 
    
    
    stn_name = tccon_data.short_location
    
    
    
    for i in np.arange(ds.time.shape[0]):
        data = ds.XCO2[i]
        
        if ds.time.data[i] in df_month.index.values:
            
            if np.isnan(df_month.loc[ds.time.data[i]]["Xco2"]) ==False:
                df.loc[ds.time.data[i]]["TCCON"] = df_month.loc[ds.time.data[i]]["Xco2"]
                
                df.loc[ds.time.data[i]]["OCO2"] = calculate_regional_means(ds=data, 
                                                                           lon_target=lon, 
                                                                           lat_target=lat, 
                                                                           radius_deg=100,
                                                                           method="old")
    plot_data(df=df.dropna(), title=stn_name, fig_path=path_to_plots)

    df = df.dropna()
    
    dfs_oco2.append(df["OCO2"])
    dfs_tccon.append(df["TCCON"])
    
result_oco2 = pd.concat(dfs_oco2, ignore_index=True)
results_tccon = pd.concat(dfs_tccon, ignore_index=True)


from scipy import stats  
import seaborn as sns 
import statsmodels.api as sm
from sklearn.metrics import mean_squared_error

x = pd.to_numeric(result_oco2, errors='coerce') 
y= pd.to_numeric(results_tccon, errors='coerce')

slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
y_pred = intercept + slope * x

mse = mean_squared_error(y_true=y, y_pred=y_pred, squared=True)
rmse = mean_squared_error(y_true=y, y_pred=y_pred, squared=False)


xx = sm.add_constant(x)
model = sm.OLS(y, xx)
r = model.fit()
predictions = r.get_prediction(exog=xx, transform=False).summary_frame(alpha=0.05)


label = "RSME: {:.2f}, MSE: {:.2f}, r²={:.2f}".format(rmse, mse,r_value)

from pyClimat.plot_utils import apply_style
apply_style(fontsize=25, style="seaborn-talk", linewidth=4,)

            
fig,ax = plt.subplots(1,1, sharex=False, figsize=(20, 15))
sns.regplot(x=x, y=y, marker="*", scatter_kws={"color":"black", "s":200, "alpha":0.7},
            color="black",ax=ax, label=label)

ax.plot(x, predictions["obs_ci_lower"], linestyle="-",
                linewidth=0.3, color="black")
ax.plot(x, predictions["obs_ci_upper"], linestyle="-",
                linewidth=0.3, color="black")

ax.set_xlabel("OCO2-OCO3 L2 gridded monthly", fontweight="bold", fontsize=28)
ax.set_ylabel("TCCON sites measurements", fontweight="bold", fontsize=28)

ax.legend(frameon=True, fontsize=24,
                  loc="upper right", borderaxespad=0,
                  )
ax.set_box_aspect(1)
plt.tight_layout()
plt.savefig(os.path.join(path_to_plots, "compare_tccon_oco2.png"), bbox_inches="tight", format= "png")





