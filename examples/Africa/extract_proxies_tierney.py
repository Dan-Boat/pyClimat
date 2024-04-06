# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 16:48:07 2024

@author: dboateng
"""

import os 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

from pyClimat.data import read_from_path
from pyClimat.analysis import compute_lterm_mean, extract_var, extract_profile, compute_lterm_diff, extract_transect, extract_profile
from pyClimat.utils import haversine, extract_indices_around
from pyClimat.plots import plot_annual_mean
from pyClimat.plot_utils import *


path_to_data = "C:/Users/dboateng/Desktop/Datasets/WAM reconstructions/"
tierney_data = "tierney_data.csv"

#df = pd.read_csv(os.path.join(path_to_data, tierney_data))

main_path = "D:/Datasets/Model_output_pst/"
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")
filename_lterm = "1003_1017_1m_mlterm.nc"


PI_data = read_from_path(pi_path, filename_lterm)
pr_pi = extract_var(Dataset=PI_data, varname="prec", units="mm/a")
pr_pi_alt = compute_lterm_mean(data=pr_pi, time="annual")


def calculate_regional_means(ds, lon_target, lat_target, radius_deg,):
    """
    Calculate regional means around a specific longitude and latitude location
    with a given radius for a NetCDF dataset using xarray.
    """
    # Find indices of the nearest grid point to the target location
    
    if hasattr(ds, "longitude"):
        ds = ds.rename({"longitude":"lon", "latitude":"lat"})
        
    ds = ds.assign_coords({"lon": (((ds.lon + 180) % 360) - 180)})
    
    indices = extract_indices_around(ds, lat_target, lon_target, radius_deg)
    
    regional_mean = ds.isel(lat=indices[0], lon=indices[1]).mean(dim=("lon", "lat")).data
        
    return np.float64(regional_mean)

mean = calculate_regional_means(ds=pr_pi_alt, lon_target=-17.282, lat_target=19.363, radius_deg=50)
print(mean)

def read_tierney():
    filename = "tierney2017gc68.txt"
    record = os.path.join(path_to_data, filename)
    
    # Read the file line by line and filter out lines starting with '#'
    lines = []
    with open(record, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                lines.append(line)
    
    # Use pandas to read the filtered lines as a dataframe
    data = pd.DataFrame([line.strip().split() for line in lines])
    
    data.columns = data.loc[0]
    data = data[1:]
    
    
    data['age_BP'] = data['age_BP'].astype('int')
    data['precip'] = data['precip'].astype('int')
    data['precip_1s_lower'] = data['precip_1s_lower'].astype('int')
    data_mid_holocene = data.loc[(data['age_BP'] >= 5000) & (data['age_BP'] <= 7000)]
    print(data_mid_holocene['precip'].mean())
    
    print(data_mid_holocene['precip_1s_lower'].mean())
    
    print(data_mid_holocene['precip'].std())
    
    
    

