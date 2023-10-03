# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 09:00:22 2023

@author: dboateng
"""

# import modules 
import os
import numpy as np
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt 
import netCDF4 as nc 
import glob
from datetime import datetime

# set paths 

main_path_to_data_oco3 = "D:/Datasets/OCO3/Data"
glob_name_oco3 = "oco3_LtCO2_20*"

main_path_to_data_oco2 = "D:/Datasets/OCO2/Data"
glob_name_oco2 = "oco2_LtCO2_20*"

# function to detect the number of files to load

path_to_csv_files_oco3 = os.path.join(main_path_to_data_oco3, "csv_files")

path_to_csv_files_oco2 = os.path.join(main_path_to_data_oco2, "csv_files")

if not os.path.exists(path_to_csv_files_oco3):
    os.makedirs(path_to_csv_files_oco3)
    
if not os.path.exists(path_to_csv_files_oco2):
    os.makedirs(path_to_csv_files_oco2)

def count_files_in_directory(path_to_data, glob_pattern="*.nc4"):
    files = glob.glob(path_to_data + "/" + glob_pattern)
    print(len(files))
    return files


# function to convert the nc file to csv

def conv_date(d):
    return datetime.strptime(str(d), '%Y%m%d%H%M%S%f')

def convHdf(path_file, frames_folder):
    data = nc.Dataset(path_file)

    df_xco2 = pd.DataFrame({
        'Xco2': data.variables['xco2'][:],
        'Latitude': data.variables['latitude'][:],
        'Longitude': data.variables['longitude'][:],
        'quality_flag': data.variables['xco2_quality_flag'][:],
        'DateTime': [conv_date(d) for d in data.variables['sounding_id'][:]],
    })

    df_xco2['DateTime'] = pd.to_datetime(df_xco2['DateTime'])
    df_xco2['Year'] = df_xco2['DateTime'].dt.year
    df_xco2['Month'] = df_xco2['DateTime'].dt.month
    df_xco2['Day'] = df_xco2['DateTime'].dt.day

    df_xco2 = df_xco2[df_xco2['quality_flag'] == 0]

    date = str(data.variables['sounding_id'][0])
    output_file = os.path.join(frames_folder, f'_xco2_{date}_.csv')
    df_xco2.to_csv(output_file, index=False)




if __name__=="__main__":
    files_oco3 = count_files_in_directory(main_path_to_data_oco3, glob_name_oco3)
    files_oco2 = count_files_in_directory(main_path_to_data_oco2, glob_name_oco2)
    
    
    for file in files_oco3:
        convHdf(path_file=file, frames_folder=path_to_csv_files_oco3)
        
    for file in files_oco2:
        convHdf(path_file=file, frames_folder=path_to_csv_files_oco2)
        
        
    
    