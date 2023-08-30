# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 12:48:40 2023

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
# define paths
main_path_to_data = "D:/Datasets/OCO2/TCCON/"
file = "br20090106_20210624.public.qc.nc"



# function to detect the number of files to load

def count_files_in_directory(path_to_data, glob_pattern="*.nc"):
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
        'Latitude': data.variables['lat'][:],
        'Longitude': data.variables['long'][:],
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
    
    
    for i in np.arange(1, 13):
        if i <10:
            id_num = "0"+str(i)
        else:
            id_num = str(i)
            
        files = count_files_in_directory(main_path_to_data)
        
        for file in files:
            tccon_data  = xr.open_dataset(os.path.join(main_path_to_data, file))
            
            lat = tccon_data["lat"][0].data
            lon = tccon_data["long"][0].data
        
        
            df_xco2 = pd.DataFrame({
                'Xco2': tccon_data['xco2'],
                "time": pd.to_datetime(tccon_data["time"]),
                "year": tccon_data["year"],
                })
        
            df = df_xco2.set_index("time")
        
            df_2015 = df.loc[df["year"]==2015]
            
            #df_2015["month"] = df_2015.index.month.values 
            
            
            if df_2015.empty == False:
                
                df_2015_month = df_2015.resample("M").mean()
                
                df_2015_month["month"] = df_2015_month.index.month.values 
                
                
                # check if the data has the whole year
                if i in df_2015_month.month.values:
                    mon_value = df_2015_month.loc[df_2015_month["month"]==i]["Xco2"].values
                    print(mon_value)
                
            else:
                print("No data for "+  str(file))
            
        
        # for file in files:
        #     convHdf(path_file=file, frames_folder=os.path.join(main_path_to_data, "csv_files"))