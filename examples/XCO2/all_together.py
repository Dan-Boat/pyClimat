# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:40:56 2023

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

import dask.dataframe as dd
import dask.array as da


#utils function
def conv_date(d):
    return datetime.strptime(str(d), '%Y%m%d%H%M%S%f')

def filter_grid_cells(grouped_df, min_data_points):
    if len(grouped_df) > min_data_points:
        return grouped_df['XCO2'].mean()
    else:
        return np.nan


# function to read the monthly nc files from path to csv

def read_nc_files_for_year_to_csv(path_nc_data, glob_name, path_save_csv, data_type):
    
    nc_files_names = glob.glob(path_nc_data + "/" + glob_name)
    
    for file in nc_files_names:
        dataset = nc.Dataset(file)
        
        df_XCO2 = pd.DataFrame({
            'XCO2': dataset.variables['xco2'][:],
            'Latitude': dataset.variables['latitude'][:],
            'Longitude': dataset.variables['longitude'][:],
            'quality_flag': dataset.variables['xco2_quality_flag'][:],
            'time': [conv_date(d) for d in dataset.variables['sounding_id'][:]],
        })
        
        df_XCO2["time"] = pd.to_datetime(df_XCO2["time"])
        df_XCO2["year"] = df_XCO2["time"].dt.year
        df_XCO2["month"] = df_XCO2["time"].dt.month
        
        df_XCO2 = df_XCO2[df_XCO2["quality_flag"] == 0]
        
        if not os.path.exists(path_save_csv):
            os.makedirs(path_save_csv)
            
        filename = str(df_XCO2["time"].values[1])[:10].replace("-", "_") + "_" + data_type + ".csv"
        
        df_XCO2.to_csv(os.path.join(path_save_csv, filename), index=False)
        
        
# function that grid to 0.5 by 0.5 monthly and then save the files nc file    

def read_csv_files_and_grid_to_nc_monhtly(path_data_OCO2, path_data_OCO3, glob_name, grid_interval,
                                          min_data_points,  datetime, path_save_nc_data):
    
    df_OCO2 = dd.read_csv(path_data_OCO2 + "/" + glob_name, parse_dates=["DateTime"])
    df_OCO3 = dd.read_csv(path_data_OCO3 + "/" + glob_name, parse_dates=["DateTime"])
    
    df_OCO2_3 = dd.concat([df_OCO2, df_OCO3], axis=0, interleave_partitions=True)
    
    x_grid_edges = da.arange(-180, 181, grid_interval)
    y_grid_edges = da.arange(-90, 91, grid_interval)

    x_grid_labels = (x_grid_edges[:-1] + x_grid_edges[1:]) / 2
    y_grid_labels = (y_grid_edges[:-1] + y_grid_edges[1:]) / 2
    
    df = pd.DataFrame()
    
    df['Longitude'] = pd.cut(df_OCO2_3['Longitude'].values.compute(), bins=x_grid_edges.compute(), labels=x_grid_labels.compute())
    df['Latitude'] = pd.cut(df_OCO2_3['Latitude'].values.compute(), bins=y_grid_edges.compute(), labels=y_grid_labels.compute())
    df['XCO2'] = df_OCO2_3['Xco2'].values.compute()
    
    df_bined = df.groupby(['Longitude', 'Latitude']).apply(filter_grid_cells, min_data_points=min_data_points).reset_index(name='XCO2')
    
    
    if grid_interval == 0.5:
        
        lat = np.arange(-89.75, 90.25, 0.5)
        lon = np.arange(-179.75, 180.25, 0.5)
        lon, lat = np.meshgrid(lon, lat)
        lat = lat.flatten()
        lon = lon.flatten()
        XCO2 = [np.nan] * len(lat)
        

        df_dimension = pd.DataFrame({'Latitude': lat, 'Longitude': lon, 'XCO2': XCO2})
        
    elif grid_interval == 1:
        
        lat = np.arange(-89.5, 90.5, 1)
        lon = np.arange(-179.5, 180.5, 1)
        lon, lat = np.meshgrid(lon, lat)
        lat = lat.flatten()
        lon = lon.flatten()
        XCO2 = [np.nan] * len(lat)

        df_dimension = pd.DataFrame({'Latitude': lat, 'Longitude': lon, 'XCO2': XCO2})
    
    df_bined = df_bined.dropna()
    
    # Create a dictionary mapping (latitude, longitude) to new XCO2 values
    mapping = dict(zip(zip(df_bined['Latitude'], df_bined['Longitude']), df_bined['XCO2']))

    # Update XCO2 values in df1 using the mapping
    df_dimension['XCO2'] = df_dimension.apply(lambda row: mapping.get((row['Latitude'], row['Longitude']), row['XCO2']), axis=1)

    # Create an xarray dataset from the modified DataFrame
    # ds = xr.Dataset({
    #     'XCO2': (('latitude', 'longitude'), df_dimension['XCO2'].values.reshape(360, 720))
    #                 }, coords={'latitude': np.arange(-89.75, 90.25, 0.5), 'longitude': np.arange(-179.75, 180.25, 0.5),
    #                             'time': datetime})
                               
                               
    # Create an xarray dataset from the modified DataFrame
    ds = xr.Dataset({
        'XCO2': (('latitude', 'longitude'), df_dimension['XCO2'].values.reshape(180, 360))
                    }, coords={'latitude': np.arange(-89.5, 90.5, 1), 'longitude': np.arange(-179.5, 180.5, 1),
                                'time': datetime})
                               
                               
    ds = xr.Dataset({
        'XCO2': (('time', 'latitude', 'longitude'), df_dimension['XCO2'].values.reshape((1, 180, 360))
        )}, coords={
            'latitude': np.arange(-89.5, 90.5, 1),
            'longitude': np.arange(-179.5, 180.5, 1),
            'time': np.atleast_1d(datetime)
        })
                               
    filename = str(datetime)[:10].replace("-", "_")+"_OCO2_3.nc"
    
    
    ds.to_netcdf(os.path.join(path_save_nc_data, filename))
               
        
        
        
    
if __name__=="__main__":
    main_path_to_data_oco3 = "D:/Datasets/OCO3/Data"
    glob_name_oco3 = "oco3_LtCO2_20*"
    path_to_csv_save = os.path.join(main_path_to_data_oco3, "csv_files")
    
    path_csv_oco2 = "D:/Datasets/OCO2/Data/csv_files"

    path_csv_oco3 = "D:/Datasets/OCO3/Data/csv_files"
    
    time = pd.date_range(start="2020-01-01", end="2020-12-31", freq="MS")
    
    glob_name = "_xco2_202001*csv"
    # read_nc_files_for_year_to_csv(path_nc_data=main_path_to_data_oco3, 
    #                               glob_name=glob_name_oco3, 
    #                               path_save_csv=path_to_csv_save,
    #                               data_type="OCO3")
    
    read_csv_files_and_grid_to_nc_monhtly(path_data_OCO2=path_csv_oco2, path_data_OCO3=path_csv_oco3, 
                                          glob_name=glob_name, grid_interval=1,
                                              min_data_points=20,  datetime=time[0],
                                              path_save_nc_data=main_path_to_data_oco3)