# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 14:09:55 2023

@author: dboateng

This script deals with organizing the dataset downloaded from the WISER portal
1. extracting the station attributes
2. extracting the relevant variables
3. Grouping them for the required analysis interms of the regional means
"""

import os 
import pandas as pd 
import numpy as np 
import xlrd 
import sys
from pathlib import Path as p



path_to_data = "D:/Datasets/GNIP_data/time_series/processed"
path_to_processed = "D:/Datasets/GNIP_data/processed"
path_to_processed_all = "D:/Datasets/GNIP_data/processed_all"

#use glob to extract all the data in path

try_path = "D:/Datasets/GNIP_data/time_series/processed/AUSTRIA_GRAZ (UNIVERSITAET).csv"



def preprocess_to_long_records():
    
    df_info = pd.DataFrame(index=np.arange(200), columns=["Name", "lat", "lon", "elev", "country","years"])

    daterange = pd.date_range(start="1960-01-01", end="2021-12-01", freq="MS") # time range for ERA5
    
    
    for i,file in enumerate(p(path_to_data).glob("*.csv")):
        
        data = pd.read_csv(file)
        name = data["Site"][0].replace(" ", "_")
        lat = data["Latitude"][0]
        lon = data["Longitude"][0]
        elev = data["Altitude"][0]
        country = data["Country"][0]
    
        
        use_columns = ["Date", "O18", "H2", "H3", "Precipitation", "Air Temperature", "Vapour Pressure"]
        data_to_save = pd.read_csv(file, usecols=use_columns, parse_dates=["Date"], index_col=False)
        
        if data_to_save.Date[0].is_month_start == False:
            data_to_save["Date"] = data_to_save.Date.values.astype("datetime64[M]")
            
            data_to_save = data_to_save.set_index(["Date"], drop=True)
            
        
        if len(data_to_save)/12 >= 30:
            
            df_info["Name"][i] = name
            df_info["lat"][i] = lat
            df_info["lon"][i] = lon
            df_info["elev"][i] = elev
            df_info["country"][i] = country
            df_info["years"][i] = len(data_to_save["O18"].dropna())/12
            
            
            data_to_save = data_to_save.replace(np.nan, -9999)
            
            print("Data long enough to be save into processed")
            
            data_to_save.to_csv(os.path.join(path_to_processed_all, name + ".csv"), index=True)
    
    
    df_info = df_info.rename(columns = {"lon":"Longitude",
                                        "lat":"Latitude","elev":"Elevation"})
    
    
    df_info = df_info.sort_values(by=["Name"], ascending=True)
    df_info = df_info.reset_index().dropna().set_index('index')
    df_info.to_csv(os.path.join(path_to_processed, "station_info.csv"), index=False)


def process_data_for_variable(path_to_data, path_to_store):
    
    varname  = "Temperature"
    use_colums = ["Date","Air Temperature"]

    daterange = pd.date_range(start="1960-01-01", end="2021-12-01", freq="MS")
    
    df_to_store = pd.DataFrame(columns = ["Time", varname])
    df_to_store["Time"] = daterange
    df_to_store = df_to_store.set_index(["Time"], drop=True)
    
    
    
    glob_name = "*.csv"
    
    for csv in p(path_to_data).glob(glob_name):
    
        
        print(csv.name)
    
        # df_stats = df.describe()
        # df_nans = df.isna().sum()
    
        df = pd.read_csv(csv, usecols=use_colums,
                         parse_dates=["Date"], dayfirst=True)
        
        
        df = df.set_index(["Date"], drop=True)
        
        # select the data that is part of the dates
        df = df.loc[daterange[0]:daterange[-1]]
        
        
        df_to_store[varname][df.index] = df["Air Temperature"][df.index] 
        
        
        df_to_store = df_to_store.replace(np.nan, -9999)
        
        df_to_store.to_csv(os.path.join(path_to_store, csv.name), index=True) 
        
        

def add_info_to_data(path_to_info, path_to_data, path_to_store, glob_name="*.csv", 
                     varname="H2"):
    """
    This function locate the data info in data_to_info by using the station code and then append it to the start of the data 
    This function also stores a summary info file (station names and station loc) that can be used to interate all the stations when applying the downscaling package
    Parameters
    ----------
    path_to_info : TYPE: str
        DESCRIPTION. Path to the data containing all the station infomation
    path_to_data : TYPE: str
        DESCRIPTION. Path to data required for appending info
    glob_name : TYPE: str
        DESCRIPTION. The global pattern used to extract data (eg. *.csv)
    Returns
    -------
    None.
    """
    info_cols = ["ID", "Name", "lat", "lon", "elev", "country",]
    
    
    df_info = pd.DataFrame(columns= info_cols[1:])
    
    for i, csv_file in enumerate(p(path_to_data).glob(glob_name)):
        sep_filename = csv_file.name.split(sep="_")
        csv_id = i + 1
        
       # data = pd.read_csv(csv_file)
        data_info = pd.read_csv(path_to_info, usecols=info_cols)
        data_info = data_info.set_index(["ID"])
        
        df_info = df_info.append(data_info.loc[csv_id])
        
        csv_info = data_info.loc[csv_id]
        
        # reading data from path (csv_file)
        data_in_glob = pd.read_csv(csv_file)
        
        print(csv_file)
        
        time_len = len(data_in_glob["Time"])
        
        # creating new dataFrame to store the data in required format
        
        df = pd.DataFrame(index=np.arange(800), columns=np.arange(2))
        
        # filling the df up using csv_file
        df.at[0,0] = "Station"
        df.at[0,1] = csv_info["Name"]
        
        df.at[1,0] = "Latitude"
        df.at[1,1] = float(csv_info["lat"])
        
        df.at[2,0] = "Longitude"
        df.at[2,1] = float(csv_info["lon"])
        
        df.at[3,0] = "Elevation"
        df.at[3,1] = int(csv_info["elev"])
        
        df.at[5,0] = "Time"
        df.at[5,1] = varname
        
        #adding data from sorted file
        for i in range(time_len):
            df.loc[6+i,0] = data_in_glob["Time"][i]
            df.loc[6+i,1] = data_in_glob[varname][i]
        
        # saving the file with station name
        
        name = csv_info["Name"]
        
        # replacing special characters eg. German umlaute
        special_char_map = {ord('ä'):'ae', ord('ü'):'ue', ord('ö'):'oe', ord('ß'):'ss', 
                            ord(",") : "", ord(" ") : "_", ord("("):"", ord("."):"",
                            ord(")"):""}
        name = name.translate(special_char_map)
        
        df.to_csv(os.path.join(path_to_store, name + ".csv"), index=False, header=False)
        
    df_info = df_info.rename(columns = {"lon":"Longitude",
                                        "lat":"Latitude","elev":"Elevation"})
    
    # replace comma with dot
    
    df_info["Name"] = df_info["Name"].apply(lambda x:x.translate(special_char_map))
    
    df_info = df_info.sort_values(by=["Name"], ascending=True)
    
    df_info = df_info.reset_index(drop=True)
                                 
    
    df_info.to_csv(os.path.join(path_to_store, "stationnames.csv"))
    
    df_info.drop("Name", axis=1, inplace=True)
    
    df_info.to_csv(os.path.join(path_to_store, "stationloc.csv"))


path_to_iso = "D:/Datasets/GNIP_data/processed/Temperature"
path_processed = "D:/Datasets/GNIP_data/processed/Temperature/processed"
path_info = "D:/Datasets/GNIP_data/station_info.csv"

#process_data_for_variable(path_to_data=path_to_processed, path_to_store="D:/Datasets/GNIP_data/processed/Temperature")

# add_info_to_data(path_to_info=path_info, path_to_data=path_to_iso, path_to_store=path_processed, 
#                   glob_name="*.csv", 
#                       varname="Temperature")


preprocess_to_long_records()