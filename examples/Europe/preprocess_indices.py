# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 11:23:11 2023

@author: dboateng
"""

import pandas as pd 
import numpy as np
import os 

#reading data
main_path = "C:/Users/dboateng/Desktop/Datasets/NAO"
data_nao_cdc = np.loadtxt(os.path.join(main_path, "cdc_nao_monthly.ascii")) #from nao
data_nao_cru = np.loadtxt(os.path.join(main_path, "nao.dat")) # Gibratar
data_nao_azore = np.loadtxt(os.path.join(main_path, "nao_azo.dat"), skiprows=1) # Gibratar
data_ea_valencia = pd.read_csv(os.path.join(main_path, "WMO_003953_SLP_anomaly.tab"), sep='\t', skiprows=12, 
                               index_col=["Date/Time"], parse_dates=["Date/Time"]) #from valencia


# organize into time and data (with datetime and extract seasons)
time = pd.date_range(start="1821-01-01", end="2022-12-01", freq="MS")
time_cdc = pd.date_range(start="1950-01-01", end="2023-01-01", freq="MS")
time_ea = pd.date_range(start="1866-01-01", end="2016-12-01", freq="MS")
time_azore = pd.date_range(start="1865-01-01", end="2002-12-01", freq="MS")

df_nao_gib = pd.DataFrame(index=time, columns=["NOA"])
df_nao_azore = pd.DataFrame(index=time_azore, columns=["NOA"])
df_nao_cdc = pd.DataFrame(index=time_cdc, columns=["NAO"])
df_ea_val = data_ea_valencia.loc[time_ea[0]:time_ea[-1]]

df_nao_cdc["NAO"] =data_nao_cdc[:, 2]


# use loop sort data (from cru)
a = 0
for i in range(0, len(df_nao_gib), 12):
        
    df_nao_gib.iloc[i:i+12, 0] = data_nao_cru[a,1:13].T #print(i, a)
    a = a+1

df_nao_gib[df_nao_gib ==-99.99] = np.nan

b = 0
for i in range(0, len(df_nao_azore), 12):
        
    df_nao_azore.iloc[i:i+12, 0] = data_nao_azore[b,1:13].T #print(i, a)
    b = b+1

df_nao_azore[df_nao_azore ==-10] = np.nan

data_to_save = [df_nao_gib, df_nao_cdc, df_nao_azore, df_ea_val]
filenames = ["NAO_Gilbraltar", "NAO_CDC", "NAO_Azore", "EA_Valencia"]

for i,df in enumerate(data_to_save):
    season = ((df.index.month % 12 + 3) // 3).map({1:'DJF', 2: 'MAM', 3:'JJA', 4:'SON'})
    df_DJF = df[season == "DJF"].to_csv(os.path.join(main_path, filenames[i] + "_DJF.csv"))
    df_MAM = df[season == "MAM"].to_csv(os.path.join(main_path,filenames[i] + "_MAM.csv"))
    df_JJA = df[season == "JJA"].to_csv(os.path.join(main_path,filenames[i] + "_JJA.csv"))
    df_SON = df[season == "SON"].to_csv(os.path.join(main_path,filenames[i] + "_SON.csv"))
    df.to_csv(os.path.join(main_path, filenames[i] + ".csv"))
    