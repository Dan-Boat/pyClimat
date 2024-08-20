# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 11:59:44 2024

@author: dboateng
"""

# import modules
import pandas as pd
import numpy as np
import os
import pyleoclim as pyleo
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import patches
import xarray as xr
import pyreadr

from pyClimat.plot_utils import apply_style

path_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots/data"
path_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"
ortega = os.path.join(path_data, "ortega2015nao.txt")
trouet = os.path.join(path_data, "nao-trouet2009-noaa.txt")
cesm = os.path.join(path_data, "CESM_DJF_pcs.csv")
giss = os.path.join(path_data, "GISS_DJF_pcs.csv")
echam = os.path.join(path_data, "ECHAM5_DJF_pcs.csv")


def shift_december(ds):
    time = ds['time']
    # Create a new time coordinate where December's year is shifted by +1
    new_years = xr.where(time.dt.month == 12, time.dt.year + 1, time.dt.year)
    return new_years


def nao_index_from_models(path, factor=-1, model="CESM"):
    
    df = pd.read_csv(path, usecols=["1", "time"])
    df = df.rename(columns={"1": "NAO"})
    
    if model == "GISS" or model=="ECHAM":
        df = df[:-1]
        
        ds = xr.DataArray(df["NAO"]*factor, dims="time", coords={"time": df["time"]})
        
        dates =xr.cftime_range(start=ds.time.values[0], end=ds.time.values[-1],
                               freq="M", calendar="noleap")
        
    else:
            
        ds = xr.DataArray(df["NAO"]*factor, dims="time", coords={"time": df["time"]})
        
        dates =xr.cftime_range(start=ds.time.values[0], end=ds.time.values[-1],
                               freq="MS", calendar="noleap")
        
    
    data = xr.DataArray(np.arange(1,len(dates)+1), dims="time", coords={"time":dates})
    
    winter_months = data.sel(time=data['time.month'].isin([12, 1, 2]))
    
    if model=="HADCM3":
        winter_months = winter_months[:-13]
    
    ds["time"] = winter_months.time
    
    shifted_years = shift_december(ds)
    
    nao = ds.groupby(shifted_years).mean(dim="time").to_dataframe()
    
    
    df_series = pyleo.Series(
        time=nao.index,
        value=nao["NAO"],
        time_name="age_AD",
        time_unit="years AD",
        value_name="NAO",
        value_unit="-",
        label=model,
    )
    
    return df_series

def pause():

    # Read the file line by line and filter out lines starting with '#'
    lines = []
    with open(ortega, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                lines.append(line)
    
    # Use pandas to read the filtered lines as a dataframe
    df1 = pd.DataFrame([line.strip().split() for line in lines])
    
    df1.columns = df1.loc[0]
    df1 = df1[1:]
    df1['age_AD'] = df1['age_AD'].astype('int')
    df1['nao_mc_mean'] = df1['nao_mc_mean'].astype('float')
    df1['nao_mc_min'] = df1['nao_mc_min'].astype('float')
    df1['nao_mc_max'] = df1['nao_mc_max'].astype('float')
    df1['nao_mc_mean-2SE'] = df1['nao_mc_mean-2SE'].astype('float')
    df1['nao_mc_mean+2SE'] = df1['nao_mc_mean+2SE'].astype('float')
    
    
    df1_series = pyleo.Series(
        time=df1["age_AD"],
        value=df1["nao_mc_mean"],
        time_name="age_AD",
        time_unit="years AD",
        value_name="NAO",
        value_unit="-",
        label="Ortega et al. 2015",
    )
    
    # Read the file line by line and filter out lines starting with '#'
    lines = []
    with open(trouet, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                lines.append(line)
    
    # Use pandas to read the filtered lines as a dataframe
    df2 = pd.DataFrame([line.strip().split() for line in lines])
    
    df2.columns = df2.loc[0]
    df2 = df2[1:]
    df2['age'] = df2['age'].astype('int')
    df2['nao'] = df2['nao'].astype('float')
    
    df2_series = pyleo.Series(
        time=df2["age"],
        value=df2["nao"],
        time_name="age_AD",
        time_unit="years AD",
        value_name="NAO",
        value_unit="-",
        label="Trouet et al. 2009",
    )
    
    cesm_series = nao_index_from_models(cesm, model="CESM")
    giss_series = nao_index_from_models(giss, model="GISS", factor=1)
    echam_series = nao_index_from_models(echam, model="ECHAM", factor=-1)
       
    
    nao_series = [df2_series.standardize(), df1_series.standardize().resample("15y").mean(), 
                  cesm_series.standardize().resample("15y").mean(),
                  giss_series.standardize().resample("15y").mean(),
                  echam_series.standardize().resample("15y").mean()
                  ]
    
    
    
    nao = pyleo.MultipleSeries(nao_series, name="North Atlantic Oscillation Comparison")
    
    
    apply_style(fontsize=28, style="seaborn-paper", linewidth=3,)
    
    
    fig, ax = nao.stackplot(time_unit="year AD", xlim=[850, 2000], colors=["b", "g", "r", "#42eff5", "#a911ba"], figsize=(15,12),
                           ylabel_fontsize=25)
    for i in range(5):
        ax.get(i).axvspan(950, 1300, facecolor="grey", alpha=0.2)
        ax.get(i).axvspan(1650, 1850, facecolor="grey", alpha=0.2)
    
    plt.savefig(os.path.join(path_plots, "nao_comparison.pdf"), bbox_inches="tight", format= "pdf")
    
pause()