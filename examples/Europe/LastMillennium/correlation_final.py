# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:38:25 2024

@author: dboateng
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates

from pyClimat.plot_utils import *
from pyClimat.plots import plot_eofsAsCovariance, plot_correlation
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_transect
from pyClimat.variables import extract_var
from pyClimat.utils import extract_region


from path_to_data_lm import *

main_path = "D:/Datasets/iGCM_datasets/"
path_to_save = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots/data"


def get_indices(df, nao_mode, ea_mode, nao_factor, ea_factor):
    
    
    # extract the indices for OP and EG
    
    df_eq_index = (df[str(nao_mode)] * nao_factor > 0) & (df[str(ea_mode)] * ea_factor > 0) | (df[str(nao_mode)] * nao_factor < 0) & (df[str(ea_mode)] * ea_factor < 0)
    df_op_index = (df[str(nao_mode)] * nao_factor > 0) & (df[str(ea_mode)] * ea_factor < 0) | (df[str(nao_mode)] * nao_factor < 0) & (df[str(ea_mode)] * ea_factor > 0)
    
    # extract them from the pcs 
    
    df_eq = df[df_eq_index]
    df_op = df[df_op_index]
    
    #create xarray for the pcs
    
    EQ_indices = xr.DataArray(df_eq[str(nao_mode)] * nao_factor, dims="time", coords={"time": df_eq["time"]})
    OP_indices = xr.DataArray(df_op[str(nao_mode)] * nao_factor, dims="time", coords={"time": df_op["time"]})
    
    return EQ_indices, OP_indices


def perform_correlation_composite(atmos_index, df, t2m_season, prec_season, d18O_season, sign="OP"):   
    
    
        
    t2m_season["time"] = df.time.values
    
    t2m_index = t2m_season.sel(time = t2m_season.time.isin(atmos_index.time.values))
    
    prec_season["time"] = df.time.values
    
    prec_index = prec_season.sel(time = prec_season.time.isin(atmos_index.time.values))
    
    d18O_season["time"] = df.time.values
    
    d18O_index = d18O_season.sel(time = d18O_season.time.isin(atmos_index.time.values))
    
    
    
    nao_t2m_sval, nao_t2m_pval, nao_t2m_sig = StatCorr(x=t2m_index, y=atmos_index, dim="time",
                                                       return_sig=True, sig=0.05)
    
 
    nao_prec_sval, nao_prec_pval, nao_prec_sig = StatCorr(x=prec_index, y=atmos_index, dim="time",
                                                       return_sig=True, sig=0.05)
    
    nao_d18O_sval, nao_d18O_pval, nao_d18O_sig = StatCorr(x=d18O_index, y=atmos_index, dim="time",
                                                       return_sig=True, sig=0.05)
    
    results_dataset = xr.Dataset({"nao_reg_temp_" + sign :nao_t2m_sval,
                                  "nao_reg_prec_" + sign :nao_prec_sval, 
                                  "nao_reg_d18O_" + sign :nao_d18O_sval,
                                  "nao_sig_temp_" + sign :nao_t2m_sig, 
                                  "nao_sig_prec_" + sign :nao_prec_sig,
                                  "nao_sig_d18O_" + sign :nao_d18O_sig})
    
    return results_dataset
    
def correlation_all(NAO_indices, t2m_season, prec_season, d18O_season):
    
    nao_t2m_sval, nao_t2m_pval, nao_t2m_sig = StatCorr(x=t2m_season, y=NAO_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    
    
    nao_prec_sval, nao_prec_pval, nao_prec_sig = StatCorr(x=prec_season, y=NAO_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    
    
    nao_d18O_sval, nao_d18O_pval, nao_d18O_sig = StatCorr(x=d18O_season, y=NAO_indices, dim="time",
                                                       return_sig=True, sig=0.05)
    
    
    results_dataset = xr.Dataset({"nao_reg_temp":nao_t2m_sval, "nao_reg_prec":nao_prec_sval, "nao_reg_d18O":nao_d18O_sval,
                                  "nao_sig_temp":nao_t2m_sig, "nao_sig_prec":nao_prec_sig, "nao_sig_d18O":nao_d18O_sig})
    
    return results_dataset
    

def perform_correlation_nao_climate(pcs_path, path_to_data, nao_mode, nao_factor=-1, filename=None,
                                             select_time=False, time_range=None, use_slice=False, 
                                             slice_start=None, slice_end=None, slice_end_plus=None,
                                             save_to_nc=True, path_to_save=None, ea_mode=None, 
                                             ea_factor=None):
    
    
    df = pd.read_csv(pcs_path, parse_dates=["time"])
    
    NAO_indices = xr.DataArray(df[str(nao_mode)] * nao_factor, dims="time", coords={"time": df["time"]})
    
    EQ_indices, OP_indices = get_indices(df, nao_mode, ea_mode, nao_factor, ea_factor)
    
    
    if select_time:
        NAO_indices = NAO_indices.sel(time=slice(slice_start, slice_end_plus))
        
        
    filename_data = filename + "_wiso_vars.nc"
    
    temp = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="tsurf")  - 273.15 #°C
    
    d18O = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="d18O")  
    prec = read_from_path(path=path_to_data, filename=filename_data, decode=True,
                          varname="prec")
    
    
    
    maxlon = 150
    minlon = -140
    maxlat = 86
    minlat = 10
    
    t2m_season = extract_region(data=temp, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF", select_time=select_time, time_range=time_range, use_slice=use_slice, 
                                 slice_start=slice_start, slice_end=slice_end) 
    
    prec_season = extract_region(data=prec, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF", select_time=select_time, time_range=time_range, use_slice=use_slice, 
                                 slice_start=slice_start, slice_end=slice_end) 
    
    d18O_season = extract_region(data=d18O, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, time="season", 
                                 season="DJF", select_time=select_time, time_range=time_range, use_slice=use_slice, 
                                 slice_start=slice_start, slice_end=slice_end)
    
    
    
    results_dataset_all = correlation_all(NAO_indices, t2m_season, prec_season, d18O_season)
    
    results_dataset_OP = perform_correlation_composite(OP_indices, df, t2m_season, 
                                                          prec_season, d18O_season, sign="OP")
    
    results_dataset_EQ = perform_correlation_composite(EQ_indices, df, t2m_season, 
                                                          prec_season, d18O_season, sign="EQ")
    
    results = xr.merge([results_dataset_all, results_dataset_OP, results_dataset_EQ])
    
    
    if save_to_nc:
        
        results.to_netcdf(os.path.join(path_to_save, filename + "_corr_data.nc"))
        
    
    return results


cesm_data  = perform_correlation_nao_climate(pcs_path= CESM_DJF_PCs, path_to_data=main_path, nao_mode=1, 
                                                 nao_factor=-1, filename="CESM", save_to_nc=True,
                                                 path_to_save=path_to_save, ea_mode=3, ea_factor=1)

giss_data  = perform_correlation_nao_climate(pcs_path= GISS_DJF_PCs, path_to_data=main_path, nao_mode=1, 
                                                  nao_factor=-1, filename="GISS", save_to_nc=True,
                                                  path_to_save=path_to_save, ea_mode=3, ea_factor=1)