# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 20:12:22 2023

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
from pyClimat.plots import plot_eofsAsCovariance
from pyClimat.stats import sliding_correlation, StatCorr, GrangerCausality
from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.analysis import extract_transect
from pyClimat.variables import extract_var
from pyClimat.utils import extract_region


from paths_to_data import *
echam_pd_data_path = "D:/Datasets/Model_output_pst/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021"
        
        



def perform_causality_testing_to_NAO(Y_varname, Y_units, EA_index, NAO_index, save_pval=True,
                        path_to_save=main_path_to_data, filename=None, season = "DJF"):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data_y = extract_var(Dataset=echam_data, varname=Y_varname, units=Y_units, Dataset_wiso=echam_wiso)
    
    data_y_season = extract_region(data=data_y, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season=season)
        
    
    Granger_object = GrangerCausality(maxlag=18, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(X=NAO_index, Y=data_y_season, Z=EA_index, 
                                                       apply_standardize=True)
    if save_pval:
        pval.to_netcdf(os.path.join(path_to_save, filename + ".nc"))
    else:
        return pval
    

def perform_causality_testing_to_climate(X_varname, X_units, EA_index, NAO_index, save_pval=True,
                        path_to_save=main_path_to_data, filename=None, season = "DJF"):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data_x = extract_var(Dataset=echam_data, varname=X_varname, units=X_units, Dataset_wiso=echam_wiso)
    
    data_x_season = extract_region(data=data_x, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season=season)
        
    
    Granger_object = GrangerCausality(maxlag=18, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(Y=NAO_index, X=data_x_season, Z=EA_index, 
                                                       apply_standardize=True)
    if save_pval:
        pval.to_netcdf(os.path.join(path_to_save, filename + ".nc"))
    else:
        return pval
    
    
    
def perform_DJF_to_NAO():
    
    # read all the required datsets (for winter)
    df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

    #select the node ands convert to xarray
    nao_index_echam = xr.DataArray(df_echam_pcs[str(1)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})

    ea_index_echam = xr.DataArray(df_echam_pcs[str(2)], dims="time", coords={"time": df_echam_pcs["time"]})
    
    perform_causality_testing_to_NAO(Y_varname="d18op", Y_units="per mil", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="NAO_caused_by_d18op_EA_DJF", season = "DJF")
    
    perform_causality_testing_to_NAO(Y_varname="temp2", Y_units="°C", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="NAO_caused_by_t2m_EA_DJF", season = "DJF")
    
    perform_causality_testing_to_NAO(Y_varname="prec", Y_units="mm/month", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="NAO_caused_by_prec_EA_DJF", season = "DJF")
    
    
def perform_DJF_to_climate():
    
    # read all the required datsets (for winter)
    df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

    #select the node ands convert to xarray
    nao_index_echam = xr.DataArray(df_echam_pcs[str(1)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})

    ea_index_echam = xr.DataArray(df_echam_pcs[str(2)], dims="time", coords={"time": df_echam_pcs["time"]})
    
    perform_causality_testing_to_climate(X_varname="d18op", X_units="per mil", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="d18op_caused_by_NAO_EA_DJF", season = "DJF")
    
    perform_causality_testing_to_climate(X_varname="temp2", X_units="°C", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="t2m_caused_by_NAO_EA_DJF", season = "DJF")
    
    perform_causality_testing_to_climate(X_varname="prec", X_units="mm/month", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="prec_caused_by_NAO_EA_DJF", season = "DJF")
    
    
def perform_NAO_EA_DJF_to_climate_JJA():
    
    # read all the required datsets (for winter)
    df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

    #select the node ands convert to xarray
    nao_index_echam = xr.DataArray(df_echam_pcs[str(1)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})

    ea_index_echam = xr.DataArray(df_echam_pcs[str(2)], dims="time", coords={"time": df_echam_pcs["time"]})
    
    perform_causality_testing_to_climate(X_varname="d18op", X_units="per mil", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="d18op_JJA_caused_by_NAO_EA_DJF", season = "JJA")
    
    perform_causality_testing_to_climate(X_varname="temp2", X_units="°C", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="t2m_JJA_caused_by_NAO_EA_DJF", season = "JJA")
    
    perform_causality_testing_to_climate(X_varname="prec", X_units="mm/month", EA_index=ea_index_echam, 
                                     NAO_index=nao_index_echam, save_pval=True,
                            path_to_save=main_path_to_data, filename="prec_JJA_caused_by_NAO_EA_DJF", season = "JJA")
    
    
    
    
def perform_JJA_to_NAO():
    
    #summer 
    # read all the required datsets (for winter)
    df_echam_pcs_jja = pd.read_csv(PCs_ECHAM_JJA_path, parse_dates=["time"])

    #select the node ands convert to xarray

    nao_index_echam_jja = xr.DataArray(df_echam_pcs_jja[str(1)], dims="time", coords={"time": df_echam_pcs_jja["time"]})

    ea_index_echam_jja = xr.DataArray(df_echam_pcs_jja[str(2)] *-1, dims="time", coords={"time": df_echam_pcs_jja["time"]})
    
    perform_causality_testing_to_NAO(Y_varname="d18op", Y_units="per mil", EA_index=ea_index_echam_jja, 
                                     NAO_index=nao_index_echam_jja, save_pval=True,
                            path_to_save=main_path_to_data, filename="NAO_caused_by_d18op_EA_JJA", season = "JJA")
    
    perform_causality_testing_to_NAO(Y_varname="temp2", Y_units="°C", EA_index=ea_index_echam_jja, 
                                     NAO_index=nao_index_echam_jja, save_pval=True,
                            path_to_save=main_path_to_data, filename="NAO_caused_by_t2m_EA_JJA", season = "JJA")
    
    perform_causality_testing_to_NAO(Y_varname="prec", Y_units="mm/month", EA_index=ea_index_echam_jja, 
                                     NAO_index=nao_index_echam_jja, save_pval=True,
                            path_to_save=main_path_to_data, filename="NAO_caused_by_prec_EA_JJA", season = "JJA")
    
    
def perform_JJA_to_climate():
    
    #summer 
    # read all the required datsets (for winter)
    df_echam_pcs_jja = pd.read_csv(PCs_ECHAM_JJA_path, parse_dates=["time"])

    #select the node ands convert to xarray

    nao_index_echam_jja = xr.DataArray(df_echam_pcs_jja[str(1)], dims="time", coords={"time": df_echam_pcs_jja["time"]})

    ea_index_echam_jja = xr.DataArray(df_echam_pcs_jja[str(2)] *-1, dims="time", coords={"time": df_echam_pcs_jja["time"]})
    
    perform_causality_testing_to_climate(X_varname="d18op", X_units="per mil", EA_index=ea_index_echam_jja, 
                                     NAO_index=nao_index_echam_jja, save_pval=True,
                            path_to_save=main_path_to_data, filename="d18op_caused_by_NAO_EA_JJA", season = "JJA")
    
    perform_causality_testing_to_climate(X_varname="temp2", X_units="°C", EA_index=ea_index_echam_jja, 
                                     NAO_index=nao_index_echam_jja, save_pval=True,
                            path_to_save=main_path_to_data, filename="t2m_caused_by_NAO_EA_JJA", season = "JJA")
    
    perform_causality_testing_to_climate(X_varname="prec", X_units="mm/month", EA_index=ea_index_echam_jja, 
                                     NAO_index=nao_index_echam_jja, save_pval=True,
                            path_to_save=main_path_to_data, filename="prec_caused_by_NAO_EA_JJA", season = "JJA")
    
    
    
if __name__ == "__main__":
    #perform_DJF_to_NAO()
    #perform_JJA_to_NAO()
    #perform_DJF_to_climate()
    #perform_JJA_to_climate()
    perform_NAO_EA_DJF_to_climate_JJA()
    
    
    
    
    