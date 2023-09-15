# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 14:53:06 2023

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
from pyClimat.utils import extract_region
from pyClimat.variables import extract_var


from paths_to_data import *
echam_pd_data_path = "D:/Datasets/Model_output_pst/PD"

def perform_causality_winter_to_summer(Y_NAO_DJF=None, Z_EA_DJF=None, X_varname=None, X_units=None, save_pval=True,
                                       path_to_save=main_path_to_data, filename=None, season = "JJA"):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    
    data_x = extract_var(Dataset=echam_data, varname=X_varname, units=X_units, Dataset_wiso=echam_wiso)
    
    data_x_season = extract_region(data=data_x, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season=season)
    
        
    Granger_object = GrangerCausality(maxlag=1, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(Y=Y_NAO_DJF, X=data_x_season, Z=Z_EA_DJF, 
                                                       apply_standardize=True)
    if save_pval:
        pval.to_netcdf(os.path.join(path_to_save, filename + ".nc"))
    else:
        return pval
 
    
 
def perform_NAO_to_climate():
    
    # read all the required datsets (for winter)
    df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

    #select the node ands convert to xarray
    nao_index_echam = xr.DataArray(df_echam_pcs[str(1)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})

    
    perform_causality_winter_to_summer(Y_NAO_DJF=nao_index_echam, X_varname="d18op", X_units="per mil", save_pval=True,
                            path_to_save=main_path_to_data, filename="d18op_JJA_caused_by_NAO_DJF", 
                            season = "JJA", )
    
    perform_causality_winter_to_summer(Y_NAO_DJF=nao_index_echam, X_varname="prec", X_units="mm/month", save_pval=True,
                            path_to_save=main_path_to_data, filename="prec_JJA_caused_by_NAO_DJF", 
                            season = "JJA", )
    
    perform_causality_winter_to_summer(Y_NAO_DJF=nao_index_echam, X_varname="temp2", X_units="°C", save_pval=True,
                            path_to_save=main_path_to_data, filename="temp_JJA_caused_by_NAO_DJF", 
                            season = "JJA",)
    
    
    
    
def perform_NAO_EA_to_climate():
    
    # read all the required datsets (for winter)
    df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

    #select the node ands convert to xarray
    nao_index_echam = xr.DataArray(df_echam_pcs[str(1)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})
    
    ea_index_echam = xr.DataArray(df_echam_pcs[str(2)], dims="time", coords={"time": df_echam_pcs["time"]})

    
    perform_causality_winter_to_summer(Y_NAO_DJF=nao_index_echam, X_varname="d18op", X_units="per mil", save_pval=True,
                            path_to_save=main_path_to_data, filename="d18op_JJA_caused_by_NAO_EA_DJF", 
                            season = "JJA", Z_EA_DJF=ea_index_echam)
    
    perform_causality_winter_to_summer(Y_NAO_DJF=nao_index_echam, X_varname="prec", X_units="mm/month", save_pval=True,
                            path_to_save=main_path_to_data, filename="prec_JJA_caused_by_NAO_EA_DJF", 
                            season = "JJA", Z_EA_DJF=ea_index_echam)
    
    perform_causality_winter_to_summer(Y_NAO_DJF=nao_index_echam, X_varname="temp2", X_units="°C", save_pval=True,
                            path_to_save=main_path_to_data, filename="temp_JJA_caused_by_NAO_EA_DJF", 
                            season = "JJA", Z_EA_DJF=ea_index_echam)
    

if __name__ == "__main__":
    #perform_NAO_to_climate()
    perform_NAO_EA_to_climate()