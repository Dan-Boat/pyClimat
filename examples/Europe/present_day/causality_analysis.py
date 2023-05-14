# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 20:11:32 2023

@author: dboateng
1. perform granger causality testing between the NAO, EA, SCAND and climate variable
using echam and ERA5
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
from pyClimat.analysis import extract_var, extract_transect
from pyClimat.utils import extract_region


from paths_to_data import *
echam_pd_data_path = "D:/Datasets/Model_output_pst/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021"
        
        
# read all the required datsets (for winter)
df_era_pcs = pd.read_csv(PCs_ERA_DJF_path , parse_dates=["time"])
df_echam_pcs = pd.read_csv(PCs_ECHAM_DJF_path, parse_dates=["time"])

#select the node ands convert to xarray

nao_index_era = xr.DataArray(df_era_pcs[str(1)] * -1, dims="time", coords={"time": df_era_pcs["time"]})
nao_index_echam = xr.DataArray(df_echam_pcs[str(1)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})


ea_index_era = xr.DataArray(df_era_pcs[str(3)], dims="time", coords={"time": df_era_pcs["time"]})
ea_index_echam = xr.DataArray(df_echam_pcs[str(2)], dims="time", coords={"time": df_echam_pcs["time"]})

scan_index_era = xr.DataArray(df_era_pcs[str(2)], dims="time", coords={"time": df_era_pcs["time"]})
scan_index_echam = xr.DataArray(df_echam_pcs[str(3)] *-1, dims="time", coords={"time": df_echam_pcs["time"]})

def perform_causality_indices(filename="causal_results.csv", plot=True):
    # perform causality with indices (echam and era)
    column_names = ["NAO_by_EA", "NAO_by_EA_SCAN", "EA_by_NAO", "EA_by_NAO_SCAN",
                    "NAO_by_EA_era", "NAO_by_EA_SCAN_era", "EA_by_NAO_era", "EA_by_NAO_SCAN_era"]
    
    df_results = pd.DataFrame(columns=column_names, index=np.arange(1))
    
    G = GrangerCausality(maxlag=15, test="params_ftest")
    df_results["NAO_by_EA"] = G.perform_granger_test(X=nao_index_echam, Y=ea_index_echam,  
                                                      apply_standardize=True, dim="time")
    
    df_results["NAO_by_EA_SCAN"] = G.perform_granger_test(X=nao_index_echam, Y=ea_index_echam, Z=scan_index_echam, 
                                                      apply_standardize=True, dim="time")
    
    df_results["EA_by_NAO"] = G.perform_granger_test(Y=nao_index_echam, X=ea_index_echam,  
                                                      apply_standardize=True, dim="time")
    
    df_results["EA_by_NAO_SCAN"] = G.perform_granger_test(Y=nao_index_echam, X=ea_index_echam, Z=scan_index_echam, 
                                                      apply_standardize=True, dim="time")
    
    
    df_results["NAO_by_EA_era"] = G.perform_granger_test(X=nao_index_era, Y=ea_index_era,  
                                                      apply_standardize=True, dim="time")
    
    df_results["NAO_by_EA_SCAN_era"] = G.perform_granger_test(X=nao_index_era, Y=ea_index_era, Z=scan_index_era, 
                                                      apply_standardize=True, dim="time")
    
    df_results["EA_by_NAO_era"] = G.perform_granger_test(Y=nao_index_era, X=ea_index_era,  
                                                      apply_standardize=True, dim="time")
    
    df_results["EA_by_NAO_SCAN_era"] = G.perform_granger_test(Y=nao_index_era, X=ea_index_era, Z=scan_index_era, 
                                                      apply_standardize=True, dim="time")
    
    if plot:
        apply_style(fontsize=23, style="seaborn-talk", linewidth=3,)
        fig,ax = plt.subplots(1,1, sharex=False, figsize=(18, 15))
        df_results.iloc[0].plot(kind="barh", ax=ax)
        ax.axvline(x=0.1, linestyle="-", color="red")
        ax.axvline(x=0.33, linestyle="-", color="blue")
        ax.axvline(x=0.66, linestyle="-", color="magenta")
        plt.savefig(os.path.join(path_to_plots, "causal_analysis_indices.svg"), format= "svg", bbox_inches="tight", dpi=300)
        
    df_results.to_csv(os.path.join(path_to_plots, filename))



def perform_causality_testing_with_climate_var(Y, varname, units=None, lev=None, lev_units=None, save_pval=True,
                        path_to_save=main_path_to_data, filename=None, Z=None, revrese=True, season="DJF"):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data = extract_var(Dataset=echam_data, varname=varname, units=units, Dataset_wiso=echam_wiso,
                       lev=lev, lev_units=lev_units)
    
    data_season = extract_region(data=data, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season=season)
    
    Granger_object = GrangerCausality(maxlag=15, test="params_ftest")
    
    pval = Granger_object.perform_granger_test(X=data_season, Y=Y, Z=Z, apply_standardize=True)
    
        
    
    if save_pval:

        pval.to_netcdf(os.path.join(path_to_save, filename + ".nc"))
        
    
    else:
        return pval

def experiment1():
    # NAO
    pval_nao_d18op, pval_d18op_nao = perform_causality_testing_with_climate_var(Y=nao_index_echam, 
                                            varname="d18op", units="per mil", filename="NAO_caused_by_d18op")
    
    pval_nao_t2m, pval_t2m_nao = perform_causality_testing_with_climate_var(Y=nao_index_echam, 
                                            varname="temp2", units="°C", filename="NAO_caused_by_t2m")
    
    pval_nao_prec, pval_prec_nao = perform_causality_testing_with_climate_var(Y=nao_index_echam, 
                                            varname="prec", units="mm/month", filename="NAO_caused_by_prec")
    
    
    # EA
    # pval_ea_d18op, pval_d18op_ea = perform_causality_testing_with_climate_var(Y=ea_index_echam, 
    #                                         varname="d18op", units="per mil", filename="EA_caused_by_d18op")
    
    # pval_ea_t2m, pval_t2m_ea = perform_causality_testing_with_climate_var(Y=ea_index_echam, 
    #                                         varname="temp2", units="°C", filename="EA_caused_by_t2m")
    
    # pval_ea_prec, pval_prec_ea = perform_causality_testing_with_climate_var(Y=ea_index_echam, 
    #                                         varname="prec", units="mm/month", filename="EA_caused_by_prec")
    
    
    pval_d18op_nao_ea = perform_causality_testing_with_climate_var(Y=nao_index_echam, Z=ea_index_echam,
                                                                    varname="d18op", units="per mil", 
                                                                    filename="d18op_caused_by_NAO_EA",
                                                                    revrese=False)
    
    pval_t2m_nao_ea = perform_causality_testing_with_climate_var(Y=nao_index_echam, Z=ea_index_echam,
                                                                    varname="temp2", units="°C", 
                                                                    filename="t2m_caused_by_NAO_EA",
                                                                    revrese=False)
    
    pval_prec_nao_ea = perform_causality_testing_with_climate_var(Y=nao_index_echam, Z=ea_index_echam,
                                                                    varname="prec", units="mm/month", 
                                                                    filename="prec_caused_by_NAO_EA",
                                                                    revrese=False)



def perform_causality_testing_with_d18op(Y_varname, Y_units=None, Z_units=None, save_pval=True,
                        path_to_save=main_path_to_data, filename=None, Z_varname="NAO",):
    
    echam_data = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly.nc", decode=True)   
    echam_wiso = read_from_path(echam_pd_data_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    
    data_d18op = extract_var(Dataset=echam_data, varname="d18op", units="per mil", Dataset_wiso=echam_wiso)
    
    d18op_season = extract_region(data=data_d18op, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF")
    
    data_Y = extract_var(Dataset=echam_data, varname=Y_varname, units=Y_units, Dataset_wiso=echam_wiso,
                       )
    
    Y_season = extract_region(data=data_Y, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                 season="DJF")
    
    if Z_varname=="NAO":
        Z = nao_index_echam
        
    elif Z_varname ==None:
        Z = None
        
    else:
        data_Z = extract_var(Dataset=echam_data, varname=Z_varname, units=Z_units, Dataset_wiso=echam_wiso,
                           )
        
        Z = extract_region(data=data_Z, maxlon=40, minlon=-40, maxlat=80, minlat=30, time="season", 
                                     season="DJF")
        
    
    Granger_object = GrangerCausality(maxlag=15, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(X=d18op_season, Y=Y_season, Z=Z, 
                                                       apply_standardize=True)
    if save_pval:
        pval.to_netcdf(os.path.join(path_to_save, filename + ".nc"))
    else:
        return pval

# pval_prec_t2m_to_d18op = perform_causality_testing_with_d18op(Y_varname="temp2", Y_units="°C", Z_varname="prec", 
#                                                             Z_units="mm/month",
#                                                            filename="d18op_caused_by_t2m_prec")


# pval_nao_t2m_to_d18op = perform_causality_testing_with_d18op(Y_varname="temp2", Y_units="°C", Z_varname="NAO",
#                                                            filename="d18op_caused_by_t2m_NAO")


# pval_t2m_to_d18op = perform_causality_testing_with_d18op(Y_varname="temp2", Y_units="°C", Z_varname=None,
#                                                            filename="d18op_caused_by_t2m")

# pval_prec_to_d18op = perform_causality_testing_with_d18op(Y_varname="prec", Y_units="mm/month", Z_varname=None,
#                                                            filename="d18op_caused_by_prec")

if __name__ == "__main__":
    perform_causality_indices()