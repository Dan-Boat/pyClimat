# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:19:46 2023

@author: dboateng

define path to the eofs (DJF and JJA)
"""
import os 

main_path_to_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/MH"

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

main_path_lgm = "D:/Datasets/CMIP6/PMIP/postprocessed/LGM"
from path_to_data_lgm import *

awi_path = os.path.join(main_path_lgm, "AWI-ESM-1-1-LR")
cesm_waccm_path = os.path.join(main_path_lgm, "CESM2-WACCM-FV2")
inm_cm_path = os.path.join(main_path_lgm, "INM-CM4-8")
miroc_path = os.path.join(main_path_lgm, "MIROC-ES2L")
mpi_esm_path = os.path.join(main_path_lgm, "MPI-ESM1-2-LR")

def define_region(region, path_to_data):
    if region=="IP":
        maxlon=3, 
        minlon=-11
        maxlat=44 
        minlat=36
    elif region == "FR":
        maxlon=6.5
        minlon=-5.5 
        maxlat=50 
        minlat=44
    elif region == "ALPS":
        maxlon=17
        minlon=5
        maxlat=48
        minlat=45
    elif region == "MED":
        maxlon=26
        minlon=2.5
        maxlat=44
        minlat=36
    elif region == "EE":
        maxlon=32
        minlon=18
        maxlat=55
        minlat=46
    elif region == "SCAND":
        maxlon=34
        minlon=1
        maxlat=67
        minlat=55
    elif region == "CE":
        maxlon=16
        minlon=3
        maxlat=54
        minlat=48
    elif region == "BI":
        maxlon=2
        minlon=-14
        maxlat=59
        minlat=51
    elif region == "ICE":
        maxlon=-11
        minlon=-30
        maxlat=67
        minlat=61
        
    t2m_data = read_from_path(path=path_to_data, filename="tas_monthly.nc", decode=True,varname="tas",
                              ) - 273.15 #°C
    prec_data = read_from_path(path=path_to_data, filename="pr_monthly.nc", decode=True, 
                               varname="pr",) *60*60*24*30  #mm/month
    
    t2m_season = extract_region(data=t2m_data, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, 
                                time="season", season="DJF", regional_mean=True) 
    
    prec_season = extract_region(data=prec_data, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat,
                                 time="season", season="DJF", regional_mean=True)
        
    return t2m_season, prec_season

def get_indices(pcs_path, nao_mode, ea_mode, nao_factor, ea_factor):
    
    
    df = pd.read_csv(pcs_path, parse_dates=["time"])
    
    NAO_indices = xr.DataArray(df[str(nao_mode)] * nao_factor, dims="time", coords={"time": df["time"]})
    EA_indices = xr.DataArray(df[str(ea_mode)] * ea_factor, dims="time", coords={"time": df["time"]})
    
    return NAO_indices, EA_indices


def perform_causal_climate_to_NAO(nao_index, path_to_data, region):
    
    
    t2m_season, prec_season = define_region(region, path_to_data) 
    
    
    Granger_object = GrangerCausality(maxlag=15, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(X=nao_index, Y=t2m_season, Z=prec_season, 
                                                       apply_standardize=True)
    
    return pval


def perform_causal_NAO_EA_to_t2m(nao_index, ea_index, path_to_data, region):
    
    t2m_season, prec_season = define_region(region, path_to_data) 
    
    
    Granger_object = GrangerCausality(maxlag=15, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(X=t2m_season, Y=nao_index, Z=ea_index, 
                                                       apply_standardize=True)
    
    return pval
 

def perform_causal_NAO_EA_to_prec(nao_index, ea_index, path_to_data, region):
    t2m_season, prec_season = define_region(region, path_to_data) 
    
    
    Granger_object = GrangerCausality(maxlag=15, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(X=prec_season, Y=nao_index, Z=ea_index, 
                                                       apply_standardize=True)
    
    return pval


def perform_causal_t2m_EA_to_NAO(nao_index, ea_index, path_to_data, region):
    
    
    t2m_season, prec_season = define_region(region, path_to_data) 
    
    
    Granger_object = GrangerCausality(maxlag=15, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(X=nao_index, Y=t2m_season, Z=ea_index, 
                                                       apply_standardize=True)
    
    return pval


def perform_causal_prec_EA_to_NAO(nao_index, ea_index, path_to_data, region):
    
    
    t2m_season, prec_season = define_region(region, path_to_data) 
    
    
    Granger_object = GrangerCausality(maxlag=15, test="params_ftest")
    
    
    
    pval = Granger_object.perform_granger_test(X=nao_index, Y=prec_season, Z=ea_index, 
                                                       apply_standardize=True)
    
    return pval




def perform_for_all(pcs_path, nao_mode, ea_mode, nao_factor, ea_factor, path_to_data):
    
    column_names = ["IP", "FR", "ALPS", "MED", "EE", "SCAND", "CE", "BI", "ICE"]
    
    df_to_nao = pd.DataFrame(columns=column_names, index=np.arange(1))
    df_to_t2m = pd.DataFrame(columns=column_names, index=np.arange(1))
    df_to_prec = pd.DataFrame(columns=column_names, index=np.arange(1))
    df_from_t2m = pd.DataFrame(columns=column_names, index=np.arange(1))
    df_from_prec = pd.DataFrame(columns=column_names, index=np.arange(1))
    
    NAO_indices, EA_indices = get_indices(pcs_path, nao_mode, ea_mode, nao_factor, ea_factor)
    
    for region in column_names:
        
        df_to_nao[region] = perform_causal_climate_to_NAO(NAO_indices, path_to_data, region)
        df_to_t2m[region] = perform_causal_NAO_EA_to_t2m(NAO_indices, EA_indices, path_to_data, region)
        df_to_prec[region] = perform_causal_NAO_EA_to_prec(NAO_indices, EA_indices, path_to_data, region)
        
        df_from_t2m[region] = perform_causal_climate_to_NAO(NAO_indices, path_to_data, region)
        df_from_prec[region] = perform_causal_NAO_EA_to_t2m(NAO_indices, EA_indices, path_to_data, region)
    
        
    results = {"climate_to_NAO":df_to_nao, "NAO_EA_to_t2m":df_to_t2m, "NAO_EA_to_prec":df_to_prec,
               "t2m_EA_to_NAO":df_from_t2m, "prec_EA_to_NAO":df_from_prec}
    
    return results

awi_data = perform_for_all(pcs_path=AWI_DJF_PCS, path_to_data=awi_path,
                                                    nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)

cesm_data = perform_for_all(pcs_path=CESM_WA_DJF_PCS, path_to_data=cesm_waccm_path, 
                                                     nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)


inm_data = perform_for_all(pcs_path=INM_DJF_PCS, path_to_data=inm_cm_path, 
                                                     nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)



miroc_data = perform_for_all(pcs_path=MIROC_DJF_PCS, path_to_data=miroc_path,
                                                    nao_mode=1, ea_mode=4, nao_factor=-1, ea_factor=1)

mpi_data = perform_for_all(pcs_path=MPI_ESM_DJF_PCS, path_to_data=mpi_esm_path, 
                                                    nao_mode=1, ea_mode=3, nao_factor=-1, ea_factor=1)




labels_lgm = ["AWI-ESM-1-1-LR", "CESM2-WACCM-FV2", "INM-CM4-8",
              "MIROC-ES2L", "MPI-ESM1-2-LR"]

lgm_data = [awi_data, cesm_data, inm_data, miroc_data, mpi_data]

column_names = ["IP", "FR", "ALPS", "MED", "EE", "SCAND", "CE", "BI", "ICE"]

df_to_nao = pd.DataFrame(columns=column_names, index=labels_lgm)

df_to_t2m = pd.DataFrame(columns=column_names, index=labels_lgm)

df_to_prec = pd.DataFrame(columns=column_names, index=labels_lgm)

df_from_t2m = pd.DataFrame(columns=column_names, index=labels_lgm)

df_from_prec = pd.DataFrame(columns=column_names, index=labels_lgm)


for i,label in enumerate(labels_lgm):
    
    df_to_nao.loc[label] = lgm_data[i].get("climate_to_NAO").values
    
    df_to_t2m.loc[label] = lgm_data[i].get("NAO_EA_to_t2m").values
    
    df_to_prec.loc[label] = lgm_data[i].get("NAO_EA_to_prec").values
    
    df_from_t2m.loc[label] = lgm_data[i].get("t2m_EA_to_NAO").values
    
    df_from_prec.loc[label] = lgm_data[i].get("prec_EA_to_NAO").values
    
df_to_nao.to_csv("LGM_climate_to_NAO.csv")
df_to_t2m.to_csv("LGM_NAO_EA_to_t2m.csv")
df_to_prec.to_csv("LGM_NAO_EA_to_prec.csv")
df_from_t2m.to_csv("LGM_t2m_EA_to_NAO.csv")
df_from_prec.to_csv("LGM_prec_EA_to_NAO.csv")


print(None)
    
    
    

