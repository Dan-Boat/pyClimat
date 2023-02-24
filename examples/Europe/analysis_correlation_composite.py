# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 12:21:20 2023

@author: dboateng
This script will experiment the implementation of the stats analysis betweeen
imported csv and xarray. 
"""
import os 
import numpy as np 
import pandas as pd 
import xarray as xr
import matplotlib.pyplot as plt 
from matplotlib.dates import YearLocator
import matplotlib.dates as mdates
from xarray import DataArray
from scipy import stats

from pyClimat.data import read_ERA_processed, read_from_path
from pyClimat.stats import _get_month, StatCorr
from pyClimat.utils import extract_region
from pyClimat.analysis import extract_var

path_to_save_files = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021"
path_to_other_vars_echam = "D:/Datasets/Model_output_pst/PD"


filename_era_pcs="ERA5_varimax_eof_pcs.csv"
filename_echam_pcs="ECHAM5-wiso_varimax_eof_pd_pcs.csv"


df_era_pcs = pd.read_csv(os.path.join(path_to_save_files, filename_era_pcs), parse_dates=["time"])
df_echam_pcs = pd.read_csv(os.path.join(path_to_save_files, filename_echam_pcs), parse_dates=["time"])

#select the node ands convert to xarray

nao_index = xr.DataArray(df_era_pcs[str(2)], dims="time", coords={"time": df_era_pcs["time"]})
ea_index = xr.DataArray(df_era_pcs[str(3)], dims="time", coords={"time": df_era_pcs["time"]})

nao_index_echam = xr.DataArray(df_echam_pcs[str(1)], dims="time", coords={"time": df_echam_pcs["time"]})
ea_index_echam = xr.DataArray(df_echam_pcs[str(4)], dims="time", coords={"time": df_echam_pcs["time"]})


# read era temperature
ERA5_t2m_path = os.path.join(ERA5_path, "t2m_monthly.nc")
ERA5_tp_path = os.path.join(ERA5_path, "tp_monthly.nc")

# ERA5_t2m = read_ERA_processed(path=ERA5_t2m_path, varname="t2m") - 273.15 #°C
# ERA5_tp = read_ERA_processed(path=ERA5_tp_path, varname="tp") * 1000 * 30  #mm/month

#extract the region and time

# data_t2m = extract_region(data=ERA5_t2m, time="season", season="DJF")
# data_tp = extract_region(data=ERA5_tp, time="season", season="DJF")

echam_data = read_from_path(path_to_other_vars_echam, "PD_1980_2014_monthly.nc", decode=True)   
echam_wiso = read_from_path(path_to_other_vars_echam, "PD_1980_2014_monthly_wiso.nc", decode=True)   
t2m_echam = extract_var(Dataset=echam_data, varname="temp2", units="°C")
tp_echam = extract_var(Dataset=echam_data, varname="prec", units="mm/month")
d18op_echam = extract_var(Dataset=echam_data, varname="d18op", units="per mil", Dataset_wiso=echam_wiso)
d18op = extract_region(data=d18op_echam, time="season", season="DJF")

# test for performing causality and stationarity testing

from statsmodels.tsa.stattools import adfuller, kpss, grangercausalitytests
from pyClimat.stats import StackArray
from sklearn.preprocessing import StandardScaler


#prepare data 
# X = StackArray(x=d18op, dim="time")
# Y = nao_index_echam.expand_dims(stacked=[0]).isel(stacked=0)

# # select one (but would require loop through all the stacks)

# X1 = X.isel(stacked=1)

# # standardize X

# X1 = X1.values.reshape(-1,1)
# scaler = StandardScaler()
# scaler = scaler.fit(X1)
# X1 = scaler.transform(X1)

# # prepare data for the statsmodel test

# lag=3
# data = np.column_stack([Y, X1])
# stats = grangercausalitytests(x=data, maxlag=lag)

# pvalues = [round(stats[i+1][0]['params_ftest'][1], 4) for i in range(lag)]
# pvalue = np.min(pvalues)
    
class GrangerCausality():
    """
    Note that this class is for spatial data at the monemt (station data would be implemented 
                                                            in a differnt class or adjust this
                                                            one)
    """
    
    def __init__(self, maxlag=3, test="params_ftest"):
        self.maxlag = maxlag
        self.test = test
        
    def prepare_data(self, X, Y, dim="time", ):
        
        if isinstance(Y, xr.DataArray):
            if len(Y.dims) ==1 or dim is None:
                sy = Y.expand_dims(stacked=[0])
            else:
                sy = StackArray(x=Y, dim=dim)
        
        if isinstance(X, xr.DataArray):
            
            if dim is None or len(X.dims) == 1:
                sx = X.expand_dims(stacked=[0])
            else:
                sx = StackArray(x=X,dim=dim)
        
        return sx, sy
    
    def compute_granger(self, i, x, y, apply_standardize=False, 
                             interchange=False):
        if x.shape == y.shape:
            xi = x.isel(stacked=i)
            yi = y.isel(stacked=i)
        else:
            xi = x.isel(stacked=i)
            yi = y.isel(stacked=0)
            
        if apply_standardize == True:
            xi = xi.values.reshape(-1,1)
            yi = yi.values.reshape(-1,1)
            
            scaler = StandardScaler()
            
            scaler_x = scaler.fit(xi)
            scaler_y = scaler.fit(yi)
            
            xi = scaler_x.transform(xi)
            yi = scaler_y.transform(yi)
            
        if interchange:
            data = np.column_stack([xi, yi])  # for checking the reverse order of causing
        else:
            data = np.column_stack([yi, xi])
            
        stats = grangercausalitytests(x=data, maxlag=self.maxlag, verbose=False)
        
        return stats   
        
            
            
            
                
    def perform_granger_test(self, X, Y, dim="time", apply_standardize=False, 
                             interchange=False):
        
        x,y = self.prepare_data(X, Y, dim)
        
        nspace = len(x.stacked)
        pval = np.zeros(x.stacked.shape)
        
        for i in range(nspace):
            stats = self.compute_granger(i, x, y, apply_standardize, interchange)
            
            pvalues = [round(stats[i+1][0][self.test][1], 4) for i in range(self.maxlag)]
            pvalue = np.min(pvalues)
            pval[i] = pvalue
            
        if nspace > 1:
            pvalx = DataArray(pval, coords=[x.stacked],name='pval').unstack('stacked')
            pvalx = pvalx.sortby("lon")
        else:
            pvalx = pval[0]
            
        return pvalx
            
            
from pyClimat.stats import StatCrossCorr 
 
coefs = StatCrossCorr(x=ea_index_echam, y=nao_index_echam, plot=True, sample_rate=1,
                      apply_standardize=True)         
        
sval,pval,sig = StatCorr(x=ea_index_echam, y=nao_index_echam, dim="time")      


G = GrangerCausality(maxlag=3, test="params_ftest")
pval = G.perform_granger_test(X=ea_index_echam, Y=nao_index_echam, apply_standardize=True)
pval_e = G.perform_granger_test(X=nao_index_echam, Y=ea_index_echam, apply_standardize=True)

era_pval = G.perform_granger_test(X=ea_index, Y=nao_index, apply_standardize=True)
era_pval_e = G.perform_granger_test(X=nao_index, Y=ea_index, apply_standardize=True)

print(pval)     
        
# sval,pval,sig = StatCorr(x=data_t2m, y=nao_index, dim="time")
# sval_tp,pval_tp,sig_tp = StatCorr(x=data_t2m, y=nao_index, dim="time")



