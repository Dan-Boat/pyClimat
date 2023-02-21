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

from pyClimat.data import read_ERA_processed
from pyClimat.stats import _get_month
from pyClimat.utils import extract_region

path_to_save_files = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/PD"
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021"
path_to_other_vars_echam = "D:/Datasets/Model_output_pst/PD"


filename_era_pcs="ERA5_standard_eof_pcs.csv"
filename_echam_pcs="ECHAM5-wiso_standard_eof_era_pcs.csv"


df_era_pcs = pd.read_csv(os.path.join(path_to_save_files, filename_era_pcs), parse_dates=["time"])

#select the node ands convert to xarray

nao_index = xr.DataArray(df_era_pcs[str(1)], dims="time", coords={"time": df_era_pcs["time"]})


# read era temperature
ERA5_t2m_path = os.path.join(ERA5_path, "t2m_monthly.nc")
ERA5_tp_path = os.path.join(ERA5_path, "tp_monthly.nc")

ERA5_t2m = read_ERA_processed(path=ERA5_t2m_path, varname="t2m") - 273.15 #°C
#ERA5_tp = read_ERA_processed(path=ERA5_tp_path, varname="tp") * 1000 * 30  #mm/month

#extract the region and time

data_t2m = extract_region(data=ERA5_t2m, time="season", season="DJF")





def StackArray(x,dim):
	'''return stacked array with only one dimension left
	INPUTS:
	   x  : xarray.DataArray or Dataset to be stacked
	   dim: sole dimension to remain after stacking
	OUTPUTS:
	   stacked: same as x, but stacked
	'''
	dims = []
	for d in x.dims:
		if d != dim:
			dims.append(d)
	return x.stack(stacked=dims)

#t2m_stack = StackArray(x=data_t2m, dim="time")

def ComputeCorr(i, x, y, method="Spearmanr"):
    
    if x.shape == y.shape:
        
        if method == "Spearmanr":
            sloc,ploc = stats.spearmanr(x.isel(stacked=i), y.isel(stacked=i))
            
        elif method == "pearsonr":
            sloc,ploc = stats.pearsonr(x.isel(stacked=i), y.isel(stacked=i))
            
    else:
        
        if method == "Spearmanr":
            sloc,ploc = stats.spearmanr(x.isel(stacked=i), y.isel(stacked=0))
            
        elif method == "pearsonr":
            sloc,ploc = stats.pearsonr(x.isel(stacked=i), y.isel(stacked=0))
        
    return sloc, ploc
    
def StatCorr(x,y,dim=None, return_sig=True):
    
    if len(y.dims) ==1:
        
        sy = y.expand_dims(stacked=[0])
    else:
        if dim is not None:
            sy = StackArray(x=y, dim=dim)
    
    if dim is None or len(x.dims) == 1:
        sx = x.expand_dims(stacked=[0])
    else:
        sx = StackArray(x,dim)
        
    nspace = len(sx.stacked)
    sval, pval = np.zeros(sx.stacked.shape), np.zeros(sx.stacked.shape)
    for i in range(nspace):
        sval[i], pval[i] = ComputeCorr(i, sx, sy)
        
    if nspace > 1:
        pvalx = DataArray(pval, coords=[sx.stacked],name='pval').unstack('stacked')
        svalx = DataArray(sval,coords=[sx.stacked],name='sval').unstack('stacked')
        
    else:
        svalx, pvalx = pval[0], sval[0]
        
    sig_loc  = xr.where(pvalx < 0.05, pvalx, pvalx*np.nan)
    
    sig_loc = sig_loc.sortby("lon")
    
    svalx, pvalx = svalx.sortby("lon"), pvalx.sortby("lon")
    
    
    if return_sig:
        return svalx, pvalx, sig_loc
        
    else:
        return svalx, pvalx
        
        
sval,pval,sig = StatCorr(x=data_t2m, y=nao_index, dim="time")


def StatTest(x,y,test,dim=None,parallel=False):
	'''Compute statistical test for significance between
	   two xr.DataArrays. Testing will be done along dimension with name `dim`
	   and the output p-value will have all dimensions except `dim`.
	   INPUTS:
	      x	 : xr.DataArray for testing.
	      y	 : xr.DataArray or scalar for testing against. Or None for single-ensemble sign test.
	      dim: dimension name along which to perform the test.
	      test:which test to use:
		    'KS' -> Kolmogorov-Smirnov
		    'MW' -> Mann-Whitney
		    'WC' -> Wilcoxon
		    'T'  -> T-test 1 sample with y=mean
                    'sign'->test against sign only.
		  parallel: Run the test in parallel? Requires the parmap package.
	   OUTPUTS:
	      pvalx: xr.DataArray containing the p-values.
		     Same dimension as x,y except `dim`.
	'''
	from xarray import DataArray
	if dim is None or len(x.dims) == 1:
		sx = x.expand_dims(stacked=[0])
		parallel = False
	else:
		sx = StackArray(x,dim)
	if parallel:
		import parmap
	nspace = len(sx.stacked)
	if isinstance(y,DataArray):
		if dim is None or len(y.dims) == 1:
			sy = y.expand_dims(stacked=[0])
		else:
			sy = StackArray(y,dim)
	else:
		sy = None
	if parallel:
		pval,sval = parmap.map(ComputeStat,list(range(nspace)),sx,y,sy,test)
	else:
		pval, sval = np.zeros(sx.stacked.shape), np.zeros(sx.stacked.shape)
		for i in range(nspace):
			pval[i], sval[i] = ComputeStat(i,sx,y,sy,test)
	if nspace > 1:
		pvalx, svalx = DataArray(pval,coords=[sx.stacked],name='pval').unstack('stacked'), DataArray(sval,coords=[sx.stacked],name='sval').unstack('stacked')
        
	else:
		pvalx, svalx = pval[0], sval[0]
	return pvalx, svalx

