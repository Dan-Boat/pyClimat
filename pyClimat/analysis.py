#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 18:49:19 2021

@author: dboateng
Contains analysis routine required for calculation, extracting variables and domain, masking out areas 
and certian statistics. 

"""

# importing packages 
import xarray as xr
import os
import pandas as pd
import numpy as np
from scipy import stats
from eofs.xarray import Eof
import statsmodels.api as sm

#importing other routines

try:
    from .data import *
    from .utils import vert_coord_convertion
    from .variables import extract_var
except:
    from data import *
    from utils import vert_coord_convertion
    from variables import extract_var


def compute_spatial_means(dataset, varname, time="annual", land=False, ocean=False, extract_region=False, lon_min=None, lon_max=None,
                          lat_min=None, lat_max=None, units=None, Dataset_wiso=None):
    
    # select variable 
    data = extract_var(Dataset=dataset, varname=varname, units=units, Dataset_wiso=Dataset_wiso,
                       )
    
    if time == "season":
        month_length = data.time.dt.days_in_month

        # Calculate the weights by grouping by 'time.season'
        time_weights = month_length.groupby('time.season') / month_length.groupby('time.season').sum()

        # Test that the sum of the weights for each season is 1.0
        np.testing.assert_allclose(time_weights.groupby('time.season').sum().values, np.ones(4))

        # Calculate the weighted average
        data_ltmean = (data * time_weights).groupby('time.season').sum(dim='time')
    
    elif time == "month":
        data_ltmean = data.groupby("time.month").mean("time")
        
        if month is not None:
            print("Calculating the monthly long-term mean for the month number:", month)
            
            if month == "JJAS": # useful for WAM analysis
                data_ltmean = data_ltmean.sel(month=data_ltmean.month.isin([6, 7, 8, 9]))
                data_ltmean = data_ltmean.mean(dim="month")
                
            elif isinstance(month, int):    
                data_ltmean = data_ltmean.sel(month=month)
                
            else:
                raise ValueError("The define month parameter is not right")
                
    else:
        data_ltmean = data.mean(dim="time")
                
    if hasattr(data_ltmean, "longitude"):
        data_ltmean = data_ltmean.rename({"longitude":"lon", "latitude":"lat"})
        
        
    if extract_region==True:   
        data_ltmean = data_ltmean.assign_coords({"lon": (((data.lon + 180) % 360) - 180)})
    
        data_ltmean = data_ltmean.where((data.lat >= minlat) & (data.lat <= maxlat), drop=True)
        data_ltmean = data_ltmean.where((data.lon >= minlon) & (data.lon <= maxlon), drop=True)
        
    
    if land==True:
        land_mask_ds = dataset["slm"][0]
        data_ltmean = data_ltmean.where(land_mask_ds == 1)
        
    if ocean==True:
        land_mask_ds = dataset["slm"][0]
        data_ltmean = data_ltmean.where(land_mask_ds == 0)
        
    
    # compute area-weighted means
    weights = np.cos(np.deg2rad(data_ltmean.lat)) #np.sqrt(np.abs(np.cos(data_raw.lat*np.pi/180)))

    weights.name = "weights"

    means = data_ltmean.weighted(weights).mean(dim=("lon", "lat"), skipna=True)
    
    
    return means
        
        

    
    
    
def compute_lterm_mean(data, time="annual", month=None, season=None, season_calendar = None, time_range=None):
    """
    

    Parameters
    ----------
    data : TYPE: datarray
        DESCRIPTION. The var_data extracted from dataset 
    time : TYPE: str, optional
        DESCRIPTION. The default is "annual". or season, month can be used for long-term estimates
    month : TYPE, optional
        DESCRIPTION. The default is None.
    season : TYPE, optional
        DESCRIPTION. The default is None.
    season_calendar : TYPE, optional
        DESCRIPTION. The default is None. Use standard if you want to consider the days of the month into consideration
        
    time_range: datetime
                The time range to extract for computing the mean

    Returns
    -------
    data_ltmean : TYPE: datarray
        DESCRIPTION. Long-term means 

    """
    
    if time_range is not None:
        
        print("------selecting time range for the long-term means ---------")
        
        data = data.sel(time =time_range)
        
    
    data_ltmean = data.mean(dim="time")
    
    if time =="season":
        if season_calendar is not None:
            if season_calendar == "standard":
                # Make a DataArray with the number of days in each month, size = len(time)
                month_length = data.time.dt.days_in_month

                # Calculate the weights by grouping by 'time.season'
                weights = month_length.groupby('time.season') / month_length.groupby('time.season').sum()

                # Test that the sum of the weights for each season is 1.0
                np.testing.assert_allclose(weights.groupby('time.season').sum().values, np.ones(4))

                # Calculate the weighted average
                data_ltmean = (data * weights).groupby('time.season').sum(dim='time')
            else:
                print("Define the season calendar well eg. standard")
        else:
            print("Calculating the seasonal mean without considering the number of days in a month")

            data_ltmean = data.groupby("time.season").mean("time")
        
        if season is not None:
            print("Calculating the seasonal long-term mean for:", season)
            data_ltmean = data_ltmean.sel(season = season)
        
            
    elif time == "month":
        data_ltmean = data.groupby("time.month").mean("time")
        
        if month is not None:
            print("Calculating the monthly long-term mean for the month number:", month)
            
            if month == "JJAS": # useful for WAM analysis
                data_ltmean = data_ltmean.sel(month=data_ltmean.month.isin([6, 7, 8, 9]))
                data_ltmean = data_ltmean.mean(dim="month")
                
            elif month == "MAM": # useful for EAS
                data_ltmean = data_ltmean.sel(month=data_ltmean.month.isin([3, 4, 5]))
                data_ltmean = data_ltmean.mean(dim="month")
                
            elif isinstance(month, int):    
                data_ltmean = data_ltmean.sel(month=month)
                
            else:
                raise ValueError("The define month parameter is not right")
    else:
        print("Define the time period for long-term mean or the annual mean is computed")
        
    return data_ltmean
    




def compute_lterm_diff(data_control, data_main, time="annual", month=None, season=None, season_calendar=None,
                       time_range_main=None, time_range_control=None):
    """
    

    Parameters
    ----------
    data_control : TYPE: datarray
        DESCRIPTION. Reference or control data
    data_main : TYPE: dataarray
        DESCRIPTION. : Main module run data
    time : TYPE: STR, optional
        DESCRIPTION. The default is "annual". But can be changed to season or month
    month : TYPE:INT, optional
        DESCRIPTION. The default is None. Define the specific month number to be computed or all will be eatimated 
    season : TYPE: STR, optional
        DESCRIPTION. The default is None. Define the specific season (eg. DJF, JJA, MAM, SON) number to be computed or all will be eatimated 
    season_calendar : TYPE:STR, optional
        DESCRIPTION. The default is None. Or Standard to consider the number of days in a month

    Returns
    -------
    data_ltmean_diff : TYPE: datarray
        DESCRIPTION. Long-tem difference

    """
    # Use compute_lterm_diff to calculate the long-term mean before the the difference
    
    data_ltmean_control = compute_lterm_mean(data=data_control, time=time, month=month, 
                                             season=season, season_calendar = season_calendar, 
                                             time_range=time_range_control)
    
    data_ltmean_main = compute_lterm_mean(data=data_main, time=time, month=month, 
                                          season=season, season_calendar = season_calendar,
                                          time_range=time_range_main)
    
    data_ltmean_diff = data_ltmean_main - data_ltmean_control
    
    return data_ltmean_diff
    

def extract_transect(data, maxlon, minlon, maxlat, minlat, sea_land_mask=False, minelev=None, maxelev=None, Dataset=None):
    """
    This function extract grid points base on coordinate extents or land sea masks or max, min elevations: it can be used to estimate 
    the statistics of a selected domain like the Alps or Andes!

    Parameters
    ----------
    data : TYPE: dataarray
        DESCRIPTION: Data to extract transect from base on coordinates, elevation or land sea masks
    maxlon : TYPE: float
        DESCRIPTION.: Maximum longitude
    minlon : TYPE: float
        DESCRIPTION: Minimum longitude
    maxlat : TYPE: float
        DESCRIPTION:Maximum latitude
    minlat : TYPE: float
        DESCRIPTION: Minimum latitude
    sea_land_mask : TYPE: str, optional
        DESCRIPTION. The default is None. Yes, means that the land mask will be selected and No means the sea nask points 
        will be selected
    minelev : TYPE, optional
        DESCRIPTION. The default is None. To select data points base on the minimum elevation value
    maxelev : TYPE. float, optional
        DESCRIPTION. The default is None. To select data points base on the maximum elevation value
    Dataset : TYPE: float, optional
        DESCRIPTION. The default is None. Dataset containing geosp and slm for masking out elevation condition and continental values

    Returns
    -------
    data_extract : TYPE
        DESCRIPTION.

    """
    
    #convert lon to -180 to 180
    if hasattr(data, "longitude"):
        data = data.rename({"longitude": "lon", "latitude":"lat"})
    
    if Dataset is not None:
        
        if hasattr(Dataset, "longitude"):
            Dataset = Dataset.rename({"longitude": "lon", "latitude":"lat"}) 
        # add mask arrary (sea mask and elevation is the same for the time, therefore, extrating one index is fine)
        
        if sea_land_mask == True:
            slm = Dataset["slm"]
            data.coords["slm"] = (("lat","lon"), slm[0].data)
        if minelev or maxelev is not None:
            
            geosp = Dataset["geosp"] / 9.8
            data.coords["geosp"] = (("lat","lon"), geosp[0].data) #same for all time point
    
    
        
    data = data.assign_coords({"lon": (((data.lon + 180) % 360) - 180)})
    
    lat_range = (data.lat >= minlat) & (data.lat <= maxlat)
    lon_range = (data.lon >= minlon) & (data.lon <= maxlon)
    
    data_extract = data.where((lat_range & lon_range), drop=True)
    
    # Caution! Since the mask can only be removed when both the lat and lon cordinates are satisfied, some points might 
    # still be there UNLESS the mean is calculated with dim("lat", "lon")!!!
    
    # Alternative solution implemented!
    # I have replace all the false conditions with np.nan values and must be removed when converted to DataFrame
    
    if sea_land_mask == True:
        if hasattr(data_extract, "slm"):
            data_extract = xr.where(data_extract["slm"]==1, data_extract, data_extract*np.nan)
            
            data_extract = data_extract.drop_vars(["slm"])
    if minelev and maxelev is not None:
        if hasattr(data_extract, "geosp"):
            data_extract = xr.where((data_extract.geosp >= minelev) & (data_extract.geosp <= maxelev), data_extract, data_extract*np.nan)
            data_extract = data_extract.drop_vars(["geosp"])
    elif minelev is not None:
        if hasattr(data_extract, "geosp"):
            data_extract = xr.where(data_extract.geosp >= minelev, data_extract, data_extract*np.nan)
            data_extract = data_extract.drop_vars(["geosp"])
        
    elif maxelev is not None:
        if hasattr(data_extract, "geosp"):
            data_extract = xr.where(data_extract.geosp <= maxelev, data_extract, data_extract*np.nan)
            data_extract = data_extract.drop_vars(["geosp"])
          
    return data_extract
    

def extract_profile(data, maxlon, minlon, maxlat, minlat, dim, to_pandas=True, sea_land_mask=False, minelev=None, maxelev=None, 
                    Dataset=None, method="mean"):
    """
    

    Parameters
    ----------
    data : TYPE: dataarray
        DESCRIPTION: Data to extract transect from base on coordinates, elevation or land sea masks
    maxlon : TYPE: float
        DESCRIPTION.: Maximum longitude
    minlon : TYPE: float
        DESCRIPTION: Minimum longitude
    maxlat : TYPE: float
        DESCRIPTION:Maximum latitude
    minlat : TYPE: float
        DESCRIPTION: Minimum latitude
    sea_land_mask : TYPE: str, optional
        DESCRIPTION. The default is None. Yes, means that the land mask will be selected and No means the sea nask points 
        will be selected
    minelev : TYPE, optional
        DESCRIPTION. The default is None. To select data points base on the minimum elevation value
    maxelev : TYPE. float, optional
        DESCRIPTION. The default is None. To select data points base on the maximum elevation value
    Dataset : TYPE: float, optional
        DESCRIPTION. The default is None. Dataset containing geosp and slm for masking out elevation condition and continental values

    dim : TYPE: str
        DESCRIPTION. lat ot lon depending on the axis of the profile
    to_pandas : TYPE:str, optional (recommended for plotting)
        DESCRIPTION. The default is None. yes if you want the data to be stored in DataFrame or Pandas Series

    Returns
    -------
    data_prof : TYPE: datarray or DataFrame or Pandas Series
        DESCRIPTION. Data extracted along a profile (lat or lon)

    """
    
    # use extract transect to get the domain or points 
    data_extract = extract_transect(data=data, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, sea_land_mask=sea_land_mask, 
                                    minelev=minelev, maxelev=maxelev, Dataset=Dataset)
    if dim in ["lat", "latitude"]:
        print("Computing the mean across longitude")
        
        if method =="std":
            data_prof = data_extract.std(dim="lon")
            
        else:
            data_prof = data_extract.mean(dim="lon")
            
            
    elif dim in ["lon", "longitude"]:
        print("Computing the mean across latitude")
        
        if method =="std":
            data_prof = data_extract.std(dim="lat")
            
        else:
            data_prof = data_extract.mean(dim="lat")
        
    else:
        print("Define the dimension to extract the profile")
        
    if to_pandas ==True:
        data_prof = data_prof.to_pandas()
        data_prof = data_prof.T
        
        # sort values base on the dim 
            
        if isinstance(data_prof, pd.Series):
            data_prof = data_prof.sort_index(axis = 0, ascending = True)
        elif isinstance(data_prof, pd.DataFrame):
            data_prof = data_prof.sort_values(by=[dim])
        else:
            print("Check the instance of the data to be sorted")
    else:
        print("Profile data not stored in DataFrame")
    
            
    return data_prof

def extract_vertical_section(data, maxlon, minlon, maxlat, minlat, dim, sea_land_mask=False, minelev=None, maxelev=None, 
                    Dataset=None, season = None, month=None, lev=None):
    
    data_extract = extract_transect(data=data, maxlon=maxlon, minlon=minlon, maxlat=maxlat, minlat=minlat, sea_land_mask=sea_land_mask, 
                                    minelev=minelev, maxelev=maxelev, Dataset=Dataset)
    
    if dim in ["lat", "latitude"]:
        print("Computing the mean across longitude")
        data_sect = data_extract.mean(dim="lon")
    elif dim in ["lon", "longitude"]:
        print("Computing the mean across latitude")
        data_sect = data_extract.mean(dim="lat")
    else:
        print("Define the dimension to extract the profile")
        
    # select season if required
    if season is not None:
        data_sect = data_sect.sel(season=season)
        
    if month is not None: 
        data_sect = data_sect.sel(month=month)
    
    if lev is not None:
        data_sect = data_sect.sel(lev=lev)
    # convert to pandas 
    df = data_sect.to_pandas()
    df = df.T  # transpose axis
    
    # sortby using the columns 
    if isinstance(df, pd.Series):
        df = df.sort_values(axis=1, ascending=True)
        df = df.sort_values(axis=0, ascending=True)
    elif isinstance(df, pd.DataFrame):
        df = df.sort_values(by=[dim], axis=0)
        #df = df.sort_values(by=["lev"], axis=1)
    
    return df
        
    
    
    
    
    
def linregression(data_x, data_y, season=None, month=None, return_yhat=True, get_ci=False):
    """
    

    Parameters
    ----------
    data_x : TYPE: datarray
        DESCRIPTION. The x-axis data for fitting 
    data_y : TYPE: datarray
        DESCRIPTION. The y-axis data for fitting 
    season : TYPE, optional (of If specific season the data is required for fitting)
        DESCRIPTION. The default is None. Date must be in seasonal coordinates for time
    month : TYPE, optional (of If specific month the data is required for fitting)
        DESCRIPTION. The default is None. Date must be in monthly coordinates for time
    return_yhat : TYPE, optional or if DataFrame containing all the fitting data and predictions are required
        DESCRIPTION. The default is None.

    Returns
    -------
    TYPE: Scipy.stats output or plus DataFrame
        DESCRIPTION.

    """
    if season is not None:
        data_x = data_x.sel(season = season).to_pandas().values.ravel()
        data_y = data_y.sel(season = season).to_pandas().values.ravel()
        
    elif month is not None: 
        data_x = data_x.sel(month = month).to_pandas().values.ravel()
        data_y = data_y.sel(month = month).to_pandas().values.ravel()
    else:
        data_x = data_x.to_pandas().values.ravel()
        data_y = data_y.to_pandas().values.ravel()
    
    # remove nan values before fitting
    
    data_x = data_x[~np.isnan(data_x)]
    data_y = data_y[~np.isnan(data_y)]
    
    regression_stats = stats.linregress(data_x, data_y)
    
    # get the prediction interval
    data_xx = sm.add_constant(data_x)
    model = sm.OLS(data_y, data_xx)
    r = model.fit()
    
    predictions = r.get_prediction(exog=data_xx, transform=False).summary_frame(alpha=0.05)
    predictions["X"] = data_x
    
    print("y = {:.5f}x [‰/100m]+ {:.2f}, r²={:.2f}".format(regression_stats.slope*100, regression_stats.intercept, 
                                                           regression_stats.rvalue*-1))
    if return_yhat == True:
        
        yhat = regression_stats.slope * data_x  + regression_stats.intercept
        
        df_x_y_yhat = pd.DataFrame(columns=["X", "Y", "yhat"])
        
        df_x_y_yhat["yhat"] = yhat
        df_x_y_yhat["X"] = data_x
        df_x_y_yhat["Y"] = data_y
        
        if get_ci == True:
            return regression_stats, df_x_y_yhat, predictions 
            
        return regression_stats, df_x_y_yhat
            
    else:
        print("Only regression stats are return")
        return regression_stats


def EOF_analysis(data, maxlon, minlon, maxlat, minlat, return_variance=False, return_pcs=False,
                 season=None, standardized=None, apply_coslat_weights=None,
                 neofs=None, pcscaling=None, neigs=None, npcs=None, lev=None):
    """
    # the projectField function can be used to generate corresponding set of pseudo-PCs using different data field

    Parameters
    ----------
    data : TYPE: datarray
        DESCRIPTION. Dataset required for EOF analysis (eg. slp or geopoth500)
    maxlon : TYPE: float
        DESCRIPTION.: Maximum longitude
    minlon : TYPE: float
        DESCRIPTION: Minimum longitude
    maxlat : TYPE: float
        DESCRIPTION:Maximum latitude
    minlat : TYPE: float
        DESCRIPTION: Minimum latitude
    return_variance : TYPE: Boolean, optional
        DESCRIPTION. The default is False. If estimated varainces of the eofs are required as ouput
    return_pcs : TYPE: Boolean, optional
        DESCRIPTION. The default is False. If the extracted pca series are required as ouput
    season : TYPE:STR, optional
        DESCRIPTION. The default is None. Name of the season eg. DFJ, JJA
    standardized : TYPE: Boolean, optional
        DESCRIPTION. The default is None. True to standardized anomalies before EOF analysis 
    apply_coslat_weights : TYPE: Boolean, optional
        DESCRIPTION. The default is None. True to apply coslat area weights before the EOF analysis
    neofs : TYPE: Float, optional
        DESCRIPTION. The default is None. The no. of PCA to perform on the dataset 
    pcscaling : TYPE: Int, optional
        DESCRIPTION. The default is None. 0 : Unsclaed PCS, 1: Scaled to Unit variance, 2: PCs are 
        multiplied by the square-root of their eigen values 
    neigs : TYPE: Int, optional
        DESCRIPTION. The default is None. the no. of eigenvalues to return fraction variance
    npcs : TYPE:Int, optional
        DESCRIPTION. The default is None. The no. of pcs retrieve 
    lev : TYPE: float, optional
        DESCRIPTION. The default is None. Vertical level if the dataset is on hybrid levels (eg. 500 for geopoth)

    Returns
    -------
    TYPE
        DESCRIPTION. eofs: covariance matrix between the npcs time series and eofs input time series
        pcs: Principal Component time series
        var_frac: variance fraction of the estimated eigen values 

    """
    if lev is not None:
        data = data.sel(lev=lev)
    
    if season is not None:
        data = data.groupby("time.season")
        data = data[season]
    # extracting transects for the analysis
    # converting lon to -180 t0 180 
    
    data = data.assign_coords({"lon": (((data.lon + 180) % 360) - 180)})
    
    data = data.where((data.lat >= minlat) & (data.lat <= maxlat), drop=True)
    data = data.where((data.lon >= minlon) & (data.lon <= maxlon), drop=True)
    
    # calculation of anomalies
    
    data_mean = data.mean(dim="time")
    
    data_anomalies = data - data_mean
    
    # standardized data
    if standardized == True:
        data_anomalies_std = data_anomalies.std(dim="time")
        
        data_anomalies /= data_anomalies_std
        
    # apply weights if equal area using cosalt
    
    if apply_coslat_weights == True:
        wtgs = np.sqrt(np.abs(np.cos(data_anomalies.lat * np.pi / 180)))
        
        data_anomalies = data_anomalies * wtgs
        
    # applying the EOF function (ref: http://doi.org/10.5334/jors.122)
    
    if True in data_anomalies.isnull():
        data_anomalies = data_anomalies.dropna(dim="time")
        
    Solver = Eof(data_anomalies)
    
    if all(parameter is not None for parameter in [neofs, pcscaling]):
        eofs_cov = Solver.eofsAsCovariance(neofs=neofs, pcscaling=pcscaling)
    else:
        eofs_cov = Solver.eofsAsCovariance()
    #sort with the right lon positoin --> becuase of changing the lon to -+180
    
    eofs_cov = eofs_cov.sortby(eofs_cov.lon)
    
    if return_variance == True:
        if neigs is not None:
            var_frac = Solver.varianceFraction(neigs = neigs)
        else:
            var_frac = Solver.varianceFraction()
            
        var_frac = var_frac.to_pandas()
            
            
    if return_pcs == True:
        if all(parameter is not None for parameter in [npcs, pcscaling]):
            pcs = Solver.pcs(pcscaling = pcscaling, npcs = npcs)
        else:
            pcs = Solver.pcs()
        pcs = pcs.to_pandas()
    if return_pcs == True and return_variance == True:
        return eofs_cov, pcs, var_frac
    elif return_pcs == True and return_variance == False:
        return eofs_cov, pcs
    elif return_pcs == False and return_variance == True:
        return eofs_cov, var_frac
    else:
        return eofs_cov



def correlation(dataA, dataB, max_pvalue=0.1, use_spearmanr=False, use_pearsonr=False, 
                return_pvalue=False, maxlon=None, minlon=None, maxlat=None, minlat=None):
    """
    

    Parameters
    ----------
    dataA : TYPE: Dataarray (3D)
        DESCRIPTION. Comparison data 1
    dataB : TYPE: Dataaray (3D)
        DESCRIPTION. Comparison data 2
    max_pvalue : TYPE: Float, optional
        DESCRIPTION. The default is 0.1. The confidence interval for correlation estimation eg. 0.05 for 95%
    use_spearmanr : TYPE: Boolean, optional
        DESCRIPTION. The default is False. True to use spearman correlation
    use_pearsonr : TYPE: Boolean, optional
        DESCRIPTION. The default is False. True to use pearson correlation 
    return_pvalue : TYPE: Boolean, optional
        DESCRIPTION. The default is False. True to retrieve pvalue as an ouput variable
    maxlon : TYPE: float
        DESCRIPTION.: Maximum longitude
    minlon : TYPE: float
        DESCRIPTION: Minimum longitude
    maxlat : TYPE: float
        DESCRIPTION:Maximum latitude
    minlat : TYPE: float
        DESCRIPTION: Minimum latitude

    Raises
    ------
    ValueError
        DESCRIPTION. If the required stats module for correlation analysis is not defined

    Returns
    -------
    stats_result : TYPE: datarray
        DESCRIPTION. Contians correlation map distribution and corresponding pvalues

    """
    #extract specific domain for correlation
    
    if all(par is not None for par in [maxlon, minlon, maxlat, minlat]):
        dataA = dataA.where((dataA.lat >= minlat) & (dataA.lat <= maxlat), drop=True)
        dataA= dataA.where((dataA.lon >= minlon) & (dataA.lon <= maxlon), drop=True)
        
        dataB = dataB.where((dataB.lat >= minlat) & (dataB.lat <= maxlat), drop=True)
        dataB = dataB.where((dataB.lon >= minlon) & (dataB.lon <= maxlon), drop=True)
        
    # define empty matrix to store r
    
    corr_data = np.zeros((dataA.shape[1], dataA.shape[2]), dtype=float)
    
    pvalue_data = np.zeros((dataA.shape[1], dataA.shape[2]), dtype=float)
    
    #looping over logitude and latitude 
    for i in range(0, dataA.shape[1]):  # lat
        for j in range(0, dataB.shape[2]): #lon
            if use_spearmanr==True:
                r,p = stats.spearmanr(dataA.isel(lat=[i], lon=[j]).values.ravel(), dataB.isel(lat=[i], lon=[j]).values.ravel())
                
            elif use_pearsonr == True:
                r,p = stats.pearsonr(dataA.isel(lat=[i], lon=[j]).values.ravel(), dataB.isel(lat=[i], lon=[j]).values.ravel())
            else:
                raise ValueError("Stats module not available, use_spearman or pearsonr")
            
            if p >= max_pvalue:
                corr_data[i,j] = np.nan
            else:
                corr_data[i,j] = r
            
            pvalue_data[i,j] = p
            
    #creating dataset to store values 
    
    stats_result = xr.Dataset(coords={"lon": (["lon"], dataA.lon), "lat": (["lat"], dataA.lat)})
    
    stats_result["correlation"] = (["lat", "lon"], corr_data)
    
    if return_pvalue == True:
        stats_result["pvalue"] = (["lat", "lon"], pvalue_data)
        
        return stats_result
    else:
        return stats_result
                
        

def student_t_test_btn_datasets(dataA, dataB, max_pvalue=0.1, return_pvalue=False, maxlon=None, 
                                minlon=None, maxlat=None, minlat=None, time=None):
    """
    

    Parameters
    ----------
   dataA : TYPE: Dataarray (3D)
        DESCRIPTION. Comparison data 1
    dataB : TYPE: Dataaray (3D)
        DESCRIPTION. Comparison data 2
    max_pvalue : TYPE: Float, optional
        DESCRIPTION. The default is 0.1. The confidence interval for correlation estimation eg. 0.05 for 95%
    return_pvalue : TYPE: Boolean, optional
        DESCRIPTION. The default is False. True to retrieve pvalue as an ouput variable
    maxlon : TYPE: float
        DESCRIPTION.: Maximum longitude
    minlon : TYPE: float
        DESCRIPTION: Minimum longitude
    maxlat : TYPE: float
        DESCRIPTION:Maximum latitude
    minlat : TYPE: float
        DESCRIPTION: Minimum latitude

    Returns
    -------
    stats_result : TYPE
        DESCRIPTION.

    """
    
    if time == "JJAS":
        dataA = dataA.sel(time= dataA.time.dt.month.isin([6, 7, 8, 9]))
        dataB = dataB.sel(time= dataB.time.dt.month.isin([6, 7, 8, 9]))
    
    if all(par is not None for par in [maxlon, minlon, maxlat, minlat]):
        
        dataA = dataA.assign_coords({"lon": (((dataA.lon + 180) % 360) - 180)})
        dataB = dataB.assign_coords({"lon": (((dataB.lon + 180) % 360) - 180)})
        
        
        dataA = dataA.where((dataA.lat >= minlat) & (dataA.lat <= maxlat), drop=True)
        dataA= dataA.where((dataA.lon >= minlon) & (dataA.lon <= maxlon), drop=True)
        
        dataB = dataB.where((dataB.lat >= minlat) & (dataB.lat <= maxlat), drop=True)
        dataB = dataB.where((dataB.lon >= minlon) & (dataB.lon <= maxlon), drop=True)
        
    # define empty matrix to store r
    
    stats_data = np.zeros((dataA.shape[1], dataA.shape[2]), dtype=float)
    
    pvalue_data = np.zeros((dataA.shape[1], dataA.shape[2]), dtype=float)
    
    #looping over logitude and latitude 
    for i in range(0, dataA.shape[1]):  # lat
        for j in range(0, dataB.shape[2]): #lon
            s,p = stats.ttest_ind(dataA.isel(lat=[i], lon=[j]).values.ravel(), dataB.isel(lat=[i], lon=[j]).values.ravel())
            
            if p >= max_pvalue:
                stats_data[i,j] = np.nan
            else:
                stats_data[i,j] = s
            
            pvalue_data[i,j] = p
            
    #creating dataset to store values 
    
    stats_result = xr.Dataset(coords={"lon": (["lon"], dataA.lon.data), "lat": (["lat"], dataA.lat.data)})
    
    stats_result["t_statistic"] = (["lat", "lon"], stats_data)
    
    if return_pvalue == True:
        stats_result["pvalue"] = (["lat", "lon"], pvalue_data)
        
        return stats_result
    else:
        return stats_result

