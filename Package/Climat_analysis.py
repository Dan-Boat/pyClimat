#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 18:49:19 2021

@author: dboateng
Contains analysis routine required for calculation, extracting variables and domain, masking out areas and certian statistics. 

"""

# importing packages 
import xarray as xr
import os
import pandas as pd
import numpy as np
from scipy import stats
from eofs.xarray import Eof

#importing other routines

try:
    from .Climat_data import *
except:
    from Climat_data import *
    



def extract_var(Dataset, varname, units=None, Dataset_wiso=None, other_data=None, lev_units=None, lev=None):
    """
    

    Parameters
    ----------
    Dataset : TYPE: Dataset
        DESCRIPTION. Processed output containing all variables or certain ones
    varname : TYPE: STR
        DESCRIPTION: Name of variable (eg. temp2:Tempeature, prec:Total precipitation, d18op: d18o in precipitation)
    units : TYPE, optional 
        DESCRIPTION. The default is None. Units of the variable for eg. °C for temperature or mm/month for precipitation
    Dataset_wiso : TYPE, optional
        DESCRIPTION. The default is None. Wiso Dataset is required for calcullating d18op and d18ov   
    other_data: TYPE: datarray, optional 
        DESCRIPTION: Additional data required for calculating units eg. Pa/s--> m/s require temperature data on pressure levels
    lev_units: TYPE: str, optional 
        DESCRIPTION: Units of vertical levels eg. hPa
    lev: TYPE: int, optional 
        DESCRIPTION: The vertical height required for the variable, eg. 500 hPa for winds

    Returns
    -------
    var_data : TYPE: datarray
        DESCRIPTION: returns the variable for etracted or calculated from the Dataset (notmally monthly long-term climatologies)

    """
    
    #Temperature
    if varname == "temp2":
        var_data = Dataset[varname]
        if units is not None:
            if units == "°C":
                var_data = var_data - 273.15 # convert temperature to dec C
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kelvin")
        
    #Total precipitation amount
    elif varname == "prec":
        var_data = Dataset["aprl"] + Dataset["aprc"]
        if units is not None:
            if units == "mm/month":
                var_data = var_data *60*60*24*30 #mm/month
                var_data.attrs["units"] = units
            elif units == "mm/a":
                var_data = var_data *60*60*24*365 #mm/a
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kg/m²s")
    
    # d18o in Precipitation
    elif varname == "d18op":
        var_data_echam = Dataset["aprl"] + Dataset["aprc"]
        if Dataset_wiso is not None:
             var_data_wiso = Dataset_wiso["wisoaprl"][:,1,:,:] + Dataset_wiso["wisoaprc"][:,1,:,:]
        SMOW = 0.2228
        wiso_wtgs = var_data_wiso / var_data_echam
        var_data = ((wiso_wtgs/SMOW) -1)
        
        if units is not None:
            if units == "per mil":
                var_data = var_data *1000 # convert to permil
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is isotopic ratio")
                
    # d18o in water vapour
    elif varname == "d18ov":
        var_data_echam = Dataset["q"][:,29,:,:]
        if Dataset_wiso is not None:
             var_data_wiso = Dataset_wiso["q18o"][:,29,:,:]
        SMOW = 0.2228
        wiso_wtgs = var_data_wiso / var_data_echam
        var_data = ((wiso_wtgs/SMOW) -1)
        
        if units is not None:
            if units == "per mil":
                var_data = var_data *1000 # convert to permil
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is isotopic ratio")
   
    # Evaporation            
    elif varname == "evap":
        var_data = Dataset["evap"]
        if units is not None:
            if units == "mm/month":
                var_data = var_data *60*60*24*30 *-1  #mm/month (positive values)
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is kg/m²s")
   
    # relative humidity near surface            
    elif varname == "relhum":
        var_data = (Dataset["relhum"][:,29,:,:] + Dataset["relhum"][:,30,:,:]) / 2 #mean of the last 2 levels close to surface
        if units is not None:
            if units == "%":
                var_data = var_data *100 #%
                var_data.attrs["units"] = units
            else:
                print("Define unit well or the default is ratio")
                
    # elevation from GEOSP
    elif varname == "elev":
        var_data = Dataset["geosp"] / 9.8  # convert to m**2/s**2 --> m
        if units is not None:
            if units == "km":
                var_data = var_data /1000  #km
                var_data.attrs["units"] = units
            else:
                var_data.attrs["units"] = units
                print("The default units of elevation is used [m]")
                
    # sea land mask
    elif varname == "slm":
        var_data = Dataset["slm"]
        
    # Wind velocity zonal
    elif varname == "u10":
        var_data = Dataset["u10"]
        
    # wind velocity meridoinal
    elif varname == "v10":
        var_data = Dataset["v10"]
        
    elif varname == "wind10":
        var_data = Dataset["wind10"]
    
        #geopotential height at pressure levels
    elif varname == "geopoth":
        var_data = Dataset["geopoth"]
        if lev_units is not None:
            if lev_units == "hPa":
                var_data["lev"] = var_data.lev / 100  #Pa --> hPa
                var_data["lev"].attrs["units"] = "hPa"
            else:
                 raise ValueError("You have defined incorrect units or its not implemented")
    
    # mean sea level pressure    
    elif varname == "slp":
        var_data = Dataset["slp"]
        if lev_units is not None:
            if lev_units == "hPa":
                var_data["lev"] = var_data.lev / 100  #Pa --> hPa
                var_data["lev"].attrs["units"] = "hPa"
            else:
                 raise ValueError("You have defined incorrect units or its not implemented")
        
        if units is not None:
            if units == "hPa":
                var_data /= 100
                var_data.attrs["units"] = units
            else:
                 raise ValueError("You have defined incorrect units or its not implemented")
    # Vertical verlocity    
    elif varname =="omega":  # Pa/s
        var_data = Dataset[varname]
        
        if lev_units is not None:
            if lev_units == "hPa":
                var_data["lev"] = var_data.lev / 100  #Pa --> hPa
                var_data["lev"].attrs["units"] = "hPa"
            else:
                 raise ValueError("You have defined incorrect units or its not implemented")
        
        if units is not None:
            #Other data must be temperature on pressure levels with the same shape as omega
            if units == "m/s":
                rgas=287.058  # m²/s²K
                g=9.80665
                pa = Dataset["lev"]
                other_data["st"]
                rho = pa/(rgas*other_data)
                var_data /= (rho*g)  #Pa/s -->m/s
            else:
                raise ValueError("You have defined incorrect units or its not implemented")
        if lev is not None:
            var_data = var_data.sel(lev=lev)
        else:
            raise ValueError("You have defined incorrect vertical level")
    # meridional wind at vertical levels    
    elif varname == "v":
        var_data = Dataset[varname]
        if lev_units is not None:
            if lev_units == "hPa":
                var_data["lev"] = var_data.lev / 100  #Pa --> hPa
                var_data["lev"].attrs["units"] = "hPa"
            else:
                 raise ValueError("You have defined incorrect units or its not implemented")
        
        if lev is not None:
            var_data = var_data.sel(lev=lev)
        else:
            raise ValueError("You have defined incorrect vertical level")
    # zonal wind at vertical levels    
    elif varname == "u":
        var_data = Dataset[varname]
        if lev_units is not None:
            if lev_units == "hPa":
                var_data["lev"] = var_data.lev / 100  #Pa --> hPa
                var_data["lev"].attrs["units"] = "hPa"
            else:
                 raise ValueError("You have defined incorrect units or its not implemented")
        if lev is not None:
            var_data = var_data.sel(lev=lev)
        else:
            raise ValueError("You have defined incorrect vertical level")
    elif varname == "e/p":
        var_prec = Dataset["aprl"] + Dataset["aprc"]
        var_evap = Dataset["evap"] 
        var_data = var_evap/var_prec
        
    else:
        print("Check the Varname or it has not been implemeted")
            
    return var_data



def compute_lterm_mean(data, time="annual", month=None, season=None, season_calendar = None):
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

    Returns
    -------
    data_ltmean : TYPE: datarray
        DESCRIPTION. Long-term means 

    """
    
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
            data_ltmean = data_ltmean.sel(month=month)
    else:
        print("Define the time period for long-term mean or the annual mean is computed")
        
    return data_ltmean
    




def compute_lterm_diff(data_control, data_main, time="annual", month=None, season=None, season_calendar = None):
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
    
    data_ltmean_control = compute_lterm_mean(data=data_control, time=time, month=month, season=season, season_calendar = season_calendar)
    data_ltmean_main = compute_lterm_mean(data=data_main, time=time, month=month, season=season, season_calendar = season_calendar)
    
    data_ltmean_diff = data_ltmean_main - data_ltmean_control
    
    return data_ltmean_diff
    

def extract_transect(data, maxlon, minlon, maxlat, minlat, sea_land_mask=None, minelev=None, maxelev=None, Dataset=None):
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
    
    if Dataset is not None:
        
        # add mask arrary
        slm = Dataset["slm"]
        geosp = Dataset["geosp"] / 9.8
        data.coords["slm"] = (("lat","lon"), slm[0].data)
        data.coords["geosp"] = (("lat","lon"), geosp[0].data) #same for all time point
    
    #convert lon to -180 to 180
    data = data.assign_coords({"lon": (((data.lon + 180) % 360) - 180)})
    
    lat_range = (data.lat >= minlat) & (data.lat <= maxlat)
    lon_range = (data.lon >= minlon) & (data.lon <= maxlon)
    data_extract = data.where((lat_range & lon_range), drop=True)
    
    # Caution! Since the mask can only be removed when both the lat and lon cordinates are satisfied, some points might 
    # still be there UNLESS the mean is calculated with dim("lat", "lon")!!!
    
    # Alternative solution implemented!
    # I have replace all the false conditions with np.nan values and must be removed when converted to DataFrame
    
    if Dataset is not None:
        
        print("Using dataset for masking")
        
        if sea_land_mask is not None:
            if sea_land_mask in (["yes", "YES", "Yes", "yep"]):
                data_extract = xr.where(data_extract["slm"] == 1, data_extract, data_extract*np.nan)
            elif sea_land_mask in (["no", "NO", "No", "nope"]):
                data_extract = xr.where(data_extract.slm == 0, data_extract, data_extract*np.nan)
            else:
                print("define yes or no for sea land mask")
        if minelev is not None:
            data_extract = xr.where(data_extract.geosp >= minelev, data_extract, data_extract*np.nan)
            
        if maxelev is not None:
            data_extract = xr.where(data_extract.geosp <= maxelev, data_extract, data_extract*np.nan)
        
        
    return data_extract
    

def extract_profile(data, maxlon, minlon, maxlat, minlat, dim, to_pandas=None, sea_land_mask=None, minelev=None, maxelev=None, 
                    Dataset=None):
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
        data_prof = data_extract.mean(dim="lon")
    elif dim in ["lon", "longitude"]:
        print("Computing the mean across latitude")
        data_prof = data_extract.mean(dim="lat")
    else:
        print("Define the dimension to extract the profile")
        
    if to_pandas is not None:
        if to_pandas in ["yes", "Yes", "YES", "yep"]:
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
            print("Use yes if you want to convert to Dataframe for ploting or saving in csv")
    
            
    return data_prof


def linregression(data_x, data_y, season=None, month=None, return_yhat=None):
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
    print("y = {:.5f}x [‰/100m]+ {:.2f}, r²={:.2f}".format(regression_stats.slope*100, regression_stats.intercept, 
                                                           regression_stats.rvalue*-1))
    if return_yhat is not None:
        if return_yhat in ["Yes", "yes", "yep"]:
            yhat = regression_stats.slope * data_x  + regression_stats.intercept
            
            df_x_y_yhat = pd.DataFrame(columns=["X", "Y", "yhat"])
            
            df_x_y_yhat["yhat"] = yhat
            df_x_y_yhat["X"] = data_x
            df_x_y_yhat["Y"] = data_y
            
            return regression_stats, df_x_y_yhat
        else:
            print("Only regression stats are return")
            
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
                                minlon=None, maxlat=None, minlat=None):
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
    
    stats_result = xr.Dataset(coords={"lon": (["lon"], dataA.lon), "lat": (["lat"], dataA.lat)})
    
    stats_result["t_statistic"] = (["lat", "lon"], stats_data)
    
    if return_pvalue == True:
        stats_result["pvalue"] = (["lat", "lon"], pvalue_data)
        
        return stats_result
    else:
        return stats_result



def PCA_analysis(*args):
    pass

def remove_lapse_rate_effect(*arg):
    pass
