# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:35:52 2023

@author: dboateng

This module contains all the statistic functions used in pyClimat (eg. EOF, lingress, t-test)
"""

# To do list
#1. Implement EOF class that has the various methods
#2. Implement statistical testing class for xarray datasets


# importing packages 
import xarray as xr
import os
import pandas as pd
import numpy as np
from scipy import stats
from eofs.xarray import Eof


class EOF_standard():
    
    def __init__(self, data, weights=True, standardize=True,
                 extract_region=True, extract_season=True, neofs=4,):
        
        self.data = data
        
        if not isinstance(self.data, xr.DataArray):
            raise TypeError("The X field must be a datarray object or type")
            
        self.weights = weights
        self.standardize = standardize
        self.extract_region = extract_region
        self.extract_season = extract_season
    
    
    
    
    def select_time_and_region(self, maxlon, minlon, maxlat, minlat, time="season", 
                               month=None, season=None):
        
        if self.extract_region:
            
            data = self.data.assign_coords({"lon": (((self.data.lon + 180) % 360) - 180)})
    
            data = data.where((data.lat >= minlat) & (data.lat <= maxlat), drop=True)
            data = data.where((data.lon >= minlon) & (data.lon <= maxlon), drop=True)
            
        if self.extract_season:
            
            if time =="season":
                data = data.groupby("time.season")
                
                if season is not None:
                    data = data.sel(season=season)
                    
                else:
                    data = data.sel(season="DJF")  # set to NH Winter as default
                    
            elif time == "month":
                
                data = data.groupby("time.month")
                
                if month == "ONDJFM":
                    data = data.sel(month=data.month.isin([10,11,12,1,2,3]))
                    
                elif month == "AMJJAS":
                    data = data.sel(month=data.month.isin([14,5,6,7,8,9]))
                    
                elif isinstance(month, int):
                    data = data.sel(month=month)
                    
                elif month == "DJF":
                    data = data.sel(month=data.month.isin([12,1,2]))
                    
                else:
                    raise ValueError("The define month parameter is not recognized")
        
        self.data = data                
        
    
    def apply_standardize(self):
        pass
    
    def apply_weights():
        pass
    
    def explained_variance():
        pass
    
    def explained_variance_ratio():
        pass
    
    def eof_solver():
        pass
    
    def pcs():
        pass
    
    def eofs_as_correlation(self):
        pass
    
    def reconstuct_X():
        pass
    
    def project_X_onto_eofs():
        pass
    

class EOF_sklearn():
    pass



class EOF_rotated():
    pass



    