# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 15:15:58 2023

@author: dboateng
"""
import os
import xarray as xr 
import pandas as pd 
import numpy as np 



# paths
path_to_test_data = "E:/Datasets/Equitech/OCO2_L2/raw_OCO2_L2_Lite_FP_10r_oco2_LtCO2_140906_B10206Ar_200730214713s.nc4"
path_to_test_data_l3 = "E:/Datasets/Equitech/OCO2_L2/raw_OCO2_GEOS_L3CO2_DAY_10r_oco2_GEOS_L3CO2_day_20150101_B10206Ar.nc4"


data = xr.open_dataset(path_to_test_data_l3)

data_l2 = xr.open_dataset(path_to_test_data)