#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:28:50 2022

@author: dboateng

Reads all the model output for the late cenozoic (LGM, MH, PLIO, PI)
Present-day and future scenarios would be updated for further analysis
"""
#import modules
import os 
import xarray as xr
from pyClimat.data import read_from_path

# define paths 
main_path = "/home/dboateng/Model_output_pst"
lgm_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
plio_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
mh_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
