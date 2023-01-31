# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:46:22 2023

@author: dboateng

Perform all the require analysis (eg. EOFs, sorting of indixes, correlation, regional means)
"""

import os 
import numpy as np 
import pandas as pd 
import xarray as xr 


from pyClimat.analysis import extract_var
from pyClimat.stats import EOF_standard

from read_data import PD_data, LGM_data, PI_data

#extract variable from data 
PD_slp = extract_var(PD_data, "slp", units="hPa")
LGM_slp = extract_var(LGM_data, "slp", units="hPa")
PI_slp = extract_var(PI_data, "slp", units="hPa")


#applying EOF class

PD_EOF = EOF_standard(data=PI_slp, weights=True, standardize=True, 
                      extract_region=True, extract_season=True, neofs=4)

PD_EOF.select_time_and_region(maxlon=60, minlon=-80, maxlat=80, minlat=20, time="season", 
                              season="DJF", month="ONDJFM")
PD_EOF.calculate_anomalies()
#PD_EOF.eof_solver(method="xeofs", apply_promax=True)
PD_EOF.eof_solver(method="Eof", apply_promax=True)
eofs = PD_EOF.eofs(eofscaling=2)
pcs = PD_EOF.pcs(pscaling=0)