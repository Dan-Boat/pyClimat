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


filename_lterm = "1003_1017_1m_mlterm.nc"

# read long-term means

PI_data = read_from_path(pi_path, filename_lterm)
LGM_data = read_from_path(lgm_path, filename_lterm)
PLIO_data = read_from_path(plio_path, filename_lterm)
MH_data = read_from_path(mh_path, filename_lterm)


# read the pressure level files
filename_plev_lterm = "1003_1017_1m_mlterm_plev.nc"

PI_plev_data = read_from_path(pi_path, filename_plev_lterm)
LGM_plev_data = read_from_path(lgm_path, filename_plev_lterm)
PLIO_plev_data = read_from_path(plio_path, filename_plev_lterm)
MH_plev_data = read_from_path(mh_path, filename_plev_lterm)



