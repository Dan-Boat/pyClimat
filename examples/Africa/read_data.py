#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:28:50 2022

@author: dboateng

Reads all the model output for the late cenozoic (LGM, MH, PLIO, PI)
Present-day and future scenarios would be updated for further analysis (second paper for the warm period only 
                                                                        (Eocene, Miocene, MH, PLIO and rcp45-85 of cmip6 data))
"""
#import modules
import os 
from pyClimat.data import read_from_path, read_ECHAM_processed

# define paths 
main_path = "D:/Datasets/Model_output_pst/"
lgm_path = os.path.join(main_path, "LGM", "MONTHLY_MEANS")
plio_path = os.path.join(main_path, "PLIO", "MONTHLY_MEANS")
mh_path = os.path.join(main_path, "MH", "MONTHLY_MEANS")
pi_path = os.path.join(main_path, "PI", "MONTHLY_MEANS")
pd_path = os.path.join(main_path, )


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



exp_name = "t004_dkrz-mistral_e5w2.3_AMIP_t159l31.6h"    # simulation with present-day simulation (not different from PI simulations)

years = "1980_2000"
period = "1m"

PD_data = read_ECHAM_processed(main_path=main_path , exp_name= exp_name, years=years,
                                                  period=period, read_wiso=False)