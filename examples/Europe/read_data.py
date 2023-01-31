# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:44:35 2023

@author: dboateng

# read all the require dataset for the NAO-d18Op analysis and set thier paths
"""
import os 
import pandas as pd
from pyClimat.data import read_from_path, read_ERA_processed

main_path = "D:/Datasets/Model_output_pst/"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"



lgm_path = os.path.join(main_path, "LGM")
plio_path = os.path.join(main_path, "PLIO")
mh_path = os.path.join(main_path, "MH")
pi_path = os.path.join(main_path, "PI")
pd_path = os.path.join(main_path, "PD")



# read data
PD_data = read_from_path(pd_path, "PD_1980_2014_monthly.nc", decode=True)
PI_data = read_from_path(pi_path, "PI_1003_1017_monthly.nc", decode=True)
LGM_data = read_from_path(lgm_path,"LGM_1003_1017_monthly.nc", decode=True)
PLIO_data = read_from_path(plio_path, "PLIO_1003_1017_monthly.nc", decode=True)
MH_data = read_from_path(mh_path, "MH_1003_1017_monthly.nc", decode=True)

# read wiso

def read_wiso():
    PD_wiso = read_from_path(pd_path, "PD_1980_2014_monthly_wiso.nc", decode=True)
    PI_wiso = read_from_path(pi_path, "PI_1003_1017_monthly_wiso.nc", decode=True)
    LGM_wiso = read_from_path(lgm_path,"LGM_1003_1017_monthly_wiso.nc", decode=True)
    PLIO_wiso = read_from_path(plio_path, "PLIO_1003_1017_monthly_wiso.nc", decode=True)
    MH_wiso = read_from_path(mh_path, "MH_1003_1017_monthly_wiso.nc", decode=True)

#read ERA5_datasets (slp, tp, geopoth, relative humidity..)
ERA5_path = "C:/Users/dboateng/Desktop/Datasets/ERA5/monthly_1950_2021/"

ERA5_tp_path = os.path.join(ERA5_path, "tp_monthly.nc")
ERA5_t2m_path = os.path.join(ERA5_path, "t2m_monthly.nc")
ERA5_msl_path = os.path.join(ERA5_path, "msl_monthly.nc")

#define time of amip 
from1980to2014 = pd.date_range(start="1980-01-01", end="2014-12-31", freq="MS")


#read in postprocessed and analysed data

#def read_ERA(): 
ERA5_t2m = read_ERA_processed(path=ERA5_t2m_path, varname="t2m")   - 273.15 #°C
ERA5_tp = read_ERA_processed(path=ERA5_tp_path, varname="tp") * 1000 * 30  #mm/month
ERA5_msl = read_ERA_processed(path=ERA5_msl_path, varname="msl") / 100 #Pa --> hPa


# load GNIP datasets 




