#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 17:04:01 2022

@author: dboateng
"""

# defining paths 

import sys
sys.path.append("/home/dboateng/Python_scripts/ClimatPackage_repogit") 

from Package import *


def extract_and_save(main_path, exp_name, path_to_save, filename_surf, filename_lev, 
                     years= "1003_1017", period = "1m"):
    
    data, data_wiso = read_ECHAM_processed(main_path=main_path , exp_name= exp_name, years=years,
                                                      period=period)
    
    data_omega = read_ECHAM_processed(main_path=main_path , exp_name= exp_name, years=years,
                                                      period=period, add_name="omega", read_wiso=False)
    
    data_plev = read_ECHAM_processed(main_path=main_path , exp_name= exp_name, years=years,
                                                      period=period, add_name="plev", read_wiso=False)
    
    
    temp = extract_var(Dataset=data , varname="temp2", units="°C")
    prec = extract_var(Dataset= data , varname="prec", units="mm/month")
    d18Op = extract_var(Dataset=data , varname="d18op", units="per mil", Dataset_wiso= data_wiso)
    u10 = extract_var(Dataset=data , varname="u10")
    v10 = extract_var(Dataset=data , varname="v10")
    elev = extract_var(Dataset=data , varname="elev", units="m")
    
    omega = extract_var(Dataset=data_omega, varname="omega", units="Pa/s", lev_units="hPa")
    cloud = extract_var(Dataset=data_plev, varname="aclcac", lev_units="hPa")
    
    
    surf_varnames = ["temp", "prec", "d18Op", "u10", "v10", "elev"]
    surf_vardata = [temp, prec, d18Op, u10, v10, elev]
    
    
    plev_varnames = ["omega", "cloud"]
    plev_vardata = [omega, cloud]
    
    data_surf = xr.Dataset()
    data_lev = xr.Dataset()
    
    
    for i in range(len(surf_varnames)):
        data_surf[surf_varnames[i]] = surf_vardata[i]
        
    for i in range(len(plev_varnames)):
        data_lev[plev_varnames[i]] = plev_vardata[i]
        
        
    data_surf.to_netcdf(os.path.join(path_to_save, filename_surf))
    data_lev.to_netcdf(os.path.join(path_to_save, filename_lev))
        

    
path_to_save = "/home/dboateng/Model_output_pst/zenodo"

# Path to experiments
main_path = "/home/dboateng/Model_output_pst"
exp_name_CTL = "a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h"
exp_name_W1E0 = "a003_hpc-bw_e5w2.3_t159_PI_Alps_east_0_t159l31.6h"
exp_name_W1E2= "a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h"

# west set up 
exp_name_W2E1 = "a009_hpc-bw_e5w2.3_t159_PI_AW200E100_t159l31.6h"
exp_name_W2E0 = "a010_hpc-bw_e5w2.3_t159_PI_AW200E0_t159l31.6h"
exp_name_W2E2 = "t017_hpc-bw_e5w2.3_PI_Alps150_t159l31.6h"
exp_name_W1E1_5 = "a008_hpc-bw_e5w2.3_t159_PI_Alps_east_150_t159l31.6h"
exp_name_W0E0 = "t015_dkrz-mistral_e5w2.3_PI_Alps0_t159l31.6h"


# extract_and_save(main_path=main_path, exp_name=exp_name_CTL, path_to_save=path_to_save, 
#                  filename_surf="CTL_surf.nc", filename_lev="CTL_plev.nc")

# extract_and_save(main_path=main_path, exp_name=exp_name_W1E0, path_to_save=path_to_save, 
#                  filename_surf="W1E0_surf.nc", filename_lev="W1E0_plev.nc")

# extract_and_save(main_path=main_path, exp_name=exp_name_W1E2, path_to_save=path_to_save, 
#                  filename_surf="W1E2_surf.nc", filename_lev="W1E2_plev.nc")

# extract_and_save(main_path=main_path, exp_name=exp_name_W2E1, path_to_save=path_to_save, 
#                  filename_surf="W2E1_surf.nc", filename_lev="W2E1_plev.nc")

# extract_and_save(main_path=main_path, exp_name=exp_name_W2E0, path_to_save=path_to_save, 
#                  filename_surf="W2E0_surf.nc", filename_lev="W2E0_plev.nc")

# extract_and_save(main_path=main_path, exp_name=exp_name_W2E2, path_to_save=path_to_save, 
#                  filename_surf="W2E2_surf.nc", filename_lev="W2E2_plev.nc")

# extract_and_save(main_path=main_path, exp_name=exp_name_W1E1_5, path_to_save=path_to_save, 
#                  filename_surf="W1E1_5_surf.nc", filename_lev="W1E1_5_plev.nc")

extract_and_save(main_path=main_path, exp_name=exp_name_W0E0, path_to_save=path_to_save, 
                 filename_surf="W0E0_surf.nc", filename_lev="W0E0_plev.nc")
