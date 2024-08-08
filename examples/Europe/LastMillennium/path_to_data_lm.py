# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:57:48 2023

@author: dboateng
"""
import os

main_path_to_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/Last_Millennium/save_files"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/LastMillennium/plots"


filenames = ["CESM", "ECHAM5", "GISS", "HADCM3", "CCSM"]
#path to covariance
CESM_DJF_EOFs = os.path.join(main_path_to_data, "CESM_DJF_eofsAsCovariance.nc")
CCSM_DJF_EOFs = os.path.join(main_path_to_data, "CCSM_DJF_eofsAsCovariance.nc")
GISS_DJF_EOFs = os.path.join(main_path_to_data, "GISS_DJF_eofsAsCovariance.nc")
ECHAM_DJF_EOFs = os.path.join(main_path_to_data, "ECHAM5_DJF_eofsAsCovariance.nc")
HADESM_DJF_EOFs = os.path.join(main_path_to_data, "HADCM3_DJF_eofsAsCovariance.nc")


# PCS
CESM_DJF_PCs = os.path.join(main_path_to_data, "CESM_DJF_pcs.csv")
CCSM_DJF_PCs = os.path.join(main_path_to_data, "CCSM_DJF_pcs.csv")
ECHAM_DJF_PCs = os.path.join(main_path_to_data, "ECHAM5_DJF_pcs.csv")
GISS_DJF_PCs = os.path.join(main_path_to_data, "GISS_DJF_pcs.csv")
HADESM_DJF_PCs = os.path.join(main_path_to_data, "HADCM3_DJF_pcs.csv")


# PCS
CESM_DJF_VARs = os.path.join(main_path_to_data, "CESM_DJF_pcs_variance.csv")
CCSM_DJF_VARs = os.path.join(main_path_to_data, "CCSM_DJF_pcs_variance.csv")
ECHAM_DJF_VARs = os.path.join(main_path_to_data, "ECHAM5_DJF_pcs_variance.csv")
GISS_DJF_VARs = os.path.join(main_path_to_data, "GISS_DJF_pcs_variance.csv")
HADESM_DJF_VARs = os.path.join(main_path_to_data, "HADCM3_DJF_pcs_variance.csv")
