# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:19:03 2023

@author: dboateng
"""

import os 

#defining paths
main_path_to_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/LGM"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/LGM"


#path to eofs 

AWI_DJF_EOFs = os.path.join(main_path_to_data, "AWI-ESM-1-1-LR_standard_eof_lgmDJF_eofsAsCovariance.nc")
AWI_DJF_VARS = os.path.join(main_path_to_data, "AWI-ESM-1-1-LR_standard_eof_lgmDJF_pcs_variance.csv")
AWI_DJF_PCS = os.path.join(main_path_to_data, "AWI-ESM-1-1-LR_standard_eof_lgmDJF_pcs.csv")

CESM_WA_DJF_EOFs = os.path.join(main_path_to_data, "CESM2-WACCM-FV2_standard_eof_lgmDJF_eofsAsCovariance.nc")
CESM_WA_DJF_VARS = os.path.join(main_path_to_data, "CESM2-WACCM-FV2_standard_eof_lgmDJF_pcs_variance.csv")
CESM_WA_DJF_PCS = os.path.join(main_path_to_data, "CESM2-WACCM-FV2_standard_eof_lgmDJF_pcs.csv")


INM_DJF_EOFs = os.path.join(main_path_to_data, "INM-CM4-8_standard_eof_lgmDJF_eofsAsCovariance.nc")
INM_DJF_VARS = os.path.join(main_path_to_data, "INM-CM4-8_standard_eof_lgmDJF_pcs_variance.csv")
INM_DJF_PCS = os.path.join(main_path_to_data, "INM-CM4-8_standard_eof_lgmDJF_pcs.csv")


MIROC_DJF_EOFs = os.path.join(main_path_to_data, "MIROC-ES2L_standard_eof_lgmDJF_eofsAsCovariance.nc")
MIROC_DJF_VARS = os.path.join(main_path_to_data, "MIROC-ES2L_standard_eof_lgmDJF_pcs_variance.csv")
MIROC_DJF_PCS = os.path.join(main_path_to_data, "MIROC-ES2L_standard_eof_lgmDJF_pcs.csv")

MPI_ESM_DJF_EOFs = os.path.join(main_path_to_data, "MPI-ESM1-2-LR_standard_eof_lgmDJF_eofsAsCovariance.nc")
MPI_ESM_DJF_VARS = os.path.join(main_path_to_data, "MPI-ESM1-2-LR_standard_eof_lgmDJF_pcs_variance.csv")
MPI_ESM_DJF_PCS = os.path.join(main_path_to_data, "MPI-ESM1-2-LR_standard_eof_lgmDJF_pcs.csv")


