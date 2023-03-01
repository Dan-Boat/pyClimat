# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:54:06 2023

@author: dboateng

Contain all the variables that lead to the directory of the used data
"""

import os

#defining paths
main_path_to_data = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots/save_files/PD"
main_path_to_obs = "C:/Users/dboateng/Desktop/Datasets/NAO"
path_to_plots = "C:/Users/dboateng/Desktop/Python_scripts/ClimatPackage_repogit/examples/Europe/plots"


EA_valencia_DJF_path = os.path.join(main_path_to_obs, "EA_Valencia_DJF.csv")
EA_valencia_JJA_path = os.path.join(main_path_to_obs, "EA_Valencia_JJA.csv")

NAO_Gilbraltar_JJA_path = os.path.join(main_path_to_obs, "NAO_Gilbraltar_JJA.csv")
NAO_Gilbraltar_DJF_path = os.path.join(main_path_to_obs, "NAO_Gilbraltar_DJF.csv")

NAO_CDC_JJA_path = os.path.join(main_path_to_obs, "NAO_CDC_JJA.csv")
NAO_CDC_DJF_path = os.path.join(main_path_to_obs, "NAO_CDC_DJF.csv")

EOFs_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_eofsAsCovariance.nc") # standard method
EOFs_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_eofsAsCovariance.nc")

PCs_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_pcs.csv") # standard method
PCs_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_pcs.csv")

Vars_ERA_JJA_path = os.path.join(main_path_to_data, "ERA5_standard_eof_JJA_pcs_variance.csv") # standard method
Vars_ERA_DJF_path = os.path.join(main_path_to_data, "ERA5_standard_eof_DJF_pcs_variance.csv")


EOFs_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_eofsAsCovariance.nc") # standard method
EOFs_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_eofsAsCovariance.nc")

PCs_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_pcs.csv") # standard method
PCs_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_pcs.csv")

Vars_ECHAM_JJA_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_JJA_pcs_variance.csv") # standard method
Vars_ECHAM_DJF_path = os.path.join(main_path_to_data, "ECHAM5-wiso_standard_eof_pd_DJF_pcs_variance.csv")
