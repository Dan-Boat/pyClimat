#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 16:55:45 2021

@author: Daniel Boateng

Reading data routine for Climat (required user declarations of paths to datasets eg. Reanalysis
                                 , ECHAM, DWD stations and Gtopo files etc)
                                 
Note: All User specifications must be declared in the control script which will import all the functions defined here                                
"""
# Importing modules

import xarray as xr
import os
import pandas as pd
import numpy as np


def read_ECHAM_processed(main_path, exp_name, years="1003_1017", period="1m", add_name=None, read_wiso=True):
    """
    Reads output processed from ECHAM

    Parameters
    ----------
    main_path : TYPE: STR
        DESCRIPTION. directory to the main path for all module output (eg. esd02-->ESD, or local path)
    exp_name : TYPE : STR
        DESCRIPTION. Name of experiment output (eg. a003_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h)
    years : TYPE, optional
        DESCRIPTION. The default is "1003_1017". or range of year you have processed
    period : TYPE, optional
        DESCRIPTION. The default is "1m". or 1d, 1y if implemted here!
    add_name : TYPE, optional
        DESCRIPTION. The default is None. or eg. _msl (for particular variable)

    Returns
    -------
    data : TYPE: Dataset
        DESCRIPTION. Dataset of echam ouput will some or all variables 
    data_wiso : TYPE: Dataset
        DESCRIPTION. Dataset of echam wiso ouput will some or all variables 

    """
    
    module_path = os.path.join(main_path, exp_name)
    processed_path = os.path.join(module_path, "output_processed")
    if period == "1m":
        
        processed_mmeans = os.path.join(processed_path, "MONTHLY_MEANS")
        processed_mwisomeans = os.path.join(processed_path, "MONTHLY_MEANS_WISO")
        if add_name is None:
            data_path = os.path.join(processed_mmeans, years + "_" + period + "_mlterm.nc")
            data_wiso_path = os.path.join(processed_mwisomeans, years + "_wiso_" + period + "_mlterm.nc")
        else:
            data_path = os.path.join(processed_mmeans, years + "_" + period + "_" + add_name +"_mlterm.nc")
            data_wiso_path = os.path.join(processed_mwisomeans, years + "_wiso_" + period + "_" + add_name + "_mlterm.nc")
    elif period == "1a":
        print("Other periods are yet to be implemented or reading annual long-term means")
        if add_name is None:
            data_path = os.path.join(processed_path, years + "_lterm.nc")
            data_wiso_path = os.path.join(processed_path, years + "_wiso_" + "lterm.nc")
        else:
            data_path = os.path.join(processed_path, years + "_" + add_name +"_lterm.nc")
            data_wiso_path = os.path.join(processed_path, years + "_wiso_" + add_name + "_lterm.nc")
    else:
        raise ValueError("define the period of the processed output")
        
        
    
    data = xr.open_dataset(data_path, decode_cf=True, use_cftime=True)
    if read_wiso==False:
        return data
    else:
        
         data_wiso = xr.open_dataset(data_wiso_path, decode_cf=True, use_cftime=True)
    
         return data, data_wiso

def read_ECHAM_input(main_path, exp_name, filename, read_var=False, varname=None):
    """
    

    Parameters
    ----------
    main_path : TYPE: str
        DESCRIPTION. Path containing all module outputs 
    exp_name : TYPE:str
        DESCRIPTION. Name of the experiment 
    varname : TYPE: str
        DESCRIPTION. Name of the input file (jan_surf file)

    Returns
    -------
    data : TYPE: dataset
        DESCRIPTION.

    """
    path_to_exp = os.path.join(main_path, exp_name)
    path_data = os.path.join(path_to_exp, "input", filename)
    
    data = xr.open_dataset(path_data, decode_cf=True, use_cftime=True)
    if read_var==True:
        if varname == "elev":
            data = data["GEOSP"] / 9.8  # convert to m**2/s**2 --> m
        else:
            data = data[varname]
    
    return data
    
def read_ERA_processed(path, varname):
    """
    

    Parameters
    ----------
    path : TYPE: STR
        DESCRIPTION. path of the ERA dataset
    varname : TYPE:STR
        DESCRIPTION. Variable name for ERA (eg. t2m for temperature, tp:precipitation)

    Returns
    -------
    data : TYPE: datarray
        DESCRIPTION.

    """
    
    dataset = xr.open_dataset(path)
    data = dataset[varname]
    
    return data
    


def read_Gtopo(path, tile_name, extract_var=None):
    """
    

    Parameters
    ----------
    path : TYPE: str
        DESCRIPTION. directrory to all the tile files (or path to tiles)
    tile_name : TYPE: str
        DESCRIPTION. Which tile to use for modification (check the image in the original files folder)
    extract_var : TYPE, optional : or yes
        DESCRIPTION. The default is None or To extract only the values to datarray

    Returns
    -------
    TYPE: Dataset or dataarray
        DESCRIPTION. It reads a particular tile file

    """
    dataset = xr.open_dataset(os.path.join(path, tile_name))
    
    if extract_var is not None:
        if extract_var in ["yes", "Yes", "YES"]:
            data = dataset["z"]
            return data
        else:
            None
    else:
        return dataset
    
def read_GNIP_data(path, filename):
    """
    

    Parameters
    ----------
    path : TYPE: str
        DESCRIPTION. The directory holding all the data
    filename : TYPE: str
        DESCRIPTION. The name of the file

    Returns
    -------
    df : TYPE: DataFrame
        DESCRIPTION. Data containing lat, lon and d18op

    """
    
    data = np.loadtxt(os.path.join(path, filename), dtype=float, unpack=True)
    
    df = pd.DataFrame(columns=["lat", "lon", "d18op"])
    df["lat"] = data[0]
    df ["d18op"] = data[2]
    df["lon"] = data[1]
    
    return df


def read_DWD_processed(*args):
    pass
